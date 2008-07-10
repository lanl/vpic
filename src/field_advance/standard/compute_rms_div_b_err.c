#define IN_sfa
#include "sfa_private.h"

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

typedef struct pipeline_args {
  field_t      * ALIGNED(128) f;
  const grid_t *              g;
  double err[MAX_PIPELINE+1];
} pipeline_args_t;

static void
pipeline( pipeline_args_t * args,
          int pipeline_rank,
          int n_pipeline ) {
  field_t      * ALIGNED(128) f = args->f;
  const grid_t *              g = args->g;
                             
  field_t * ALIGNED(16) f0;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  double err;

  // Process voxels assigned to this pipeline

  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );
  
  f0 = &f(x,y,z);

  err = 0;
  for( ; n_voxel; n_voxel-- ) {
    err += f0->div_b_err*f0->div_b_err;
    f0++;
    
    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      f0 = &f(x,y,z);
    }
  }
    
  args->err[pipeline_rank] = err;
}

double
compute_rms_div_b_err( field_t      * ALIGNED(128) f,
                       const grid_t *              g ) {
  pipeline_args_t args[1];
  int p;
  
  double err = 0, local[2], global[2];

  if( f==NULL ) ERROR(("Bad field"));
  if( g==NULL ) ERROR(("Bad grid"));

# if 0 // Original non-pipelined version
  err = 0;
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      for( x=1; x<=nx; x++ ) {
        err += f0->div_b_err*f0->div_b_err;
        f0++;
      }
    }
  }
# endif

  args->f = f;
  args->g = g;

  EXEC_PIPELINES( pipeline, args, 0 );
  WAIT_PIPELINES();

  err = 0;
  for( p=0; p<=N_PIPELINE; p++ ) err += args->err[p];
  local[0] = err*g->dx*g->dy*g->dz;
  local[1] = g->nx*g->ny*g->nz*g->dx*g->dy*g->dz;
  mp_allsum_d( local, global, 2, g->mp );
  return g->eps0*sqrt(global[0]/global[1]);
}
