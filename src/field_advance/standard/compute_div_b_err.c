#define IN_sfa
#include "sfa_private.h"

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

typedef struct pipeline_args {
  field_t      * ALIGNED(128) f;
  const grid_t *              g;
} pipeline_args_t;

static void
pipeline( pipeline_args_t * args,
          int pipeline_rank,
          int n_pipeline ) {
  field_t      * ALIGNED(128) f = args->f;
  const grid_t *              g = args->g;
  
  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? g->rdx : 0;
  const float py = (ny>1) ? g->rdy : 0;
  const float pz = (nz>1) ? g->rdz : 0;

  // Process the voxels assigned to this pipeline
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x+1,y,  z  ); \
  fy = &f(x,  y+1,z  ); \
  fz = &f(x,  y,  z+1)

  LOAD_STENCIL();

  for( ; n_voxel; n_voxel-- ) {
    f0->div_b_err = px*( fx->cbx - f0->cbx ) +
                    py*( fy->cby - f0->cby ) +
                    pz*( fz->cbz - f0->cbz );
    f0++; fx++; fy++; fz++;
    
    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      LOAD_STENCIL();
    }
  }

# undef LOAD_STENCIL

}

void
compute_div_b_err( field_t      * ALIGNED(128) f,
                   const grid_t *              g ) {
  pipeline_args_t args[1];
  
  if( f==NULL ) ERROR(("Bad field"));
  if( g==NULL ) ERROR(("Bad grid"));
  
# if 0 // Original non-pipelined version
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(2,y,z);
      fy = &f(1,y+1,z);
      fz = &f(1,y,z+1);
      for( x=1; x<=nx; x++ ) {
	f0->div_b_err = px*( fx->cbx - f0->cbx ) +
	                py*( fy->cby - f0->cby ) +
                        pz*( fz->cbz - f0->cbz );
	f0++; fx++; fy++; fz++;
      }
    }
  }
# endif

  args->f = f;
  args->g = g;

  EXEC_PIPELINES( pipeline, args, 0 );
  WAIT_PIPELINES();
}
