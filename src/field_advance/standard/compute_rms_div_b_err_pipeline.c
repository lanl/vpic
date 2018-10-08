#define IN_sfa

#include "sfa_private.h"

#include "../../util/pipelines/pipelines_exec.h"

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

typedef struct pipeline_args
{
  const field_t * ALIGNED(128) f;
  const grid_t  *              g;
  double err[MAX_PIPELINE+1];
} pipeline_args_t;

static void
compute_rms_div_b_err_pipeline_scalar( pipeline_args_t * args,
                                       int pipeline_rank,
                                       int n_pipeline )
{
  const field_t * ALIGNED(128) f = args->f;
  const grid_t  *              g = args->g;
                             
  const field_t * ALIGNED(16) f0;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  double err;

  // Process voxels assigned to this pipeline

  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );
  
  f0 = &f(x,y,z);

  err = 0;
  for( ; n_voxel; n_voxel-- )
  {
    err += f0->div_b_err*f0->div_b_err;
    f0++;

    x++;
    if ( x > nx )
    {
      x=1, y++;
      if( y>ny ) y=1, z++;
      f0 = &f(x,y,z);
    }
  }
    
  args->err[pipeline_rank] = err;
}

double
compute_rms_div_b_err_pipeline( const field_array_t * fa )
{
  pipeline_args_t args[1];
  int p;
  
  double err = 0, local[2], global[2];

  if ( !fa )
  {
    ERROR( ( "Bad args") );
  }

# if 0 // Original non-pipelined version
  field_t * ALIGNED(16) f0;
  int z, y, x;
  int nx = g->nx;
  int ny = g->ny;
  int nz = g->nz;

  err = 0;
  for( z=1; z<=nz; z++ )
  {
    for( y=1; y<=ny; y++ )
    {
      f0 = &f(1,y,z);
      for( x=1; x<=nx; x++ )
      {
        err += f0->div_b_err*f0->div_b_err;
        f0++;
      }
    }
  }
# endif

  args->f = fa->f;
  args->g = fa->g;

  EXEC_PIPELINES( compute_rms_div_b_err, args, 0 );
  WAIT_PIPELINES();

  err = 0;
  for( p=0; p<=N_PIPELINE; p++ ) err += args->err[p];

  local[0] = err*fa->g->dV;
  local[1] = (fa->g->nx*fa->g->ny*fa->g->nz)*fa->g->dV;
  mp_allsum_d( local, global, 2 );
  return fa->g->eps0*sqrt(global[0]/global[1]);
}
