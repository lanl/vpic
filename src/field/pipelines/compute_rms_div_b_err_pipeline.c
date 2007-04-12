#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
compute_rms_div_b_err_pipeline( compute_rms_div_b_err_pipeline_args_t * args,
                                int pipeline_rank,
                                int n_pipeline ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;
  int n_voxel;
                             
  double err;
  field_t *f0;
  int x, y, z, nx, ny, nz;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;

# if 0 /* Original non-pipelined version */
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

  /* Process voxels assigned to this pipeline */

  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
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
