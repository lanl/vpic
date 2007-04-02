#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
compute_rms_div_e_err_pipeline( compute_rms_div_e_err_pipeline_args_t * args,
                                int pipeline_rank ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;
  int n_voxel;
  
  double err;
  int x, y, z, nx, ny, nz;
  field_t *f0;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;

# if 0 /* Original non-pipelined version */
  err = 0;
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      for( x=2; x<=nx; x++ ) {
        err += f0->div_e_err*f0->div_e_err;
        f0++;
      }
    }
  }
# endif

  /* Process voxels assigned to this pipeline */

  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  f0 = &f(x,y,z);

  err = 0;
  for( ; n_voxel; n_voxel-- ) {
    err += f0->div_e_err*f0->div_e_err;
    f0++;

    x++;
    if( x>nx ) {
      x=2, y++;
      if( y>ny ) y=2, z++;
      f0 = &f(x,y,z);
    }
  }

  args->err[pipeline_rank] = err;
}

