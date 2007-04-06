#if 0 // Original non-pipelined non-vectorized version 
  err = 0;
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      for( x=2; x<=nx; x++ ) {
        err += (double)f0->div_e_err*(double)f0->div_e_err;
        f0++;
      }
    }
  }
#endif

#include <field_pipelines.h>
#include <v4.h>

using namespace v4;

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
compute_rms_div_e_err_pipeline( compute_rms_div_e_err_pipeline_args_t * args,
                                int pipeline_rank ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;

  field_t *f0;
  int x, y, z, n_voxel;
  
  double err;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  // Process voxels assigned to this pipeline 

  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  f0 = &f(x,y,z);

  err = 0;
  for( ; n_voxel; n_voxel-- ) {

    // FIXME: I stared at this for a long time and came to the
    // conclusion that this is not worth vectorizing even on the
    // SPUs.  Namely, this operation is not executed very often,
    // uses mixed precision arithmetic (making it V4 unfriendly),
    // does not have a V4 friendly horizontal SIMD implementation
    // and the best way to do a vertical SIMD variant would have
    // slightly different round-off properties of the scalar and
    // would require vector scalar gathing. 

    err += (double)f0->div_e_err*(double)f0->div_e_err;
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

