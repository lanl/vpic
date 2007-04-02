#include <field_pipelines.h>
#include CONCAT3(<,V4VERSION,>)

using namespace v4;

#define a(x,y,z) a[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
reduce_accumulators_pipeline( reduce_accumulators_pipeline_args_t * args,
                              int pipeline_rank ) {
  accumulator_t * ALIGNED a = args->a;
  const grid_t  *         g = args->g;

  int p, x, y, z, n_voxel;
  float * ALIGNED pa, * ALIGNED paa;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;
  const int stride = 12*(nx+2)*(ny+2)*(nz+2);

  v4float vjxx, vjyy, vjzz;
  v4float vjx,  vjy,  vjz;
  
  // Process voxels assigned to this pipeline

  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );


  pa = &a(x,y,z).jx[0];
  
  for( ; n_voxel; n_voxel-- ) {
  
    // More evil than usual---Duff inspired V4 accelerated loop.
    // As n_pipeline is generally compile time determined, this will
    // reduce down to the desired branchless assembly at compile
    // compile.

    load_4x1( pa,   vjx );
    load_4x1( pa+4, vjy );
    load_4x1( pa+8, vjz );
    paa = pa + stride;
    p   = n_pipeline;

#   define REDUCE_ACCUM()                 \
    load_4x1( paa,   vjxx ); vjx += vjxx; \
    load_4x1( paa+4, vjyy ); vjy += vjyy; \
    load_4x1( paa+8, vjzz ); vjz += vjzz; \
    paa += stride

    switch( p ) {
    default: for(;p>8;p--) { REDUCE_ACCUM(); }
    case 8:  REDUCE_ACCUM();
    case 7:  REDUCE_ACCUM();
    case 6:  REDUCE_ACCUM();
    case 5:  REDUCE_ACCUM();
    case 4:  REDUCE_ACCUM();
    case 3:  REDUCE_ACCUM();
    case 2:  REDUCE_ACCUM();
    case 1:  REDUCE_ACCUM();
    case 0:  break;
    }

    store_4x1( vjx, pa   );
    store_4x1( vjy, pa+4 );
    store_4x1( vjz, pa+8 );

    pa += 12;

    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      pa = &a(x,y,z).jx[0];
    }
  }
}

