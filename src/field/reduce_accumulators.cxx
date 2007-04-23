#include <field.h>

#ifndef V4_ACCELERATION
#define REDUCE_ACCUMULATORS_PIPELINE (pipeline_func_t)reduce_accumulators_pipeline
#else
#define REDUCE_ACCUMULATORS_PIPELINE (pipeline_func_t)reduce_accumulators_pipeline_v4
#endif

#define a(x,y,z) a[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

typedef struct reduce_accumulators_pipeline_args {
  accumulator_t * ALIGNED(128) a;
  const grid_t  *              g;
} reduce_accumulators_pipeline_args_t;

static void
reduce_accumulators_pipeline( reduce_accumulators_pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline ) {
  accumulator_t * ALIGNED(128) a = args->a;
  const grid_t  *              g = args->g;

  float * ALIGNED(16) pa, * ALIGNED(16) paa;
  int p, x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;
  const int stride = 12*(nx+2)*(ny+2)*(nz+2);

  // Process voxels assigned to this pipeline

  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  pa = &a(x,y,z).jx[0];

  for( ; n_voxel; n_voxel-- ) {

    paa = pa + stride;
    for( p=0; p<n_pipeline; p++ ) {
      pa[0]  += paa[0];
      pa[1]  += paa[1];
      pa[2]  += paa[2];
      pa[3]  += paa[3];
      pa[4]  += paa[4];
      pa[5]  += paa[5];
      pa[6]  += paa[6];
      pa[7]  += paa[7];
      pa[8]  += paa[8];
      pa[9]  += paa[9];
      pa[10] += paa[10];
      pa[11] += paa[11];
      paa    += stride;
    }
    pa += 12;

    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      pa = &a(x,y,z).jx[0];
    }
  }
}

#ifdef V4_ACCELERATION

using namespace v4;

static void
reduce_accumulators_pipeline_v4( reduce_accumulators_pipeline_args_t * args,
                                 int pipeline_rank,
                                 int n_pipeline ) {
  accumulator_t * ALIGNED(128) a = args->a;
  const grid_t  *              g = args->g;

  float * ALIGNED(16) pa, * ALIGNED(16) paa;
  int p, x, y, z, n_voxel;

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

#   undef REDUCE_ACCUM

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

#endif

void
reduce_accumulators( accumulator_t * ALIGNED(128) a,
                     const grid_t  *              g ) {
  reduce_accumulators_pipeline_args_t args[1];
  
  if( a==NULL ) ERROR(("Bad accumulator"));
  if( g==NULL ) ERROR(("Bad grid"));
  
  args->a = a;
  args->g = g;

  PSTYLE.dispatch( REDUCE_ACCUMULATORS_PIPELINE, args, 0 );
  reduce_accumulators_pipeline( args, PSTYLE.n_pipeline, PSTYLE.n_pipeline );
  PSTYLE.wait();
}
