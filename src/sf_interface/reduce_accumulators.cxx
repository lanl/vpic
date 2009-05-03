#define IN_sf_interface
#define HAS_V4_PIPELINE
#include "sf_interface_private.h"

// FIXME: THIS NEEDS TO REDUCE THE ACTUAL NUMBER OF ACCUMULATORS PRESENT
// AS THE NUMBER OF PIPELINES EXECUTED HERE MAY NOT MATCH THE NUMBER OF
// ACCUMULATORS USED DURING OTHER OPERATIONS (E.G. 8 SPU PIPELINES AND
// 2 PPU PIPELINES!)

#define a(x,y,z) a[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
reduce_accumulators_pipeline( reduce_accumulators_pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline ) {
  accumulator_t * ALIGNED(128) a = args->a;
  const grid_t  *              g = args->g;

  const accumulator_t * ALIGNED(16) sa;
  accumulator_t       * ALIGNED(16) da;
  int n, x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;
  const int na = args->na;
  const int stride = POW2_CEIL((nx+2)*(ny+2)*(nz+2),2);

  // Process voxels assigned to this pipeline

  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  da = &a(x,y,z);

  for( ; n_voxel; n_voxel-- ) {

    sa = da + stride;
    for( n=1; n<na; n++ ) {
      da->jx[0] += sa->jx[0]; da->jx[1] += sa->jx[1]; da->jx[2] += sa->jx[2]; da->jx[3] += sa->jx[3];
      da->jy[0] += sa->jy[0]; da->jy[1] += sa->jy[1]; da->jy[2] += sa->jy[2]; da->jy[3] += sa->jy[3];
      da->jz[0] += sa->jz[0]; da->jz[1] += sa->jz[1]; da->jz[2] += sa->jz[2]; da->jz[3] += sa->jz[3];
      sa += stride;
    }

    da++;

    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      da = &a(x,y,z);
    }
  }
}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

#error "SPU version not hooked up yet!"

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

using namespace v4;

void
reduce_accumulators_pipeline_v4( reduce_accumulators_pipeline_args_t * args,
                                 int pipeline_rank,
                                 int n_pipeline ) {
  accumulator_t * ALIGNED(128) a = args->a;
  const grid_t  *              g = args->g;

  const accumulator_t * ALIGNED(16) sa;
  accumulator_t       * ALIGNED(16) da;
  int n, x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;
  const int na = args->na;
  const int stride = POW2_CEIL((nx+2)*(ny+2)*(nz+2),2);

  v4float vjxx, vjyy, vjzz;
  v4float vjx,  vjy,  vjz;

  // Process voxels assigned to this pipeline

  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  da = &a(x,y,z);

  for( ; n_voxel; n_voxel-- ) {

    // More evil than usual---Duff inspired V4 accelerated loop.

    load_4x1( da->jx, vjx );
    load_4x1( da->jy, vjy );
    load_4x1( da->jz, vjz );
    sa = da + stride;
    n = na - 1;

#   define REDUCE_ACCUM()                  \
    load_4x1( sa->jx, vjxx ); vjx += vjxx; \
    load_4x1( sa->jy, vjyy ); vjy += vjyy; \
    load_4x1( sa->jz, vjzz ); vjz += vjzz; \
    sa += stride

    switch( n ) {
    default: for(;n>8;n--) { REDUCE_ACCUM(); }
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

    store_4x1( vjx, da->jx );
    store_4x1( vjy, da->jy );
    store_4x1( vjz, da->jz );

    da++;

    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      da = &a(x,y,z);
    }
  }
}

#endif

void
reduce_accumulators( accumulator_t * ALIGNED(128) a,
                     const grid_t  *              g ) {
  reduce_accumulators_pipeline_args_t args[1];
  int na;
  
  if( a==NULL ) ERROR(("Bad accumulator"));
  if( g==NULL ) ERROR(("Bad grid"));

  /**/                       na = serial.n_pipeline;
  if( na<thread.n_pipeline ) na = thread.n_pipeline;
# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  if( na<spu.n_pipeline    ) na = spu.n_pipeline;
# endif
  na++; /* na = 1 + max( {serial,thread,spu}.n_pipeline ) */

  args->a  = a;
  args->g  = g;
  args->na = na;

  EXEC_PIPELINES( reduce_accumulators, args, 0 );
  WAIT_PIPELINES();
}

