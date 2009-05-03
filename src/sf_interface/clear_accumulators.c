#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include "sf_interface_private.h"

void
clear_accumulators_pipeline( accumulators_pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline ) {
  int i, n;
  DISTRIBUTE( args->n_array*args->stride, accumulators_n_block,
              pipeline_rank, n_pipeline, i, n );
  CLEAR( args->a+i, n );
}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

// SPU pipeline is defined in a different compilation unit

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "V4 version not hooked up yet!"

#endif

void
clear_accumulators( accumulator_t * ALIGNED(128) a,
                    const grid_t  *              g ) {
  accumulators_pipeline_args_t args[1];
  int n_array, stride;

  if( a==NULL ) ERROR(("Invalid accumulator"));
  if( g==NULL ) ERROR(("Invalid grid"));

  /**/                            n_array = serial.n_pipeline;
  if( n_array<thread.n_pipeline ) n_array = thread.n_pipeline;
# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  if( n_array<spu.n_pipeline    ) n_array = spu.n_pipeline;
# endif
  n_array++; /* n_array = 1 + max( {serial,thread,spu}.n_pipeline ) */

  stride = POW2_CEIL((g->nx+2)*(g->ny+2)*(g->nz+2),2);

  if( stride<accumulators_n_block ) CLEAR( a, n_array*stride );
  else {
    args->a = a, args->n_array = n_array, args->stride = stride;
    EXEC_PIPELINES( clear_accumulators, args, 0 );
    WAIT_PIPELINES();
  }
}

