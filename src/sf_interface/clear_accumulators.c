#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include "sf_interface_private.h"

void
clear_accumulators_pipeline( accumulators_pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline ) {
  accumulator_t * ALIGNED(16) a = args->a;
  int n = args->n, n_array = args->n_array, s_array = args->s_array, i;
  DISTRIBUTE(n, accumulators_n_block, pipeline_rank, n_pipeline, i, n); a += i;
  for( ; n_array; n_array--, a+=s_array ) CLEAR( a, n );
}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

// SPU pipeline is defined in a different compilation unit

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "V4 version not hooked up yet!"

#endif

#define VOX(x,y,z) VOXEL(x,y,z, g->nx,g->ny,g->nz)

void
clear_accumulators( accumulator_t * ALIGNED(128) a,
                    const grid_t  *              g ) {
  DECLARE_ALIGNED_ARRAY( accumulators_pipeline_args_t, 128, args, 1 );
  int i0, na;

  if( a==NULL ) ERROR(("Invalid accumulator"));
  if( g==NULL ) ERROR(("Invalid grid"));

  i0 = (VOX(1,1,1)/2)*2; /* Round i0 down to even for 128B align on Cell */

  /**/                       na = serial.n_pipeline;
  if( na<thread.n_pipeline ) na = thread.n_pipeline;
# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  if( na<spu.n_pipeline    ) na = spu.n_pipeline;
# endif
  na++; /* 1 + max( {serial,thread,spu}.n_pipeline ) */

  args->a       = a + i0;
  args->n       = ((( VOX(g->nx,g->ny,g->nz) - i0 + 1 )+1)/2)*2;
  args->n_array = na;
  args->s_array = POW2_CEIL((g->nx+2)*(g->ny+2)*(g->nz+2),2);

  EXEC_PIPELINES( clear_accumulators, args, 0 );
  WAIT_PIPELINES();
}

