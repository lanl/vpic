// Any compiler worth its salt will have an optimized memset
#define IN_field_pipeline
#include <field_pipelines.h>

// FIXME: THE CLEAR FUNCTION MUST TAKE INTO ACCOUNT THAT IT MAKE NEED TO
// RUN HERE ON A DIFFERENT NUMBER OF THREADS THAN THERE ARE ACCUMULATORS!

static void
clear_accumulators_pipeline( clear_accumulators_pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline ) {
  // Have each pipeline clear its personal accumulator
  memset( args->a + pipeline_rank*args->n_voxel, 0,
          args->n_voxel*sizeof(accumulator_t) );
}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

#error "SPU version not hooked up yet!"

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "V4 version not hooked up yet!"

#endif

void
clear_accumulators( accumulator_t * ALIGNED(128) a,
                    const grid_t  *              g ) {
  clear_accumulators_pipeline_args_t args[1];

  if( a==NULL ) ERROR(("Invalid accumulator"));
  if( g==NULL ) ERROR(("Invalid grid"));

  args->a       = a;
  args->n_voxel = (g->nx+2)*(g->ny+2)*(g->nz+2);

  EXEC_PIPELINES( clear_accumulators, args, 0 );
  WAIT_PIPELINES();
}
