// Any compiler worth its salt will have an optimized memset
#define IN_field_pipeline
#include <field_pipelines.h>

static void
clear_accumulators_pipeline( clear_accumulators_pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline ) {

  double target;
  accumulator_t * a;
  int n;

  // Rotate the pipeline rank so that the host thread clears
  //   accumulator(:,:,:,0)
  // and the bulk processing threads clear pipelines
  //   accumulator(:,:,:,1:n_pipeline)
  // in the typical case.  This improves cache prefetching when
  // running multithreaded as the host thread reduces to
  // accumulator 0 and so forth in the particle module. 

  pipeline_rank++; if( pipeline_rank>n_pipeline ) pipeline_rank = 0;

  // Clear a fraction of the voxels
  // FIXME: Should the host clear an equal fraction of the accumulators
  // or a clean up fraction??

  target = (double)args->n_voxel / (double)(1+n_pipeline);
  n = (int)( target*(double) pipeline_rank    + 0.5 );
  a = args->a + n;
  n = (int)( target*(double)(pipeline_rank+1) + 0.5 ) - n;

  memset( a, 0, n*sizeof(accumulator_t) );
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
  int n;

  if( a==NULL ) ERROR(("Invalid accumulator"));
  if( g==NULL ) ERROR(("Invalid grid"));

  /**/                      n = serial.n_pipeline;
  if( n<thread.n_pipeline ) n = thread.n_pipeline;
# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  if( n<spu.n_pipeline    ) n = spu.n_pipeline;
# endif
  n++; /* n = 1 + max( {serial,thread,spu}.n_pipeline ) */
  n *= POW2_CEIL((g->nx+2)*(g->ny+2)*(g->nz+2),2);

  args->a       = a;
  args->n_voxel = n;

  EXEC_PIPELINES( clear_accumulators, args, 0 );
  WAIT_PIPELINES();
}

