#include <field_pipelines.h>

void
clear_accumulators_pipeline( clear_accumulators_pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline ) {
  /* Have each pipeline clear its personal accumulator */
  memset( args->a + (pipeline_rank+1)*args->n_voxel, 0, args->n_voxel*sizeof(accumulator_t) );
}

