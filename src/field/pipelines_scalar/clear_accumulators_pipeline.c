#include <field_pipelines.h>
#include <string.h> /* For memset */

void
clear_accumulators_pipeline( clear_accumulators_pipeline_args_t * args,
                             int pipeline_rank ) {
  /* Have each pipeline clear its personal accumulator */
  memset( args->a + (pipeline_rank+1)*args->n_voxel, 0, args->n_voxel*sizeof(accumulator_t) );
}

