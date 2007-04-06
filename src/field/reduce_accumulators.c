#include <field_pipelines.h>

void
reduce_accumulators( accumulator_t * ALIGNED a,
                     const grid_t  *         g ) {
  reduce_accumulators_pipeline_args_t args[1];
  
  if( a==NULL ) ERROR(("Bad accumulator"));
  if( g==NULL ) ERROR(("Bad grid"));
  
  args->a = a;
  args->g = g;
  dispatch_pipelines( reduce_accumulators_pipeline, args, 0 );
  wait_for_pipelines();
}

