#include <field_pipelines.h>

void
reduce_accumulators( accumulator_t * ALIGNED a,
                     const grid_t  *         g ) {
  reduce_accumulators_pipeline_args_t args[1];
  pipeline_request_t request[1];
  
  if( a==NULL ) { ERROR(("Bad accumulator")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));        return; }
  
  args->a = a;
  args->g = g;
  dispatch_pipelines( reduce_accumulators_pipeline, args, 0, request );
  wait_for_pipelines( request );
}

