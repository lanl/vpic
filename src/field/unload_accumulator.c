#include <field_pipelines.h>

void
unload_accumulator( field_t * ALIGNED f, 
                    const accumulator_t * ALIGNED a,
                    const grid_t * g ) {
  unload_accumulator_pipeline_args_t args[1];
  pipeline_request_t request[1];
  
  if( f==NULL ) { ERROR(("Bad field"));       return; }
  if( a==NULL ) { ERROR(("Bad accumulator")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));        return; }

  args->f = f;
  args->a = a;
  args->g = g;
  dispatch_pipelines( unload_accumulator_pipeline, args, 0, request );
  wait_for_pipelines( request );
}

