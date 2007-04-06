#include <field_pipelines.h>

void
unload_accumulator( field_t * ALIGNED f, 
                    const accumulator_t * ALIGNED a,
                    const grid_t * g ) {
  unload_accumulator_pipeline_args_t args[1];
  
  if( f==NULL ) ERROR(("Bad field"));
  if( a==NULL ) ERROR(("Bad accumulator"));
  if( g==NULL ) ERROR(("Bad grid"));

  args->f = f;
  args->a = a;
  args->g = g;
  dispatch_pipelines( unload_accumulator_pipeline_v4, args, 0 );
  wait_for_pipelines();
}

