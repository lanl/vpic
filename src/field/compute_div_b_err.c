#include <field_pipelines.h>

void compute_div_b_err( field_t * ALIGNED f,
                        const grid_t * g ) {
  compute_div_b_err_pipeline_args_t args[1];
  
  if( f==NULL ) ERROR(("Bad field"));
  if( g==NULL ) ERROR(("Bad grid"));
  
  args->f = f;
  args->g = g;  
  dispatch_pipelines( compute_div_b_err_pipeline, args, 0 );
  wait_for_pipelines();
}
