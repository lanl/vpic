#include <field_pipelines.h>

void
clear_accumulators( accumulator_t * ALIGNED a,
                    const grid_t * g ) {
  clear_accumulators_pipeline_args_t args[1];
  int n_voxel;

  if( a==NULL ) ERROR(("Invalid accumulator"));
  if( g==NULL ) ERROR(("Invalid grid"));

  n_voxel = (g->nx+2)*(g->ny+2)*(g->nz+2);

  /* Have the host clear its personal accumulator while the pipelines
     are clearing theirs */
  args->a       = a;
  args->n_voxel = n_voxel;
  dispatch_pipelines( clear_accumulators_pipeline_v4, args, 0 );
  memset( a, 0, n_voxel*sizeof(accumulator_t) );
  wait_for_pipelines();
}

