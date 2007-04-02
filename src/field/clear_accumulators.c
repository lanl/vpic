#include <field_pipelines.h>
#include <string.h> /* For memset */

void
clear_accumulators( accumulator_t * ALIGNED a,
                    const grid_t * g ) {
  clear_accumulators_pipeline_args_t args[1];
  pipeline_request_t request[1];
  int n_voxel;

  if( a==NULL ) { ERROR(("Invalid accumulator")); return; }
  if( g==NULL ) { ERROR(("Invalid grid"));        return; }

  n_voxel = (g->nx+2)*(g->ny+2)*(g->nz+2);

  /* Have the host clear its personal accumulator while the pipelines
     are clearing theirs */
  args->a       = a;
  args->n_voxel = n_voxel;
  dispatch_pipelines( clear_accumulators_pipeline, args, 0, request );
  memset( a, 0, n_voxel*sizeof(accumulator_t) );
  wait_for_pipelines( request );
}

