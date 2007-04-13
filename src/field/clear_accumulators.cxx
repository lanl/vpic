#include <field.h>

/* Any compiler worth its salt has highly optimized version for the
   target platform */

#undef V4_ACCELERATION

#ifndef V4_ACCELERATION
#define CLEAR_ACCUMULATORS_PIPELINE (pipeline_func_t)clear_accumulators_pipeline
#else
#define CLEAR_ACCUMULATORS_PIPELINE (pipeline_func_t)clear_accumulators_pipeline_v4
#endif

typedef struct clear_accumulators_pipeline_args {
  accumulator_t * ALIGNED a; /* Base of all the accumulators */
  int n_voxel;
} clear_accumulators_pipeline_args_t;

static void
clear_accumulators_pipeline( clear_accumulators_pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline ) {
  /* Have each pipeline clear its personal accumulator */
  memset( args->a + (pipeline_rank+1)*args->n_voxel, 0, args->n_voxel*sizeof(accumulator_t) );
}

#ifdef V4_ACCELERATION
#error "V4 version not implemented"
#endif

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
  PMETHOD.dispatch( CLEAR_ACCUMULATORS_PIPELINE, args, 0 );
  memset( a, 0, n_voxel*sizeof(accumulator_t) );
  PMETHOD.wait();
}
