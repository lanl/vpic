#include <particle_pipelines.h>

void
center_p( particle_t           * ALIGNED p,
          const int                      n,
          const float                    q_m,
          const interpolator_t * ALIGNED f,
          const grid_t         *         g ) {
  center_p_pipeline_args_t args[1];

  if( n<0     ) ERROR(("Bad number of particles"));
  if( f==NULL ) ERROR(("Bad interpolator"));
  if( g==NULL ) ERROR(("Bad grid"));

  /* Have the pipelines do the bulk of particles in quads and
     have the host do the final incomplete quad. */

  args->p   = p;
  args->n   = n;
  args->q_m = q_m;
  args->f   = f;
  args->g   = g;

  dispatch_pipelines( center_p_pipeline_v4, args, 0 );
  center_p_pipeline( args, _n_pipeline, _n_pipeline );
  wait_for_pipelines();
}
