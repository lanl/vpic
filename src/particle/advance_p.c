/* FIXME: PARTICLE MOVERS NEED TO BE OVERALLOCATED IN STRUCTORS TO
   ACCOUNT FOR SPLITTING THE MOVER ARRAY BETWEEN HOST AND PIPELINES */

#include <particle_pipelines.h>

int
advance_p( particle_t           * ALIGNED p,
           const int                      n,
           const float                    q_m,
           particle_mover_t     * ALIGNED pm,
           int                            nm,       
           accumulator_t        * ALIGNED a,
           const interpolator_t * ALIGNED f,
           const grid_t         *         g ) {
  advance_p_pipeline_args_t args[1];
  int rank;

  if( p==NULL  ) ERROR(("Bad particle array"));
  if( n<0      ) ERROR(("Bad number of particles"));
  if( pm==NULL ) ERROR(("Bad particle mover"));
  if( nm<0     ) ERROR(("Bad number of movers"));
  if( a==NULL  ) ERROR(("Bad accumulator"));
  if( f==NULL  ) ERROR(("Bad interpolator"));
  if( g==NULL  ) ERROR(("Bad grid"));

  args->p   = p;
  args->n   = n;
  args->q_m = q_m;
  args->pm  = pm;
  args->nm  = nm;
  args->a   = a;
  args->f   = f;
  args->g   = g;

  dispatch_pipelines( advance_p_pipeline_v4, args, 0 );

  /* Have the host processor do the incomplete quad if necessary.
     Note: This is overlapped with the pipelined processing.  As such,
     it uses an entire accumulator.  Reserving an entirely accumulator
     for the host processor to handle at most 3 particles is wasteful
     of memory.  It is anticipated that it may be useful at some point
     in the future have pipelines accumulating currents while the host
     processor is doing other more substantive work (e.g. accumulating
     currents from particles received from neighboring nodes).
     However, it is worth reconsidering this at some point in the
     future. */

  advance_p_pipeline( args, _n_pipeline, _n_pipeline );

  wait_for_pipelines();

  /* FIXME: HIDEOUS HACK UNTIL BETTER PARTICLE MOVER SEMANTICS
     INSTALLED FOR DEALING WITH PIPELINES.  COMPACT THE PARTICLE
     MOVERS TO ELIMINATE HOLES IN THE ALLOCATION. */

  nm = 0;
  for( rank=0; rank<=_n_pipeline; rank++ ) {
    if( pm+nm!=args->seg[rank].pm )
      memmove( pm+nm, args->seg[rank].pm, args->seg[rank].nm*sizeof(*pm) );
    nm += args->seg[rank].nm;
  }

  return nm;
}
