#define IN_collision

/* #define HAS_V4_PIPELINE */

#include "unary.h"
#include "collision_private.h"

/* Private interface *********************************************************/

void
unary_pipeline_scalar( unary_collision_model_t * RESTRICT cm,
                       int pipeline_rank,
                       int n_pipeline )
{
  if ( pipeline_rank == n_pipeline )
  {
    return; /* No host straggler cleanup */
  }

  unary_rate_constant_func_t rate_constant = cm->rate_constant;
  unary_collision_func_t     collision     = cm->collision;

  /**/  void       * RESTRICT params = cm->params;
  const species_t  * RESTRICT sp     = cm->sp;
  /**/  particle_t * RESTRICT p      = cm->sp->p;
  /**/  rng_t      * RESTRICT rng    = cm->rp->rng[ pipeline_rank ];

  const float dt = sp->g->dt * (float) cm->interval;

  double n_target = (double) sp->np / (double) n_pipeline;

  /**/  int i  = (int) ( 0.5 + n_target * (double)  pipeline_rank    );
  const int i1 = (int) ( 0.5 + n_target * (double) (pipeline_rank+1) );

  float pr_coll;
  int n_large_pr = 0;

  /* For each computational particle assigned to this pipeline, compute
     the probability a comoving physical particle had collision with
     the background.  If this "probability" is greater than one, make
     a note for diagnostic purposes.  Then flip a bias coin of that
     probability to decide if this particle should undergo a collision. */

  for( ; i < i1; i++ )
  {
    pr_coll = dt * rate_constant( params, sp, &p[i] );

    if ( pr_coll > 1 )
    {
      n_large_pr++;
    }

    /* Yes, strictly < (so that 0 rate constants guarantee no collision,
       and, yes, _c0, so that 1 probabilities guarantee a collision  */
    if ( frand_c0( rng ) < pr_coll )
    {
      collision( params, sp, &p[i], rng );
    }
  }

  cm->n_large_pr[ pipeline_rank ] = n_large_pr;
}

void
apply_unary_collision_model_pipeline( unary_collision_model_t * cm )
{
  int p, n_large_pr = 0;

  if ( cm->interval < 1                   ||
       ( cm->sp->g->step % cm->interval ) )
  {
    return;
  }

  EXEC_PIPELINES( unary, cm, 0 );

  WAIT_PIPELINES();

  for( p = 0; p < N_PIPELINE; p++ )
  {
    n_large_pr += cm->n_large_pr[p];
  }

  if ( n_large_pr )
  {
    WARNING( ( "%i particles in species \"%s\" encountered a large collision "
               "probability in collision model \"%s\".  The collision rate for "
               "such particles will be lower than it should be physically.  "
               "Consider lowering the collision operator interval or reducing "
               "the timestep.", n_large_pr, cm->sp->name, cm->name ) );
  }
}
