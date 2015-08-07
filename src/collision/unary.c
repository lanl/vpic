#define IN_collision
/*#define HAS_V4_PIPELINE*/
#include "collision_private.h"

/* Private interface *********************************************************/

typedef struct unary_collision_model {
  char * name;
  unary_rate_constant_func_t rate_constant;
  unary_collision_func_t collision;
  void * params;
  species_t * sp;
  rng_pool_t * rp;
  int interval;
  int n_large_pr[ MAX_PIPELINE ];
} unary_collision_model_t;

void
unary_pipeline( unary_collision_model_t * RESTRICT cm,
                 int pipeline_rank,
                 int n_pipeline ) {
  if( pipeline_rank==n_pipeline ) return; /* No host straggler cleanup */

  unary_rate_constant_func_t rate_constant = cm->rate_constant;
  unary_collision_func_t     collision     = cm->collision;

  /**/  void       * RESTRICT params = cm->params;
  const species_t  * RESTRICT sp     = cm->sp;
  /**/  particle_t * RESTRICT p      = cm->sp->p;
  /**/  rng_t      * RESTRICT rng    = cm->rp->rng[ pipeline_rank ];
  const float dt = sp->g->dt * (float)cm->interval;

  double n_target = (double)sp->np / (double)n_pipeline;
  /**/  int i  = (int)( 0.5 + n_target*(double) pipeline_rank    );
  const int i1 = (int)( 0.5 + n_target*(double)(pipeline_rank+1) );

  float pr_coll;
  int n_large_pr = 0;

  /* For each computational particle assigned to this pipeline, compute
     the probability a comoving physical particle had collision with
     the background.  If this "probability" is greater than one, make
     a note for diagnostic purposes.  Then flip a bias coin of that
     probability to decide if this particle should undergo a collision. */

  for( ; i<i1; i++ ) {
    pr_coll = dt * rate_constant( params, sp, &p[i] );
    if( pr_coll>1 ) n_large_pr++;

    /* Yes, strictly < (so that 0 rate constants guarantee no collision,
       and, yes, _c0, so that 1 probabilities guarantee a collision  */
    if( frand_c0(rng) < pr_coll ) collision( params, sp, &p[i], rng );
  }

  cm->n_large_pr[pipeline_rank] = n_large_pr;
}

void
apply_unary_collision_model( unary_collision_model_t * cm ) {
  int p, n_large_pr = 0;
  if( cm->interval<1 || (cm->sp->g->step % cm->interval) ) return;
  EXEC_PIPELINES( unary, cm, 0 );
  WAIT_PIPELINES();
  for( p=0; p<N_PIPELINE; p++ ) n_large_pr += cm->n_large_pr[p];
  if( n_large_pr )
    WARNING(( "%i particles in species \"%s\" encountered a large collision "
              "probability in collision model \"%s\".  The collision rate for "
              "such particles will be lower than it should be physically.  "
              "Consider lowering the collision operator interval or reducing "
              "the timestep.", n_large_pr, cm->sp->name, cm->name ));
}

void
checkpt_unary_collision_model( const collision_op_t * cop ) {
  const unary_collision_model_t * cm =
    (const unary_collision_model_t *)cop->params;
  CHECKPT( cm, 1 );
  CHECKPT_STR( cm->name );
  CHECKPT_SYM( cm->rate_constant );
  CHECKPT_SYM( cm->collision );
  CHECKPT_PTR( cm->params );
  CHECKPT_PTR( cm->sp );
  CHECKPT_PTR( cm->rp );
  checkpt_collision_op_internal( cop );
}

collision_op_t *
restore_unary_collision_model( void ) {
  unary_collision_model_t * cm;
  RESTORE( cm );
  RESTORE_STR( cm->name );
  RESTORE_SYM( cm->rate_constant );
  RESTORE_SYM( cm->collision );
  RESTORE_PTR( cm->params );
  RESTORE_PTR( cm->sp );
  RESTORE_PTR( cm->rp );
  return restore_collision_op_internal( cm );
}

void
delete_unary_collision_model( collision_op_t * cop ) {
  unary_collision_model_t * cm = (unary_collision_model_t *)cop->params;
  FREE( cm->name );
  FREE( cm );
  delete_collision_op_internal( cop );
}

/* Public interface **********************************************************/

collision_op_t *
unary_collision_model( const char       * RESTRICT name,
                       unary_rate_constant_func_t  rate_constant,
                       unary_collision_func_t      collision,
                       /**/  void       * RESTRICT params,
                       /**/  species_t  * RESTRICT sp,
                       /**/  rng_pool_t * RESTRICT rp,
                       int                         interval ) {
  unary_collision_model_t * cm;
  size_t len = name ? strlen(name) : 0;

  if( !rate_constant || !collision || !sp || !rp || rp->n_rng<N_PIPELINE )
    ERROR(( "Bad args" ));
  if( !len ) ERROR(( "Cannot specify a nameless collision model" ));
  if( params && !object_id( params ) )
    ERROR(( "collision model parameters must be checkpoint registered" ));

  MALLOC( cm, 1 );
  MALLOC( cm->name, len+1 );
  strcpy( cm->name, name ); 
  cm->rate_constant = rate_constant;
  cm->collision     = collision;
  cm->params        = params;
  cm->sp            = sp;
  cm->rp            = rp;
  cm->interval      = interval;
  return new_collision_op_internal( cm,
                                    (collision_op_func_t)apply_unary_collision_model,
                                    delete_unary_collision_model,
                                    (checkpt_func_t)checkpt_unary_collision_model,
                                    (restore_func_t)restore_unary_collision_model,
                                    NULL );
}

