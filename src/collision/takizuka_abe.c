#define IN_collision
#include "takizuka_abe.h"

/* Private interface *********************************************************/

void
apply_takizuka_abe( takizuka_abe_t * cm ){
  if( cm->interval<1 || (cm->spi->g->step % cm->interval) ) return;
  if( cm->spi->last_sorted!=cm->spi->g->step ) sort_p( cm->spi );
  if( cm->spj->last_sorted!=cm->spi->g->step ) sort_p( cm->spj );
  apply_takizuka_abe_pipeline(cm);
}

void
checkpt_takizuka_abe( const collision_op_t * cop ) {
  const takizuka_abe_t * cm =
    (const takizuka_abe_t *)cop->params;
  CHECKPT( cm, 1 );
  CHECKPT_STR( cm->name );
  CHECKPT_PTR( cm->spi );
  CHECKPT_PTR( cm->spj );
  CHECKPT_PTR( cm->rp );
  checkpt_collision_op_internal( cop );
}

collision_op_t *
restore_takizuka_abe( void ) {
  takizuka_abe_t * cm;
  RESTORE( cm );
  RESTORE_STR( cm->name );
  RESTORE_PTR( cm->spi );
  RESTORE_PTR( cm->spj );
  RESTORE_PTR( cm->rp );
  return restore_collision_op_internal( cm );
}

void
delete_takizuka_abe( collision_op_t * cop ) {
  takizuka_abe_t * cm = (takizuka_abe_t *)cop->params;
  // UNREGISTER_OBJECT( cm );
  FREE( cm->name );
  FREE( cm );
  delete_collision_op_internal( cop );
}

/* Public interface **********************************************************/

collision_op_t *
takizuka_abe( const char       * RESTRICT name,
              /**/  species_t  * RESTRICT spi,
              /**/  species_t  * RESTRICT spj,
              /**/  rng_pool_t * RESTRICT rp,
              const double                nu0,
              const int                   interval ) {
  takizuka_abe_t * cm;
  size_t len = name ? strlen(name) : 0;

  if( !spi || !spj || spi->g!=spj->g || !rp || rp->n_rng<N_PIPELINE ) ERROR(( "Bad args" ));
  if( len==0 ) ERROR(( "Cannot specify a nameless collision model" ));

  MALLOC( cm, 1 );
  MALLOC( cm->name, len+1 );
  strcpy( cm->name, name );
  cm->spi      = spi;
  cm->spj      = spj;
  cm->rp       = rp;
  cm->nu0      = nu0;
  cm->interval = interval;

  return new_collision_op_internal( cm,
                                    (collision_op_func_t)apply_takizuka_abe,
                                    delete_takizuka_abe,
                                    (checkpt_func_t)checkpt_takizuka_abe,
                                    (restore_func_t)restore_takizuka_abe,
                                    NULL );
}
