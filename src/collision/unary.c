#define IN_collision

#include "unary.h"

/* Private interface *********************************************************/

//----------------------------------------------------------------------------//
// Include various programming model implementation files. For now, this is
// just the pipeline model. When there are more models, probably want these
// to be conditionally included.
//----------------------------------------------------------------------------//

#include "unary_pipeline.c"

//----------------------------------------------------------------------------//
// Top level function to select and call proper apply_unary_collision_model
// function.
//----------------------------------------------------------------------------//

void
apply_unary_collision_model( unary_collision_model_t * cm )
{
  if ( cm->interval < 1                   ||
       ( cm->sp->g->step % cm->interval ) )
  {
    return;
  }

  // Conditionally execute this when more abstractions are available.
  apply_unary_collision_model_pipeline( cm );
}

void
checkpt_unary_collision_model( const collision_op_t * cop )
{
  const unary_collision_model_t * cm =
    ( const unary_collision_model_t * ) cop->params;

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
restore_unary_collision_model( void )
{
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
delete_unary_collision_model( collision_op_t * cop )
{
  unary_collision_model_t * cm = ( unary_collision_model_t * ) cop->params;

  FREE( cm->name );
  FREE( cm );

  delete_collision_op_internal( cop );
}

/* Public interface **********************************************************/

collision_op_t *
unary_collision_model( const char * RESTRICT name,
                       unary_rate_constant_func_t rate_constant,
                       unary_collision_func_t collision,
                       void * RESTRICT params,
                       species_t * RESTRICT sp,
                       rng_pool_t * RESTRICT rp,
                       int interval )
{
  unary_collision_model_t * cm;

  size_t len = name ? strlen(name) : 0;

  if ( !rate_constant         ||
       !collision             ||
       !sp                    ||
       !rp                    ||
       rp->n_rng < N_PIPELINE )
  {
    ERROR( ( "Bad args" ) );
  }

  if ( !len )
  {
    ERROR( ( "Cannot specify a nameless collision model." ) );
  }

  if ( params               &&
       !object_id( params ) )
  {
    ERROR( ( "collision model parameters must be checkpoint registered." ) );
  }

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
                                    ( collision_op_func_t ) apply_unary_collision_model,
                                    delete_unary_collision_model,
                                    ( checkpt_func_t ) checkpt_unary_collision_model,
                                    ( restore_func_t ) restore_unary_collision_model,
                                    NULL );
}
