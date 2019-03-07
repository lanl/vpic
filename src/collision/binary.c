#define IN_collision

#include "binary.h"

/* FIXME: ADD SAMPLE TO UNARY */

/* Private interface *********************************************************/

//----------------------------------------------------------------------------//
// Top level function to select and call proper apply_binary_collision_model
// function.
//----------------------------------------------------------------------------//

void
apply_binary_collision_model( binary_collision_model_t * cm )
{
  if ( cm->interval < 1 || ( cm->spi->g->step % cm->interval ) )
  {
    return;
  }

  if ( cm->spi->last_sorted != cm->spi->g->step )
  {
    sort_p( cm->spi );
  }

  if ( cm->spj->last_sorted != cm->spi->g->step )
  {
    sort_p( cm->spj );
  }

  // Conditionally execute this when more abstractions are available.
  apply_binary_collision_model_pipeline( cm );
}

void
checkpt_binary_collision_model( const collision_op_t * cop )
{
  const binary_collision_model_t * cm =
    ( const binary_collision_model_t * ) cop->params;

  CHECKPT( cm, 1 );
  CHECKPT_STR( cm->name );
  CHECKPT_SYM( cm->rate_constant );
  CHECKPT_SYM( cm->collision );
  CHECKPT_PTR( cm->params );
  CHECKPT_PTR( cm->spi );
  CHECKPT_PTR( cm->spj );
  CHECKPT_PTR( cm->rp );

  checkpt_collision_op_internal( cop );
}

collision_op_t *
restore_binary_collision_model( void )
{
  binary_collision_model_t * cm;

  RESTORE( cm );
  RESTORE_STR( cm->name );
  RESTORE_SYM( cm->rate_constant );
  RESTORE_SYM( cm->collision );
  RESTORE_PTR( cm->params );
  RESTORE_PTR( cm->spi );
  RESTORE_PTR( cm->spj );
  RESTORE_PTR( cm->rp );

  return restore_collision_op_internal( cm );
}

void
delete_binary_collision_model( collision_op_t * cop )
{
  binary_collision_model_t * cm
    = (binary_collision_model_t *) cop->params;

  FREE( cm->name );
  FREE( cm );

  delete_collision_op_internal( cop );
}

/* Public interface **********************************************************/

collision_op_t *
binary_collision_model( const char * RESTRICT name,
                        binary_rate_constant_func_t rate_constant,
                        binary_collision_func_t collision,
                        void * RESTRICT params,
                        species_t * RESTRICT spi,
                        species_t * RESTRICT spj,
                        rng_pool_t * RESTRICT rp,
                        double sample,
                        int interval )
{
  binary_collision_model_t * cm;

  size_t len = name ? strlen(name) : 0;

  if ( !rate_constant         ||
       !collision             ||
       !spi                   ||
       !spj                   ||
       spi->g != spj->g       ||
       !rp                    ||
       rp->n_rng < N_PIPELINE )
  {
    ERROR( ( "Bad args" ) );
  }

  if ( len == 0 )
  {
    ERROR( ( "Cannot specify a nameless collision model" ) );
  }

  if ( params               &&
       !object_id( params ) )
  {
    ERROR( ( "collision model parameters must be checkpoint registered" ) );
  }

  MALLOC( cm, 1 );
  MALLOC( cm->name, len+1 );

  strcpy( cm->name, name ); 

  cm->rate_constant = rate_constant;
  cm->collision     = collision;
  cm->params        = params;
  cm->spi           = spi;
  cm->spj           = spj;
  cm->rp            = rp;
  cm->sample        = sample;
  cm->interval      = interval;

  return new_collision_op_internal( cm,
                                    ( collision_op_func_t ) apply_binary_collision_model,
                                    delete_binary_collision_model,
                                    ( checkpt_func_t ) checkpt_binary_collision_model,
                                    ( restore_func_t ) restore_binary_collision_model,
                                    NULL );
}
