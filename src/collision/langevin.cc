#define IN_collision

#include "langevin.h"

/* Private interface *********************************************************/

//----------------------------------------------------------------------------//
// Top level function to select and call the proper apply_langevin function.
//----------------------------------------------------------------------------//

void
apply_langevin( langevin_t * l )
{
  if ( l->interval < 1                  ||
       ( l->sp->g->step % l->interval ) )
  {
    return;
  }

  // Conditionally execute this when more abstractions are available.
  apply_langevin_pipeline( l );
}

void
checkpt_langevin( const collision_op_t * cop )
{
  const langevin_t * l = ( const langevin_t * ) cop->params;

  CHECKPT( l, 1 );
  CHECKPT_PTR( l->sp );
  CHECKPT_PTR( l->rp );

  checkpt_collision_op_internal( cop );
}

collision_op_t *
restore_langevin( void )
{
  langevin_t * l;

  RESTORE( l );
  RESTORE_PTR( l->sp );
  RESTORE_PTR( l->rp );

  return restore_collision_op_internal( l );
}

void
delete_langevin( collision_op_t * cop )
{
  FREE( cop->params );

  delete_collision_op_internal( cop );
}

/* Public interface **********************************************************/

collision_op_t *
langevin( float kT,
          float nu,
          species_t * RESTRICT sp,
          rng_pool_t * RESTRICT rp,
          int interval )
{
  langevin_t * l;

  if ( !sp    ||
       !rp    ||
       kT < 0 ||
       nu < 0 )
  {
    ERROR( ( "Bad args" ) );
  }

  MALLOC( l, 1 );

  l->sp       = sp;
  l->rp       = rp;
  l->kT       = kT;
  l->nu       = nu;
  l->interval = interval;

  return new_collision_op_internal( l,
                                    ( collision_op_func_t ) apply_langevin,
                                    delete_langevin,
                                    ( checkpt_func_t ) checkpt_langevin,
                                    ( restore_func_t ) restore_langevin,
                                    NULL );
}
