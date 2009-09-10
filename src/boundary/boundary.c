#include "boundary.h"

/* Private interface *********************************************************/

void
checkpt_particle_bc_internal( const particle_bc_t * RESTRICT pbc ) {
  CHECKPT( pbc, 1 );
  CHECKPT_SYM( pbc->interact );
  CHECKPT_SYM( pbc->delete_pbc );
  CHECKPT_PTR( pbc->next );
}

particle_bc_t *
restore_particle_bc_internal( void * params ) {
  particle_bc_t * pbc;
  RESTORE( pbc );
  RESTORE_SYM( pbc->interact );
  pbc->params = params;
  RESTORE_ALIGNED( pbc->delete_pbc );
  RESTORE_PTR( pbc->next );
  return pbc;
}

particle_bc_t *
new_particle_bc_internal( particle_bc_func_t interact,
                          void * params,
                          delete_particle_bc_func_t delete_pbc,
                          checkpt_func_t checkpt,
                          restore_func_t restore,
                          reanimate_func_t reanimate ) {
  particle_bc_t * pbc;

  if( !interact || !delete_pbc ) ERROR(( "Bad args" ));

  MALLOC( pbc, 1 );
  CLEAR( pbc, 1 );

  pbc->interact   = interact;
  pbc->params     = params;
  pbc->delete_pbc = delete_pbc;
  pbc->id         = 0;    /* Set by define_particle_bc */
  pbc->next       = NULL; /* Set by define_particle_bc */

  REGISTER_OBJECT( pbc, checkpt, restore, reanimate );
  return pbc;
}

void
delete_particle_bc_internal( particle_bc_t * pbc ) {
  if( !pbc ) return;
  UNREGISTER_OBJECT( pbc );
  FREE( pbc );
}

/* Public interface **********************************************************/

int
num_particle_bc( const particle_bc_t * RESTRICT pbc_list ) {
  return pbc_list ? (-pbc_list->id-3) : 0;
}

void
delete_particle_bc_list( particle_bc_t * pbc_list ) {
  particle_bc_t * pbc;
  while( pbc_list ) {
    pbc = pbc_list;
    pbc_list = pbc_list->next;
    pbc->delete_pbc( pbc );
  }
}

