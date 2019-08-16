#define IN_boundary
#include "boundary_private.h"

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
  pbc->params = params;
  RESTORE_SYM( pbc->interact );
  RESTORE_SYM( pbc->delete_pbc );
  RESTORE_PTR( pbc->next );
  return pbc;
}

particle_bc_t *
new_particle_bc_internal( void * params,
                          particle_bc_func_t interact,
                          delete_particle_bc_func_t delete_pbc,
                          checkpt_func_t checkpt,
                          restore_func_t restore,
                          reanimate_func_t reanimate ) {
  particle_bc_t * pbc;
  MALLOC( pbc, 1 );
  CLEAR( pbc, 1 );
  pbc->params     = params;
  pbc->interact   = interact;
  pbc->delete_pbc = delete_pbc;
  /* id, next set by append_particle_bc */
  REGISTER_OBJECT( pbc, checkpt, restore, reanimate );
  return pbc;
}

void
delete_particle_bc_internal( particle_bc_t * pbc ) {
  UNREGISTER_OBJECT( pbc );
  FREE( pbc );
}

/* Public interface **********************************************************/

int
num_particle_bc( const particle_bc_t * RESTRICT pbc_list ) {
  return pbc_list ? (-pbc_list->id-2) : 0;
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

particle_bc_t *
append_particle_bc( particle_bc_t * pbc,
                    particle_bc_t ** pbc_list ) {
  if( !pbc || !pbc_list ) ERROR(( "Bad args" ));
  if( pbc->next ) ERROR(( "Particle boundary condition already in a list" ));
  // Assumes reflective/absorbing are -1, -2
  pbc->id   = -3-num_particle_bc( *pbc_list );
  pbc->next = *pbc_list;
  *pbc_list = pbc;
  return pbc;
}

int64_t
get_particle_bc_id( particle_bc_t * pbc ) {
  if( !pbc ) return 0;
  if( pbc==(particle_bc_t *) absorb_particles ) return  absorb_particles;
  if( pbc==(particle_bc_t *)reflect_particles ) return reflect_particles;
  return pbc->id;
}

