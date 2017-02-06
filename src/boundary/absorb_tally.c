#define IN_boundary
#include "boundary_private.h"

// Special absorbing boundary condition that tallies the number of
// particles absorbed of each species.
//
// Written by:  Brian J. Albright, X-1, LANL   July, 2005
// Major revamp KJB: Sep 2009

/* Private interface ********************************************************/

typedef struct absorb_tally {
  /**/  species_t     * sp_list;
  const field_array_t * fa;
  /**/  int           * tally;
} absorb_tally_t;

int
interact_absorb_tally( absorb_tally_t      * RESTRICT at,
                       species_t           * RESTRICT sp,
                       particle_t          * RESTRICT p,
                       particle_mover_t    * RESTRICT pm,
                       particle_injector_t * RESTRICT pi,
                       int                            max_pi,
                       int                            face ) {
  at->tally[ sp->id ]++;
  accumulate_rhob( at->fa->f, p, at->fa->g, sp->q );
  return 0;
}

void
checkpt_absorb_tally( const particle_bc_t * RESTRICT pbc ) {
  const absorb_tally_t * RESTRICT at = (const absorb_tally_t *)pbc->params;
  CHECKPT( at, 1 );
  CHECKPT_PTR( at->sp_list );
  CHECKPT_PTR( at->fa );
  CHECKPT( at->tally, num_species( at->sp_list ) );
  checkpt_particle_bc_internal( pbc );
}

particle_bc_t *
restore_absorb_tally( void ) {
  absorb_tally_t * at;
  RESTORE( at );
  RESTORE_PTR( at->sp_list );
  RESTORE_PTR( at->fa );
  RESTORE( at->tally );
  return restore_particle_bc_internal( at );
}

void
delete_absorb_tally( particle_bc_t * RESTRICT pbc ) {
  absorb_tally_t * at = (absorb_tally_t *)pbc->params;
  FREE( at->tally );
  FREE( at );
  delete_particle_bc_internal( pbc );
}

/* Publc interface **********************************************************/

particle_bc_t *
absorb_tally( /**/  species_t      * RESTRICT sp_list,
              const field_array_t  * RESTRICT fa ) {
  if( !sp_list || !fa ) ERROR(( "Bad args" ));
  absorb_tally_t * at;
  MALLOC( at, 1 );
  at->sp_list = sp_list;
  at->fa      = fa;
  MALLOC( at->tally, num_species( sp_list ) );
  CLEAR( at->tally, num_species( sp_list ) );
  return new_particle_bc_internal( at,
                                   (particle_bc_func_t)interact_absorb_tally,
                                   delete_absorb_tally,
                                   (checkpt_func_t)checkpt_absorb_tally,
                                   (restore_func_t)restore_absorb_tally,
                                   NULL );
}

int *
get_absorb_tally( particle_bc_t * pbc ) {
  if( !pbc ) ERROR(( "Bad args" ));
  return ((absorb_tally_t *)pbc->params)->tally;
}

