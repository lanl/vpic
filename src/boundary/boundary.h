#ifndef _boundary_h_
#define _boundary_h_

#include "../species_advance/species_advance.h"

struct particle_bc;

/* The interacting particle is removed from the particle array after
   an interact function is called.  So, if you want to preserve the
   interacting particle after the interaction, you must copy it
   over the to particle injector buffer.   It is the responsibility
   of these handlers update rhob according to net charge added and
   removed from the simulation by these functions. */

typedef int /* Number of particles injected */
(*particle_bc_func_t)(                   /* The boundary whose ... */
  void                * RESTRICT b,      /* parameters are b was hit by ...  */
  species_t           * RESTRICT sp,     /* a particle from this species ... */
  particle_t          * RESTRICT p,      /* this particle in fact
                                            (position is hit location, momentum
                                            is at time of the hit) ... */
  particle_mover_t    * RESTRICT pm,     /* who had this much displacement
                                            remaining when it hit */
  particle_injector_t * RESTRICT pi,     /* Injectors for particles created by
                                            the interaction */
  int                            max_pi, /* Max number injections allowed */
  int                            face ); /* CONVENIENCE: Which face of the
                                            the voxel containing the above
                                            particle was hit */

typedef void
(*delete_particle_bc_func_t)( struct particle_bc * RESTRICT pbc );

typedef struct particle_bc {
  particle_bc_func_t interact;
  void * params;
  delete_particle_bc_func_t delete_pbc;
  int64_t id;
  struct particle_bc * next;
} particle_bc_t;

BEGIN_C_DECLS

/* In boundary.c */

int
num_particle_bc( const particle_bc_t * RESTRICT pbc_list );

void
delete_particle_bc_list( particle_bc_t * RESTRICT pbc_list );

/* In boundary_p.cxx */

void
boundary_p( particle_bc_t       * RESTRICT pbc_list,
            species_t           * RESTRICT sp_list,
            field_array_t       * RESTRICT fa,
            accumulator_array_t * RESTRICT aa );

/* In maxwellian_reflux.c */

/* DO NOT CALL DIRECTLY ... CALL THROUGH define_particle_bc */
particle_bc_t *
maxwellian_reflux( species_t * RESTRICT sp_list,
                   mt_rng_t  * RESTRICT rng );

void
set_reflux_temp( /**/  particle_bc_t * RESTRICT mr,
                 const species_t     * RESTRICT sp,
                 float ut_para,
                 float ut_perp );

/* In absorb_tally.c */

/* DO NOT CALL DIRECTLY ... CALL THROUGH define_particle_bc */
particle_bc_t *
absorb_tally( /**/  species_t      * RESTRICT sp_list,
              const field_array_t  * RESTRICT fa );

int *
get_absorb_tally( particle_bc_t * pbc );

END_C_DECLS

#endif /* _boundary_h_ */

