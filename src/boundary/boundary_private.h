#ifndef _boundary_private_h_
#define _boundary_private_h_

#ifndef IN_boundary
#error "Do not include boundary_private.h; use boundary.h"
#endif

#include "boundary.h"

/* The interacting particle is removed from the particle array after
   an interact function is called.  So, if you want to preserve the
   interacting particle after the interaction, you must copy it
   over the to particle injector buffer.   It is the responsibility
   of these handlers update rhob according to net charge added and
   removed from the simulation by these functions. */

typedef int                  /* Number of particles injected */
    ( *particle_bc_func_t )( /* The boundary whose ... */
                             void* RESTRICT
                                 b, /* parameters are b was hit by ...  */
                             species_t* RESTRICT
                                 sp, /* a particle from this species ... */
                             particle_t* RESTRICT
                                 p, /* this particle in fact
                                       (position is hit location, momentum
                                       is at time of the hit) ... */
                             particle_mover_t* RESTRICT
                                 pm, /* who had this much displacement
                                        remaining when it hit */
                             particle_injector_t* RESTRICT
                                 pi,     /* Injectors for particles created by
                                            the interaction */
                             int max_pi, /* Max number injections allowed */
                             int face ); /* CONVENIENCE: Which face of the
                                            the voxel containing the above
                                            particle was hit */

typedef void ( *delete_particle_bc_func_t )( particle_bc_t* RESTRICT pbc );

struct particle_bc
{
    void* params;
    particle_bc_func_t interact;
    delete_particle_bc_func_t delete_pbc;
    int64_t id;
    particle_bc_t* next;
};

BEGIN_C_DECLS

void checkpt_particle_bc_internal( const particle_bc_t* pbc );

particle_bc_t* restore_particle_bc_internal( void* params );

particle_bc_t* new_particle_bc_internal( void* params,
                                         particle_bc_func_t interact,
                                         delete_particle_bc_func_t delete_pbc,
                                         checkpt_func_t checkpt,
                                         restore_func_t restore,
                                         reanimate_func_t reanimate );

void delete_particle_bc_internal( particle_bc_t* pbc );

END_C_DECLS

#endif /* _boundary_h_ */
