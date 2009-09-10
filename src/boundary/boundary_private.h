#ifndef _boundary_private_h_
#define _boundary_private_h_

#ifndef IN_boundary
#error "Do not include boundary_private.h; use boundary.h"
#endif

#include "boundary.h"

BEGIN_C_DECLS

void
checkpt_particle_bc_internal( const particle_bc_t * pbc );

particle_bc_t *
restore_particle_bc_internal( void * params );

particle_bc_t *
new_particle_bc_internal( particle_bc_func_t interact,
                          void * params,
                          delete_particle_bc_func_t delete_pbc,
                          checkpt_func_t checkpt,
                          restore_func_t restore,
                          reanimate_func_t reanimate );

void
delete_particle_bc_internal( particle_bc_t * pbc );

END_C_DECLS

#endif /* _boundary_h_ */

