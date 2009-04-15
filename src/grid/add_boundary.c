#include "grid.h"

// Define custom boundary handler for grid g.  Return integer tag in
// range [-3...-inf] associated with boundary handler which is used
// with set_pbc in the same context as absorb_particles,
// reflect_particles, etc.  to associate boundary handler with
// boundary.  Space for handlers is automatically allocated.

int
IUO_add_boundary( grid_t *g,
                  boundary_handler_t handler,
                  const void * initial_params,
                  int sz ) {
  boundary_t * gb;

  if( g==NULL || handler==NULL || sz<0 || sz>MAX_BOUNDARY_DATA_SIZE ||
     ( initial_params==NULL && sz!=0 ) ) {
	MESSAGE(("Add boundary encountered invalid boundary!!!"));
    return INVALID_BOUNDARY;
  } // if

  MALLOC( gb, g->nb+1 );
  COPY( gb, g->boundary, g->nb );
  FREE( g->boundary );
  g->boundary = gb;

  g->boundary[g->nb].handler = handler;
  CLEAR( (char *)g->boundary[g->nb].params, MAX_BOUNDARY_DATA_SIZE );
  COPY(  (char *)g->boundary[g->nb].params, initial_params, sz );

  return -(g->nb++) - 3;
}
