#include <grid.h>

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

  if( g==NULL || handler==NULL || sz<0 || sz>MAX_BOUNDARY_DATA_SIZE ||
     ( initial_params==NULL && sz!=0 ) )
    return INVALID_BOUNDARY;

  if( g->boundary==NULL )
    g->boundary = (boundary_t *)malloc( (g->nb+1)*sizeof(boundary_t) );
  else
    g->boundary = (boundary_t *)realloc( g->boundary,
                                         (g->nb+1)*sizeof(boundary_t) );

  if( g->boundary==NULL ) ERROR(( "g->boundary alloc failed" ));

  g->boundary[g->nb].handler = handler;
  memset( g->boundary[g->nb].params, 0, MAX_BOUNDARY_DATA_SIZE );
  memcpy( g->boundary[g->nb].params, initial_params, sz );

  return -(g->nb++) - 3;
}
