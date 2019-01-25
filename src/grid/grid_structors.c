/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "grid.h"
#include "grid_structors.h"

/* Though these functions are not part of grid's public API, they must
   not be declared as static */

void
checkpt_grid( const grid_t * g ) {
  CHECKPT( g, 1 );
  if( g->range    ) CHECKPT_ALIGNED( g->range, world_size+1, 16 );
  if( g->neighbor ) CHECKPT_ALIGNED( g->neighbor, 6*g->nv, 128 );
  CHECKPT_PTR( g->mp );
}

grid_t *
restore_grid( void ) {
  grid_t * g;
  RESTORE( g );
  if( g->range    ) RESTORE_ALIGNED( g->range );
  if( g->neighbor ) RESTORE_ALIGNED( g->neighbor );
  RESTORE_PTR( g->mp );
  return g;
}

grid_t *
new_grid( void ) {
  grid_t *g;
  int i;
  MALLOC( g, 1 );
  CLEAR( g, 1 );
  for( i=0; i<27; i++ ) g->bc[i] = anti_symmetric_fields;
  g->bc[BOUNDARY(0,0,0)] = world_rank;
  g->mp = new_mp( 27 );
  REGISTER_OBJECT( g, checkpt_grid, restore_grid, NULL );
  return g;
}

void
delete_grid( grid_t * g ) {
  if( !g ) return;
  UNREGISTER_OBJECT( g );
  FREE_ALIGNED( g->neighbor );
  FREE_ALIGNED( g->range );
  delete_mp( g->mp );
  FREE( g );
}

