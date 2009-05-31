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

grid_t *
new_grid( void ) {
  grid_t *g;
  int i;
  MALLOC( g, 1 );
  CLEAR( g, 1 );
  g->mp = new_mp();
  for( i=0; i<27; i++ ) g->bc[i] = anti_symmetric_fields;
  g->bc[BOUNDARY(0,0,0)] = mp_rank(g->mp);
  return g;
}

void
delete_grid( grid_t * g ) {
  if( g==NULL ) return;
  delete_mp( &(g->mp) );
  FREE_ALIGNED( g->range );
  FREE_ALIGNED( g->neighbor );
  FREE( g->boundary );
  FREE( g );
}

