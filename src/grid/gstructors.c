/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <grid.h>

grid_t *new_grid(void) {
  grid_t *g;
  int i;

  g = (grid_t *)malloc(sizeof(grid_t));
  if( g==NULL ) {
    ERROR(("Could not allocate grid"));
    return NULL;
  }

  g->mp = new_mp();
  if( g->mp==NULL ) {
    ERROR(("Could not create message passer"));
    free(g);
    return NULL;
  }
  g->dt = 0;
  g->cvac = 0;
  g->eps0 = 0;
  g->damp = 0;
  g->x0 = 0;    
  g->y0 = 0;   
  g->z0 = 0; 
  g->dx = 0;
  g->dy = 0;
  g->dz = 0;
  g->nx = 0;
  g->ny = 0;
  g->nz = 0;
  for( i=0; i<27; i++ ) g->bc[i] = anti_symmetric_fields;
  g->bc[BOUNDARY(0,0,0)] = mp_rank(g->mp);
  g->range = NULL;
  g->neighbor = NULL;
  g->rangel = 0;
  g->rangeh = 0;
  g->nb = 0;
  g->boundary = NULL; 
  return g;
}

void delete_grid( grid_t **g ) {
  if( g==NULL ) return;
  delete_mp(&((*g)->mp));
  if( (*g)->range!=NULL ) free_aligned( (*g)->range );
  if( (*g)->neighbor!=NULL ) free_aligned( (*g)->neighbor );
  if ( (*g)->boundary ) free( (*g)->boundary );
  (*g)->range    = NULL;
  (*g)->neighbor = NULL;
  (*g)->boundary = NULL;
  free(*g);
  (*g) = NULL;
}

