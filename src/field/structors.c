/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <string.h> /* For memset */
#include <field.h>

field_t * ALIGNED new_field( const grid_t * RESTRICT g ) {
  field_t * ALIGNED f;
  int req;

  if( g==NULL    ) { ERROR(("Bad grid."));       return NULL; }
  if( g->nx<1 ||
      g->ny<1 ||
      g->nz<1    ) { ERROR(("Bad resolution.")); return NULL; }

  req = (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(field_t);
  f = (field_t * ALIGNED)malloc_aligned( req, preferred_alignment );
  if( f==NULL ) {
    ERROR(("Failed to allocate field."));
    return NULL;
  }
  clear_field( f, g );

  return f;
}

void delete_field( field_t ** ALIGNED f ) {
  if( f==NULL ) return;
  if( *f!=NULL ) free_aligned(*f);
  *f = NULL;
}

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void clear_field( field_t * RESTRICT ALIGNED f,
	          const grid_t * RESTRICT g ) {
  if( f==NULL ) { ERROR(("Bad field")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));  return; }
  memset( f, 0, (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(field_t) );
}

void clear_jf( field_t * RESTRICT ALIGNED f, const grid_t * RESTRICT g ) {
  int nx, ny, nz, x, y, z;
  field_t *f0;

  if( f==NULL ) { ERROR(("Bad field")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));  return; }

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  for( z=0; z<=nz+1; z++ ) {
    for( y=0; y<=ny+1; y++ ) {
      f0 = &f(0,y,z);
      for( x=0; x<=nx+1; x++ ) {
	f0->jfx = 0;
	f0->jfy = 0;
	f0->jfz = 0;
	f0++;
      }
    }
  }
}

void clear_rhof( field_t * RESTRICT ALIGNED f, const grid_t * RESTRICT g ) {
  int nx, ny, nz, x, y, z;
  field_t *f0;

  if( f==NULL ) { ERROR(("Bad field")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));  return; }

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  for( z=0; z<=nz+1; z++ ) {
    for( y=0; y<=ny+1; y++ ) {
      f0 = &f(0,y,z);
      for( x=0; x<=nx+1; x++ ) {
	f0->rhof = 0;
	f0++;
      }
    }
  }
}

/*****************************************************************************/

hydro_t * ALIGNED new_hydro( const grid_t * RESTRICT g ) {
  hydro_t * ALIGNED h;
  int req;

  if( g==NULL    ) { ERROR(("Bad grid."));       return NULL; }
  if( g->nx<1 ||
      g->ny<1 ||
      g->nz<1    ) { ERROR(("Bad resolution.")); return NULL; }

  req = (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(hydro_t);
  h = (hydro_t * ALIGNED)malloc_aligned( req, preferred_alignment );
  if( h==NULL ) {
    ERROR(("Failed to allocate hydro."));
    return NULL;
  }
  clear_hydro(h,g);

  return h;
}

void delete_hydro( hydro_t ** ALIGNED h ) {
  if( h==NULL ) return;
  if( *h!=NULL ) free_aligned(*h);
  *h = NULL;
}

void clear_hydro( hydro_t * RESTRICT ALIGNED h,
	          const grid_t * RESTRICT g ) {
  if( h==NULL ) { ERROR(("Bad hydro")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));  return; }
  memset( h, 0, (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(hydro_t) );
}

/*****************************************************************************/

interpolator_t * ALIGNED new_interpolator( const grid_t * RESTRICT g ) {
  interpolator_t *fi;
  int req;

  if( g==NULL ) {
    ERROR(("Invalid grid."));
    return NULL;
  }
  if( g->nx<1 || g->ny<1 || g->nz<1 ) {
    ERROR(("Invalid grid resolution."));
    return NULL;
  }
  
  req = (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(interpolator_t);
  fi = (interpolator_t * ALIGNED)malloc_aligned( req, preferred_alignment );
  if( fi==NULL ) {
    ERROR(("Failed to allocate interpolator."));
    return NULL;
  }
  clear_interpolator(fi,g);

  return fi;
}

void delete_interpolator( interpolator_t ** ALIGNED fi ) {
  if( fi==NULL ) return;
  if( *fi!=NULL ) free_aligned( fi );
  *fi = NULL;
}

void clear_interpolator( interpolator_t * RESTRICT ALIGNED fi,
                         const grid_t * RESTRICT g ) {
  if( fi==NULL ) { ERROR(("Bad interpolator")); return; }
  if( g==NULL  ) { ERROR(("Bad grid"));         return; }
  memset( fi, 0, (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(interpolator_t) );
}

/*****************************************************************************/

accumulator_t * ALIGNED new_accumulator( const grid_t * RESTRICT g ) {
  accumulator_t *a;
  int req;

  if( g==NULL    ) { ERROR(("Bad grid."));       return NULL; }
  if( g->nx<1 ||
      g->ny<1 ||
      g->nz<1    ) { ERROR(("Bad resolution.")); return NULL; }
  
  req = (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(accumulator_t);
  a = (accumulator_t * ALIGNED)malloc_aligned( req, preferred_alignment );
  if( a==NULL ) {
    ERROR(("Failed to allocate accumulator."));
    return NULL;
  }
  clear_accumulator( a, g );

  return a;
}

void delete_accumulator( accumulator_t ** ALIGNED a ) {
  if( a==NULL ) return;
  if( *a!=NULL ) free_aligned( a );
  *a = NULL;
}

void clear_accumulator( accumulator_t * RESTRICT ALIGNED a,
			const grid_t * RESTRICT g ) {
  if( a==NULL ) { ERROR(("Invalid accumulator")); return; }
  if( g==NULL ) { ERROR(("Invalid grid"));        return; }
  memset( a, 0, (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(accumulator_t) );
}
