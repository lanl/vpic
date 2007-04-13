/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <field.h>
#include <pipelines.h>

field_t * ALIGNED new_field( const grid_t * g ) {
  field_t * ALIGNED f;
  int req;

  if( g==NULL ) ERROR(("Bad grid."));
  if( g->nx<1 || g->ny<1 || g->nz<1 ) ERROR(("Bad resolution."));

  req = (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(field_t);
  f = (field_t * ALIGNED)malloc_aligned( req, preferred_alignment );
  if( f==NULL ) ERROR(("Failed to allocate field."));
  clear_field( f, g );

  return f;
}

void delete_field( field_t ** ALIGNED f ) {
  if( f==NULL ) return;
  if( *f!=NULL ) free_aligned(*f);
  *f = NULL;
}

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void clear_field( field_t * ALIGNED f,
	          const grid_t * g ) {
  if( f==NULL ) ERROR(("Bad field"));
  if( g==NULL ) ERROR(("Bad grid"));
  /* FIXME: SPU THIS */
  memset( f, 0, (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(field_t) );
}

void clear_jf( field_t * ALIGNED f, const grid_t * g ) {
  int nx, ny, nz, x, y, z;
  field_t *f0;

  if( f==NULL ) ERROR(("Bad field")); 
  if( g==NULL ) ERROR(("Bad grid"));

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

void clear_rhof( field_t * ALIGNED f, const grid_t * g ) {
  int nx, ny, nz, x, y, z;
  field_t *f0;

  if( f==NULL ) ERROR(("Bad field"));
  if( g==NULL ) ERROR(("Bad grid"));

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

hydro_t * ALIGNED new_hydro( const grid_t * g ) {
  hydro_t * ALIGNED h;
  int req;

  if( g==NULL ) ERROR(("Bad grid."));
  if( g->nx<1 || g->ny<1 || g->nz<1 ) ERROR(("Bad resolution."));

  req = (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(hydro_t);
  h = (hydro_t * ALIGNED)malloc_aligned( req, preferred_alignment );
  if( h==NULL ) ERROR(("Failed to allocate hydro."));
  clear_hydro(h,g);

  return h;
}

void delete_hydro( hydro_t ** ALIGNED h ) {
  if( h==NULL ) return;
  if( *h!=NULL ) free_aligned(*h);
  *h = NULL;
}

void clear_hydro( hydro_t * ALIGNED h,
	          const grid_t * g ) {
  if( h==NULL ) ERROR(("Bad hydro"));
  if( g==NULL ) ERROR(("Bad grid"));
  memset( h, 0, (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(hydro_t) );
}

/*****************************************************************************/

interpolator_t * ALIGNED new_interpolator( const grid_t * g ) {
  interpolator_t *fi;
  int req;

  if( g==NULL ) ERROR(("Invalid grid."));
  if( g->nx<1 || g->ny<1 || g->nz<1 ) ERROR(("Invalid grid resolution."));
  
  req = (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(interpolator_t);
  fi = (interpolator_t * ALIGNED)malloc_aligned( req, preferred_alignment );
  if( fi==NULL ) ERROR(("Failed to allocate interpolator."));
  clear_interpolator(fi,g);

  return fi;
}

void delete_interpolator( interpolator_t ** ALIGNED fi ) {
  if( fi==NULL ) return;
  if( *fi!=NULL ) free_aligned( fi );
  *fi = NULL;
}

void clear_interpolator( interpolator_t * ALIGNED fi,
                         const grid_t * g ) {
  if( fi==NULL ) ERROR(("Bad interpolator"));
  if( g==NULL  ) ERROR(("Bad grid"));
  memset( fi, 0, (g->nx+2)*(g->ny+2)*(g->nz+2)*sizeof(interpolator_t) );
}

/*****************************************************************************/

#ifndef PMETHOD
#define PMETHOD serial
#endif

accumulator_t * ALIGNED
new_accumulators( const grid_t * g ) {
  accumulator_t *a;
  size_t req;

  if( g==NULL ) ERROR(("Bad grid."));
  if( g->nx<1 || g->ny<1 || g->nz<1 ) ERROR(("Bad resolution."));
  
  req = (size_t)(g->nx+2)*(size_t)(g->ny+2)*(size_t)(g->nz+2)*(size_t)(1+PMETHOD.n_pipeline)*sizeof(accumulator_t);
  a = (accumulator_t * ALIGNED)malloc_aligned( req, preferred_alignment );
  if( a==NULL ) ERROR(("Failed to allocate accumulator."));
  clear_accumulators( a, g );

  return a;
}

void
delete_accumulators( accumulator_t ** ALIGNED a ) {
  if( a==NULL ) return;
  if( *a!=NULL ) free_aligned( a );
  *a = NULL;
}

