/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "sf_interface.h"

hydro_t * ALIGNED(128)
new_hydro( grid_t * g ) {
  hydro_t * ALIGNED(128) h;

  if( g==NULL ) ERROR(("Bad grid."));
  if( g->nx<1 || g->ny<1 || g->nz<1 ) ERROR(("Bad resolution."));
  MALLOC_ALIGNED( h, (g->nx+2)*(g->ny+2)*(g->nz+2), 128 );
  clear_hydro(h,g);
  return h;
}

void
delete_hydro( hydro_t * ALIGNED(128) h ) {
  FREE_ALIGNED(h);
}

void
clear_hydro( hydro_t * ALIGNED(128) h,
             const grid_t * g ) {
  if( h==NULL ) ERROR(("Bad hydro"));
  if( g==NULL ) ERROR(("Bad grid"));
  // FIXME: SPU THIS?
  CLEAR( h, (g->nx+2)*(g->ny+2)*(g->nz+2) );
}

/*****************************************************************************/

interpolator_t * ALIGNED(128)
new_interpolator( grid_t * g ) {
  interpolator_t * ALIGNED(128) fi;

  if( g==NULL ) ERROR(("Invalid grid."));
  if( g->nx<1 || g->ny<1 || g->nz<1 ) ERROR(("Invalid grid resolution.")); 
  MALLOC_ALIGNED( fi, (g->nx+2)*(g->ny+2)*(g->nz+2), 128 );
  CLEAR( fi, (g->nx+2)*(g->ny+2)*(g->nz+2) );
  return fi;
}

void
delete_interpolator( interpolator_t * ALIGNED(128) fi ) {
  FREE_ALIGNED( fi );
}

/*****************************************************************************/

accumulator_t * ALIGNED(128)
new_accumulators( grid_t * g ) {
  accumulator_t * ALIGNED(128) a;
  size_t req;

  if( g==NULL ) ERROR(("Bad grid."));
  if( g->nx<1 || g->ny<1 || g->nz<1 ) ERROR(("Bad resolution."));
  
  /**/                        req = serial.n_pipeline;
  if( req<thread.n_pipeline ) req = thread.n_pipeline;
# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  if( req<spu.n_pipeline    ) req = spu.n_pipeline;
# endif
  req++; /* req = 1 + max( {serial,thread,spu}.n_pipeline ) */
  req *= POW2_CEIL((g->nx+2)*(g->ny+2)*(g->nz+2),2);

  MALLOC_ALIGNED( a, req, 128 );
  CLEAR( a, req );

  return a;
}

void
delete_accumulators( accumulator_t * ALIGNED(128) a ) {
  FREE_ALIGNED( a );
}
