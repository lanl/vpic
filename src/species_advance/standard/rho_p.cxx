/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised and extended from earlier V4PIC versions
 *
 */

#define IN_spa
#include "spa_private.h"

using namespace v4;

// accumulate_rho_p adds the charge density associated with the
// supplied particle array to the rhof of the fields.  Trilinear
// interpolation is used.  rhof is known at the nodes at the same time
// as particle positions.  No effort is made to fix up edges of the
// computational domain; see note in synchronize_rhob about why this
// is done this way.  All particles on the list must be inbounds.

void
accumulate_rho_p( field_t          * RESTRICT ALIGNED(16)  f,
                  const particle_t * RESTRICT ALIGNED(128) p,
                  int                                      n,
                  const grid_t     * RESTRICT              g ) {
  v4float q, wl, wh, rl, rh;
  float r8V;
  int i, v, sy, sz;

  if( f==NULL ) ERROR(("Bad field"));
  if( p==NULL ) ERROR(("Bad particle array"));
  if( n<0     ) ERROR(("Bad number of particles"));
  if( g==NULL ) ERROR(("Bad grid"));

  // Load the grid data
  r8V = g->r8V;
  sy  = g->sy;
  sz  = g->sz;

  for( i=0; i<n; i++ ) {
    // Gather rhof for this voxel
    v = p[i].i;
    rl = v4float( f[v      ].rhof, f[v      +1].rhof,
                  f[v   +sy].rhof, f[v   +sy+1].rhof);
    rh = v4float( f[v+sz   ].rhof, f[v+sz   +1].rhof,
                  f[v+sz+sy].rhof, f[v+sz+sy+1].rhof);

    // Compute the trilinear weights
    load_4x1( &p[i].dx, wl );
    trilinear( wl, wh );
    
    // Reduce the particle charge to rhof and scatter the result
    q = v4float( r8V*p[i].q );
    store_4x1_tr( fma(q,wl,rl), &f[v      ].rhof, &f[v      +1].rhof,
                                &f[v   +sy].rhof, &f[v   +sy+1].rhof );
    store_4x1_tr( fma(q,wh,rh), &f[v+sz   ].rhof, &f[v+sz   +1].rhof,
                                &f[v+sz+sy].rhof, &f[v+sz+sy+1].rhof );
  }
}

