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

# if 1
  float w0, w1, w2, w3, w4, w5, w6, w7, dz;
# else
  using namespace v4;
  v4float q, wl, wh, rl, rh;
# endif

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

#   if 1
    // After detailed experiments and studying of assembly dumps, it was
    // determined that if the platform does not support efficient 4-vector
    // SIMD memory gather/scatter operations, the savings from using
    // "trilinear" are slightly outweighed by the overhead of the
    // gather/scatters.
 
    // Load the particle data

    w0 = p[i].dx;
    w1 = p[i].dy;
    dz = p[i].dz;
    v  = p[i].i;
    w7 = p[i].q*r8V;

    // Compute the trilinear weights
    // Though the PPE should have hardware fma/fmaf support, it was
    // measured to be more efficient _not_ to use it here.  (Maybe the
    // compiler isn't actually generating the assembly for it.

#   define FMA( x,y,z) ((z)+(x)*(y))
#   define FNMS(x,y,z) ((z)-(x)*(y))
    w6=FNMS(w0,w7,w7);                    // q(1-dx)
    w7=FMA( w0,w7,w7);                    // q(1+dx)
    w4=FNMS(w1,w6,w6); w5=FNMS(w1,w7,w7); // q(1-dx)(1-dy), q(1+dx)(1-dy)
    w6=FMA( w1,w6,w6); w7=FMA( w1,w7,w7); // q(1-dx)(1+dy), q(1+dx)(1+dy)
    w0=FNMS(dz,w4,w4); w1=FNMS(dz,w5,w5); w2=FNMS(dz,w6,w6); w3=FNMS(dz,w7,w7);
    w4=FMA( dz,w4,w4); w5=FMA( dz,w5,w5); w6=FMA( dz,w6,w6); w7=FMA( dz,w7,w7);
#   undef FNMS
#   undef FMA

    // Reduce the particle charge to rhof

    f[v      ].rhof += w0; f[v      +1].rhof += w1;
    f[v   +sy].rhof += w2; f[v   +sy+1].rhof += w3;
    f[v+sz   ].rhof += w4; f[v+sz   +1].rhof += w5;
    f[v+sz+sy].rhof += w6; f[v+sz+sy+1].rhof += w7;

#   else

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

#   endif

  }
}

