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

#include "../species_advance.h"

// accumulate_rho_p adds the charge density associated with the
// supplied particle array to the rhof of the fields.  Trilinear
// interpolation is used.  rhof is known at the nodes at the same time
// as particle positions.  No effort is made to fix up edges of the
// computational domain; see note in synchronize_rhob about why this
// is done this way.  All particles on the list must be inbounds.

void
accumulate_rho_p( /**/  field_array_t * RESTRICT fa,
                  const species_t     * RESTRICT sp ) {
  if( !fa || !sp || fa->g!=sp->g ) ERROR(( "Bad args" ));

  /**/  field_t    * RESTRICT ALIGNED(128) f = fa->f;
  const particle_t * RESTRICT ALIGNED(128) p = sp->p;

  const float q_8V = sp->q*sp->g->r8V;
  const int np = sp->np;
  const int sy = sp->g->sy;
  const int sz = sp->g->sz;

# if 1
  float w0, w1, w2, w3, w4, w5, w6, w7, dz;
# else
  using namespace v4;
  v4float q, wl, wh, rl, rh;
# endif

  int n, v;

  // Load the grid data

  for( n=0; n<np; n++ ) {

#   if 1
    // After detailed experiments and studying of assembly dumps, it was
    // determined that if the platform does not support efficient 4-vector
    // SIMD memory gather/scatter operations, the savings from using
    // "trilinear" are slightly outweighed by the overhead of the
    // gather/scatters.
 
    // Load the particle data

    w0 = p[n].dx;
    w1 = p[n].dy;
    dz = p[n].dz;
    v  = p[n].i;
    w7 = p[n].w*q_8V;

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

    v = p[n].i;
    rl = v4float( f[v      ].rhof, f[v      +1].rhof,
                  f[v   +sy].rhof, f[v   +sy+1].rhof);
    rh = v4float( f[v+sz   ].rhof, f[v+sz   +1].rhof,
                  f[v+sz+sy].rhof, f[v+sz+sy+1].rhof);

    // Compute the trilinear weights

    load_4x1( &p[n].dx, wl );
    trilinear( wl, wh );
    
    // Reduce the particle charge to rhof and scatter the result

    q = v4float( p[n].w*q_8V );
    store_4x1_tr( fma(q,wl,rl), &f[v      ].rhof, &f[v      +1].rhof,
                                &f[v   +sy].rhof, &f[v   +sy+1].rhof );
    store_4x1_tr( fma(q,wh,rh), &f[v+sz   ].rhof, &f[v+sz   +1].rhof,
                                &f[v+sz+sy].rhof, &f[v+sz+sy+1].rhof );

#   endif

  }
}

#if 0
using namespace v4;
// Note: If part of the body of accumulate_rhob, under the hood
// there is a check for initialization that occurs everytime
// accumulate_rhob is called!
static const v4float ax[4] = { v4float(1,1,1,1), v4float(2,1,2,1),
                               v4float(1,2,1,2), v4float(2,2,2,2) };
static const v4float ay[4] = { v4float(1,1,1,1), v4float(2,2,1,1),
                               v4float(1,1,2,2), v4float(2,2,2,2) };
#endif

void
accumulate_rhob( field_t          * RESTRICT ALIGNED(128) f,
                 const particle_t * RESTRICT ALIGNED(32)  p,
                 const grid_t     * RESTRICT              g,
                 const float                              qsp ) {
# if 1

  // See note in rhof for why this variant is used.
  float w0 = p->dx, w1 = p->dy, w2, w3, w4, w5, w6, w7, dz = p->dz;
  int v = p->i, x, y, z, sy = g->sy, sz = g->sz;
  w7 = (qsp*g->r8V)*p->w;

  // Compute the trilinear weights
  // See note in rhof for why FMA and FNMS are done this way.

# define FMA( x,y,z) ((z)+(x)*(y))
# define FNMS(x,y,z) ((z)-(x)*(y))
  w6=FNMS(w0,w7,w7);                    // q(1-dx)
  w7=FMA( w0,w7,w7);                    // q(1+dx)
  w4=FNMS(w1,w6,w6); w5=FNMS(w1,w7,w7); // q(1-dx)(1-dy), q(1+dx)(1-dy)
  w6=FMA( w1,w6,w6); w7=FMA( w1,w7,w7); // q(1-dx)(1+dy), q(1+dx)(1+dy)
  w0=FNMS(dz,w4,w4); w1=FNMS(dz,w5,w5); w2=FNMS(dz,w6,w6); w3=FNMS(dz,w7,w7);
  w4=FMA( dz,w4,w4); w5=FMA( dz,w5,w5); w6=FMA( dz,w6,w6); w7=FMA( dz,w7,w7);
# undef FNMS
# undef FMA

  // Adjust the weights for a corrected local accumulation of rhob.
  // See note in synchronize_rho why we must do this for rhob and not
  // for rhof.

  x  = v;    z = x/sz;
  if( z==1     ) w0 += w0, w1 += w1, w2 += w2, w3 += w3;
  if( z==g->nz ) w4 += w4, w5 += w5, w6 += w6, w7 += w7;
  x -= sz*z; y = x/sy;
  if( y==1     ) w0 += w0, w1 += w1, w4 += w4, w5 += w5;
  if( y==g->ny ) w2 += w2, w3 += w3, w6 += w6, w7 += w7;
  x -= sy*y;
  if( x==1     ) w0 += w0, w2 += w2, w4 += w4, w6 += w6;
  if( x==g->nx ) w1 += w1, w3 += w3, w5 += w5, w7 += w7;

  // Reduce the particle charge to rhob

  f[v      ].rhob += w0; f[v      +1].rhob += w1;
  f[v   +sy].rhob += w2; f[v   +sy+1].rhob += w3;
  f[v+sz   ].rhob += w4; f[v+sz   +1].rhob += w5;
  f[v+sz+sy].rhob += w6; f[v+sz+sy+1].rhob += w7;

# else

  v4float q, wl, wh, rl, rh;
  int v, sy = g->sy, sz = g->sz;
  int i, j;

  // Gather rhob for this voxel

  v = p->i;
  rl = v4float( f[v      ].rhob, f[v      +1].rhob,
                f[v   +sy].rhob, f[v   +sy+1].rhob);
  rh = v4float( f[v+sz   ].rhob, f[v+sz   +1].rhob,
                f[v+sz+sy].rhob, f[v+sz+sy+1].rhob);

  // Compute the trilinear weights

  load_4x1( &p->dx, wl );
  trilinear( wl, wh );

  // Adjust the weights for a corrected local accumulation of rhob.
  // See note in synchronize_rho why we must do this for rhob and not
  // for rhof.  Why yes, this code snippet is branchless and evil.

  i = v;
  j = i/sz; i -= sz*j;       load_4x1( &ax[(j==1    )?3:0], q ); wl *= q;
  /**/                       load_4x1( &ax[(j==g->nz)?3:0], q ); wh *= q;
  j = i/sy; i -= sy*j;
  j = (j==1) + 2*(j==g->ny); load_4x1( &ay[j], q ); wl *= q; wh *= q;
  i = (i==1) + 2*(i==g->nx); load_4x1( &ax[i], q ); wl *= q; wh *= q;

  // Reduce the particle charge to rhof and scatter the result

  q = v4float( (qsp*g->r8V)*p->w );
  store_4x1_tr( fma(q,wl,rl), &f[v      ].rhob, &f[v      +1].rhob,
                              &f[v   +sy].rhob, &f[v   +sy+1].rhob );
  store_4x1_tr( fma(q,wh,rh), &f[v+sz   ].rhob, &f[v+sz   +1].rhob,
                              &f[v+sz+sy].rhob, &f[v+sz+sy+1].rhob );

# endif
}
