// FIXME: THREAD THIS! HYDRO MEM SEMANTICS WILL NEED UPDATING.
// FIXME: V4 ACCELERATE THIS.  COULD BE BASED OFF ENERGY_P.

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

// accumulate_hydro_p adds the hydrodynamic fields associated with the
// supplied particle_list to the hydro array.  Trilinear interpolation
// is used.  hydro is known at the nodes at the same time as particle
// positions. No effort is made to fix up edges of the computational
// domain.  All particles on the list must be inbounds.  Note, the
// hydro jx,jy,jz are for diagnostic purposes only; they are not
// accumulated with a charge conserving algorithm.

void
accumulate_hydro_p( hydro_array_t              * RESTRICT ha,
                    const species_t            * RESTRICT sp,
                    const interpolator_array_t * RESTRICT ia ) {
  /**/  hydro_t        * RESTRICT ALIGNED(128) h;
  const particle_t     * RESTRICT ALIGNED(128) p;
  const interpolator_t * RESTRICT ALIGNED(128) f;
  float c, qsp, mspc, qdt_2mc, qdt_4mc2, r8V;
  int np, stride_10, stride_21, stride_43;

  float dx, dy, dz, ux, uy, uz, w, vx, vy, vz, ke_mc;
  float w0, w1, w2, w3, w4, w5, w6, w7, t;
  int i, n;

  if( !ha || !sp || !ia || ha->g!=sp->g || ha->g!=ia->g )
    ERROR(( "Bad args" ));

  h = ha->h;
  p = sp->p;
  f = ia->i;

  c        = sp->g->cvac;
  qsp      = sp->q;
  mspc     = sp->m*c;
  qdt_2mc  = (qsp*sp->g->dt)/(2*mspc);
  qdt_4mc2 = qdt_2mc / (2*c);
  r8V      = sp->g->r8V;

  np        = sp->np;
  stride_10 = VOXEL(1,0,0, sp->g->nx,sp->g->ny,sp->g->nz) -
              VOXEL(0,0,0, sp->g->nx,sp->g->ny,sp->g->nz);
  stride_21 = VOXEL(0,1,0, sp->g->nx,sp->g->ny,sp->g->nz) -
              VOXEL(1,0,0, sp->g->nx,sp->g->ny,sp->g->nz);
  stride_43 = VOXEL(0,0,1, sp->g->nx,sp->g->ny,sp->g->nz) -
              VOXEL(1,1,0, sp->g->nx,sp->g->ny,sp->g->nz);

  for( n=0; n<np; n++ ) {

    // Load the particle
    dx = p[n].dx;
    dy = p[n].dy;
    dz = p[n].dz;
    i  = p[n].i;
    ux = p[n].ux;
    uy = p[n].uy;
    uz = p[n].uz;
    w  = p[n].w;
    
    // Half advance E
    ux += qdt_2mc*((f[i].ex+dy*f[i].dexdy) + dz*(f[i].dexdz+dy*f[i].d2exdydz));
    uy += qdt_2mc*((f[i].ey+dz*f[i].deydz) + dx*(f[i].deydx+dz*f[i].d2eydzdx));
    uz += qdt_2mc*((f[i].ez+dx*f[i].dezdx) + dy*(f[i].dezdy+dx*f[i].d2ezdxdy));

    // Boris rotation - Interpolate B field
    w5 = f[i].cbx + dx*f[i].dcbxdx;
    w6 = f[i].cby + dy*f[i].dcbydy;
    w7 = f[i].cbz + dz*f[i].dcbzdz;

    // Boris rotation - curl scalars (0.5 in v0 for half rotate) and
    // kinetic energy computation. Note: gamma-1 = |u|^2 / (gamma+1)
    // is the numerically accurate way to compute gamma-1
    ke_mc = ux*ux + uy*uy + uz*uz; // ke_mc = |u|^2 (invariant)
    vz = sqrt(1+ke_mc);            // vz = gamma    (invariant)
    ke_mc *= c/(vz+1);             // ke_mc = c|u|^2/(gamma+1) = c*(gamma-1)

#ifdef ENABLE_NON_RELATIVSTIC
    vz = 1.0f;
#else
    vz = c/vz;                     // vz = c/gamma
#endif

    w0 = qdt_4mc2*vz;
    w1 = w5*w5 + w6*w6 + w7*w7;    // |cB|^2
    w2 = w0*w0*w1;
    w3 = w0*(1+(1./3.)*w2*(1+0.4*w2));
    w4 = w3/(1 + w1*w3*w3); w4 += w4;

    // Boris rotation - uprime
    w0 = ux + w3*( uy*w7 - uz*w6 );
    w1 = uy + w3*( uz*w5 - ux*w7 );
    w2 = uz + w3*( ux*w6 - uy*w5 );

    // Boris rotation - u
    ux += w4*( w1*w7 - w2*w6 );
    uy += w4*( w2*w5 - w0*w7 );
    uz += w4*( w0*w6 - w1*w5 );

    // Compute physical velocities
    vx  = ux*vz;
    vy  = uy*vz;
    vz *= uz;

    // Compute the trilinear coefficients
    w0  = r8V*w;    // w0 = (1/8)(w/V)
    dx *= w0;       // dx = (1/8)(w/V) x
    w1  = w0+dx;    // w1 = (1/8)(w/V) + (1/8)(w/V)x = (1/8)(w/V)(1+x)
    w0 -= dx;       // w0 = (1/8)(w/V) - (1/8)(w/V)x = (1/8)(w/V)(1-x)
    w3  = 1+dy;     // w3 = 1+y
    w2  = w0*w3;    // w2 = (1/8)(w/V)(1-x)(1+y)
    w3 *= w1;       // w3 = (1/8)(w/V)(1+x)(1+y)
    dy  = 1-dy;     // dy = 1-y
    w0 *= dy;       // w0 = (1/8)(w/V)(1-x)(1-y)
    w1 *= dy;       // w1 = (1/8)(w/V)(1+x)(1-y)
    w7  = 1+dz;     // w7 = 1+z
    w4  = w0*w7;    // w4 = (1/8)(w/V)(1-x)(1-y)(1+z) = (w/V) trilin_0 *Done
    w5  = w1*w7;    // w5 = (1/8)(w/V)(1+x)(1-y)(1+z) = (w/V) trilin_1 *Done
    w6  = w2*w7;    // w6 = (1/8)(w/V)(1-x)(1+y)(1+z) = (w/V) trilin_2 *Done
    w7 *= w3;       // w7 = (1/8)(w/V)(1+x)(1+y)(1+z) = (w/V) trilin_3 *Done
    dz  = 1-dz;     // dz = 1-z
    w0 *= dz;       // w0 = (1/8)(w/V)(1-x)(1-y)(1-z) = (w/V) trilin_4 *Done
    w1 *= dz;       // w1 = (1/8)(w/V)(1+x)(1-y)(1-z) = (w/V) trilin_5 *Done
    w2 *= dz;       // w2 = (1/8)(w/V)(1-x)(1+y)(1-z) = (w/V) trilin_6 *Done
    w3 *= dz;       // w3 = (1/8)(w/V)(1+x)(1+y)(1-z) = (w/V) trilin_7 *Done

    // Accumulate the hydro fields
#   define ACCUM_HYDRO( wn)                             \
    t  = qsp*wn;        /* t  = (qsp w/V) trilin_n */   \
    h[i].jx  += t*vx;                                   \
    h[i].jy  += t*vy;                                   \
    h[i].jz  += t*vz;                                   \
    h[i].rho += t;                                      \
    t  = mspc*wn;       /* t = (msp c w/V) trilin_n */  \
    dx = t*ux;          /* dx = (px w/V) trilin_n */    \
    dy = t*uy;                                          \
    dz = t*uz;                                          \
    h[i].px  += dx;                                     \
    h[i].py  += dy;                                     \
    h[i].pz  += dz;                                     \
    h[i].ke  += t*ke_mc;                                \
    h[i].txx += dx*vx;                                  \
    h[i].tyy += dy*vy;                                  \
    h[i].tzz += dz*vz;                                  \
    h[i].tyz += dy*vz;                                  \
    h[i].tzx += dz*vx;                                  \
    h[i].txy += dx*vy

    /**/            ACCUM_HYDRO(w0); // Cell i,j,k
    i += stride_10; ACCUM_HYDRO(w1); // Cell i+1,j,k
    i += stride_21; ACCUM_HYDRO(w2); // Cell i,j+1,k
    i += stride_10; ACCUM_HYDRO(w3); // Cell i+1,j+1,k
    i += stride_43; ACCUM_HYDRO(w4); // Cell i,j,k+1
    i += stride_10; ACCUM_HYDRO(w5); // Cell i+1,j,k+1
    i += stride_21; ACCUM_HYDRO(w6); // Cell i,j+1,k+1
    i += stride_10; ACCUM_HYDRO(w7); // Cell i+1,j+1,k+1

#   undef ACCUM_HYDRO
  }
}
