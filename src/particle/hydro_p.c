/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised and extended from earlier V4PIC versions
 *
 */

#include <species.h>

/* accumulate_hydro_p adds the hydrodynamic fields associated with the
   supplied particle_list to the hydro array. Trilinear interpolation is used.
   hydro is known at the nodes at the same time as particle positions. No
   effort is made to fix up edges of the computational domain. All particles
   on the list must be inbounds. Note, the hydro jx,jy,jz are for diagnostic
   purposes only; they are not accumulated with a charge conserving
   algorithm. */

#define h(x,y,z) h0[INDEX_FORTRAN_3(x,y,z,0,g->nx+1,0,g->ny+1,0,g->nz+1)]

void accumulate_hydro_p( hydro_t * ALIGNED h0,
                         const particle_t * ALIGNED p,
                         int n,
                         float q_m,
                         const interpolator_t * ALIGNED f0,
                         const grid_t * g ) {
  float dx, dy, dz; int ii;
  float ux, uy, uz, q;
  float vx, vy, vz, ke_mc;
  float w0, w1, w2, w3, w4, w5, w6, w7;
  float qdt_2mc, qdt_4mc2, c, r8V, mc_q;
  const interpolator_t * ALIGNED f;
  hydro_t * ALIGNED h;
  int stride_10, stride_21, stride_43;

  if( h0==NULL ) ERROR(("Bad hydro"));
  if( p==NULL  ) ERROR(("Bad particle array"));
  if( n<0      ) ERROR(("Bad number of particles"));
  if( f0==NULL ) ERROR(("Bad field"));
  if( g==NULL  ) ERROR(("Bad grid"));
  
  qdt_2mc  = 0.5*q_m*g->dt/g->cvac;
  qdt_4mc2 = 0.25*q_m*g->dt/(g->cvac*g->cvac);
  c = g->cvac;
  r8V = 0.125/(g->dx*g->dy*g->dz);
  mc_q = g->cvac/q_m;

  stride_10 = &h(1,0,0) - &h(0,0,0);
  stride_21 = &h(0,1,0) - &h(1,0,0);
  stride_43 = &h(0,0,1) - &h(1,1,0);

  for(;n;n--,p++) {
    /* Load the particle */
    dx = p->dx;
    dy = p->dy;
    dz = p->dz;
    ii = p->i;
    ux = p->ux;
    uy = p->uy;
    uz = p->uz;
    q  = p->q;
    
    /* Half advance E */
    f  = f0 + ii;
    ux += qdt_2mc*((f->ex+dy*f->dexdy) + dz*(f->dexdz+dy*f->d2exdydz));
    uy += qdt_2mc*((f->ey+dz*f->deydz) + dx*(f->deydx+dz*f->d2eydzdx));
    uz += qdt_2mc*((f->ez+dx*f->dezdx) + dy*(f->dezdy+dx*f->d2ezdxdy));

    /* Boris rotation - Interpolate B field */
    w5 = f->cbx + dx*f->dcbxdx;
    w6 = f->cby + dy*f->dcbydy;
    w7 = f->cbz + dz*f->dcbzdz;

    /* Boris rotation - curl scalars (0.5 in v0 for half rotate)
       and kinetic energy computation. Note: gamma-1 = |u|^2 / (gamma+1) is
       the numerically accurate way to compute gamma-1 */
    ke_mc = ux*ux + uy*uy + uz*uz; /* ke_mc = |u|^2 (invariant) */
    vz = sqrt(1+ke_mc);            /* vz = gamma    (invariant) */
    ke_mc *= c/(vz+1);             /* ke_mc = c|u|^2/(gamma+1) = c*(gamma-1) */
    vz = c/vz;                     /* vz = c/gamma */
    w0 = qdt_4mc2*vz;
    w1 = w5*w5 + w6*w6 + w7*w7;    /* |cB|^2 */
    w2 = w0*w0*w1;
    w3 = w0*(1+(1./3.)*w2*(1+0.4*w2));
    w4 = w3/(1 + w1*w3*w3); w4 += w4;

    /* Boris rotation - uprime */
    w0 = ux + w3*( uy*w7 - uz*w6 );
    w1 = uy + w3*( uz*w5 - ux*w7 );
    w2 = uz + w3*( ux*w6 - uy*w5 );

    /* Boris rotation - u */
    ux += w4*( w1*w7 - w2*w6 );
    uy += w4*( w2*w5 - w0*w7 );
    uz += w4*( w0*w6 - w1*w5 );

    /* Compute physical velocities */
    vx  = ux*vz;
    vy  = uy*vz;
    vz *= uz;

    /* Compute the trilinear coefficients */
    w0  = r8V*q;    /* w0 = q/(8V) = w/8               */
    dx *= w0;       /* dx = wx                         */
    w1  = w0+dx;    /* w1 = w/8 + wx/8 = (w/8)(1+x)    */
    w0 -= dx;       /* w0 = w/8 - wx/8 = (w/8)(1-x)    */
    w3  = 1+dy;     /* w3 = 1+y                        */
    w2  = w0*w3;    /* w2 = (w/8)(1-x)(1+y)            */
    w3 *= w1;       /* w3 = (w/8)(1+x)(1+y)            */
    dy  = 1-dy;     /* dy = 1-y                        */
    w0 *= dy;       /* w0 = (w/8)(1-x)(1-y)            */
    w1 *= dy;       /* w1 = (w/8)(1+x)(1-y)            */
    w7  = 1+dz;     /* w7 = 1+z                        */
    w4  = w0*w7;    /* w4 = (w/8)(1-x)(1-y)(1+z) *Done */
    w5  = w1*w7;    /* w5 = (w/8)(1+x)(1-y)(1+z) *Done */
    w6  = w2*w7;    /* w6 = (w/8)(1-x)(1+y)(1+z) *Done */
    w7 *= w3;       /* w7 = (w/8)(1+x)(1+y)(1+z) *Done */
    dz  = 1-dz;     /* dz = 1-z                        */
    w0 *= dz;       /* w0 = (w/8)(1-x)(1-y)(1-z) *Done */
    w1 *= dz;       /* w1 = (w/8)(1+x)(1-y)(1-z) *Done */
    w2 *= dz;       /* w2 = (w/8)(1-x)(1+y)(1-z) *Done */
    w3 *= dz;       /* w3 = (w/8)(1+x)(1+y)(1-z) *Done */

    /* Accumulate the hydro fields */
#   define ACCUM_HYDRO(wn)           \
    h->rho += wn;    /* wn = q/V */  \
    h->jx  += wn*vx;                 \
    h->jy  += wn*vy;                 \
    h->jz  += wn*vz;                 \
    wn *= mc_q;      /* wn = mc/V */ \
    dx = wn*ux;      /* dx = px/V */ \
    dy = wn*uy;                      \
    dz = wn*uz;                      \
    h->ke  += wn*ke_mc;              \
    h->px  += dx;                    \
    h->py  += dy;                    \
    h->pz  += dz;                    \
    h->txx += dx*vx;                 \
    h->tyy += dy*vy;                 \
    h->tzz += dz*vz;                 \
    h->tyz += dy*vz;                 \
    h->tzx += dz*vx;                 \
    h->txy += dx*vy

    h = h0 + ii;    ACCUM_HYDRO(w0); /* Cell i,j,k       */
    h += stride_10; ACCUM_HYDRO(w1); /* Cell i+1,j,k     */
    h += stride_21; ACCUM_HYDRO(w2); /* Cell i,j+1,k     */
    h += stride_10; ACCUM_HYDRO(w3); /* Cell i+1,j+1,k   */
    h += stride_43; ACCUM_HYDRO(w4); /* Cell i,j,k+1     */
    h += stride_10; ACCUM_HYDRO(w5); /* Cell i+1,j,k+1   */
    h += stride_21; ACCUM_HYDRO(w6); /* Cell i,j+1,k+1   */
    h += stride_10; ACCUM_HYDRO(w7); /* Cell i+1,j+1,k+1 */

#   undef ACCUM_HYDRO
  }
}
