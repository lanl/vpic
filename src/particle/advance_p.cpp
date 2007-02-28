/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised and extended from earlier V4PIC versions
 *
 */

#include <math.h> /* For sqrt */
#include <particle.h>

static particle_mover_t *
advance_p_no_v4( particle_t * RESTRICT ALIGNED p0,
                 const int i0,
                 const int i1,
                 float qdt_2mc, /* q_m on entry */
                 particle_mover_t * RESTRICT ALIGNED pm,
                 int nm,
                 accumulator_t * RESTRICT ALIGNED a0,
                 const interpolator_t * RESTRICT ALIGNED f0,
                 const grid_t * RESTRICT g ) {
  int ii, n = i1 - i0 + 1;
  float dx, dy, dz, ux, uy, uz, q;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4, v5;
  const float cdt_dx     = g->cvac*g->dt/g->dx;
  const float cdt_dy     = g->cvac*g->dt/g->dy;
  const float cdt_dz     = g->cvac*g->dt/g->dz;
  const float one        = 1;
  const float one_third  = 1./3.;
  const float two_fifths = 0.4;
  const interpolator_t *f;
  float *a;
  particle_t * RESTRICT ALIGNED p = p0 + i0;

  qdt_2mc *= 0.5*g->dt/g->cvac;
    
  for(;n;n--,p++) {
    dx = p->dx;                              /* Load position */
    dy = p->dy;
    dz = p->dz;
    ii = p->i;
    f = f0 + ii;                             /* Interpolate E */
    hax = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                    dz*( f->dexdz + dy*f->d2exdydz ) );
    hay = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                    dx*( f->deydx + dz*f->d2eydzdx ) );
    haz = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                    dy*( f->dezdy + dx*f->d2ezdxdy ) );
    cbx = f->cbx + dx*f->dcbxdx;             /* Interpolate B */
    cby = f->cby + dy*f->dcbydy;
    cbz = f->cbz + dz*f->dcbzdz;
    ux = p->ux;                              /* Load momentum */
    uy = p->uy;
    uz = p->uz;
    q  = p->q;
    ux += hax;                               /* Half advance E */
    uy += hay;
    uz += haz;
    v0 = qdt_2mc/(float)sqrt(one + (ux*ux + uy*uy + uz*uz)); /* Boris - scalars */
    v1 = cbx*cbx + cby*cby + cbz*cbz;
    v2 = v0*v0*v1;
    v3 = v0*(one+one_third*v2*(one+two_fifths*v2));
    v4 = v3/(one + v1*v3*v3);
    v4 += v4;
    v0 = ux + v3*( uy*cbz - uz*cby );        /* Boris - uprime */
    v1 = uy + v3*( uz*cbx - ux*cbz );
    v2 = uz + v3*( ux*cby - uy*cbx );
    ux += v4*( v1*cbz - v2*cby );            /* Boris - rotation */
    uy += v4*( v2*cbx - v0*cbz );
    uz += v4*( v0*cby - v1*cbx );
    ux += hax;                               /* Half advance E */
    uy += hay;
    uz += haz;
    p->ux = ux;                              /* Store momentum */
    p->uy = uy;
    p->uz = uz;
    v0 = one/(float)sqrt(one + (ux*ux+ uy*uy + uz*uz)); /* Get norm displacement */
    ux *= cdt_dx;
    uy *= cdt_dy;
    uz *= cdt_dz;
    ux *= v0;
    uy *= v0;
    uz *= v0;
    v0 = dx + ux;                            /* Steak midpoint (inbnds) */
    v1 = dy + uy;
    v2 = dz + uz;
    v3 = v0 + ux;                            /* New position */
    v4 = v1 + uy;
    v5 = v2 + uz;
    if(  v3<=one && v4<=one && v5<=one &&          /* Check if inbnds */
        -v3<=one && -v4<=one && -v5<=one ) {


      /* Common case (inbnds)
         Note: accumulator values are 4 times the total physical charge that
         passed through the appropriate current quadrant in a time-step */

      p->dx = v3;                            /* Store new position */
      p->dy = v4;
      p->dz = v5;
      dx = v0;                               /* Streak mid */
      dy = v1;
      dz = v2;
      v5 = q*ux*uy*uz*one_third;             /* Compute correction */
      a = (float *)( a0 + ii );              /* Get accumulator */
#     define accumulate_j(X,Y,Z)                                  \
      v4  = q*u##X;   /* v2 = q ux                            */  \
      v1  = v4*d##Y;  /* v1 = q ux dy                         */  \
      v0  = v4-v1;    /* v0 = q ux (1-dy)                     */  \
      v1 += v4;       /* v1 = q ux (1+dy)                     */  \
      v4  = one+d##Z; /* v4 = 1+dz                            */  \
      v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */  \
      v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */  \
      v4  = one-d##Z; /* v4 = 1-dz                            */  \
      v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */  \
      v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */  \
      v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */  \
      v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */  \
      v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */  \
      v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */  \
      a[0] += v0;                                                 \
      a[1] += v1;                                                 \
      a[2] += v2;                                                 \
      a[3] += v3
      accumulate_j(x,y,z); a += 4;
      accumulate_j(y,z,x); a += 4;
      accumulate_j(z,x,y);
#     undef accumulate_j

    } else if( nm>0 ) {

      pm->dispx = ux;
      pm->dispy = uy;
      pm->dispz = uz;
      pm->i     = p - p0;
      if( move_p( p0, pm, a0, g ) ) pm++, nm--;

    }
  }

  if( nm==0 ) WARNING(("No particle movers remaining"));
  return pm;
}

#ifdef V4VERSION
#include CONCAT3(<,V4VERSION,>)
using namespace v4;

static particle_mover_t *
advance_p_v4( particle_t * RESTRICT ALIGNED p0,
              int i0,
              const int i1,
              const float q_m,
              particle_mover_t * RESTRICT ALIGNED pm,
              int nm,
              accumulator_t * RESTRICT ALIGNED a0,
              const interpolator_t * RESTRICT ALIGNED f0,
              const grid_t * RESTRICT g ) {
  v4float dx, dy, dz; v4int ii;
  v4float ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v0, v1, v2, v3, v4, v5;
  v4int outbnd;
  v4float qdt_2mc(0.5*q_m*g->dt/g->cvac);
  v4float cdt_dx(g->cvac*g->dt/g->dx);
  v4float cdt_dy(g->cvac*g->dt/g->dy);
  v4float cdt_dz(g->cvac*g->dt/g->dz);
  v4float one(1.);
  v4float one_third(1./3.);
  v4float two_fifths(0.4);
  v4float neg_one(-1.);
  v4float *vp0, *vp1, *vp2, *vp3;
  particle_t * RESTRICT ALIGNED p;
  int nq;
  
  p   = p0 + i0;              /* Pointer to first particle */
  nq  = ( i1 - i0 + 1 ) >> 2; /* Number of complete quads */
  i0 += nq << 2;              /* Start of straggler quad */


  for( ; nq; nq--, p+=4 ) {
    swizzle(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

    // Interpolate fields
    vp0 = (v4float *)(f0 + ii(0));
    vp1 = (v4float *)(f0 + ii(1));
    vp2 = (v4float *)(f0 + ii(2));
    vp3 = (v4float *)(f0 + ii(3));
    swizzle(vp0,vp1,vp2,vp3,hax,v0,v1,v2);
    hax = qdt_2mc*((hax+dy*v0)+dz*(v1+dy*v2));
    swizzle(vp0+1,vp1+1,vp2+1,vp3+1,hay,v3,v4,v5);
    hay = qdt_2mc*((hay+dz*v3)+dx*(v4+dz*v5));
    swizzle(vp0+2,vp1+2,vp2+2,vp3+2,haz,v0,v1,v2);
    haz = qdt_2mc*((haz+dx*v0)+dy*(v1+dx*v2));
    swizzle(vp0+3,vp1+3,vp2+3,vp3+3,cbx,v3,cby,v4); cbx += dx*v3;
    /**/                                            cby += dy*v4;
    half_swizzle(vp0+4,vp1+4,vp2+4,vp3+4,cbz,v5);   cbz += dz*v5;

    // Update momentum
    swizzle(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,q);
    ux += hax;
    uy += hay;
    uz += haz;
    v0 = qdt_2mc*rsqrt(one + ux*ux + uy*uy + uz*uz);
    v1 = cbx*cbx + cby*cby + cbz*cbz;
    v2 = v0*v0*v1;
    v3 = v0*(one+one_third*v2*(one+two_fifths*v2));
    v4 = v3*rcp(one+v1*v3*v3); v4 += v4;
    v0 = ux + v3*( uy*cbz - uz*cby );
    v1 = uy + v3*( uz*cbx - ux*cbz );
    v2 = uz + v3*( ux*cby - uy*cbx );
    ux += v4*( v1*cbz - v2*cby );
    uy += v4*( v2*cbx - v0*cbz );
    uz += v4*( v0*cby - v1*cbx );
    ux += hax;
    uy += hay;
    uz += haz;
    deswizzle(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux);
    
    // Update the position of inbnd particles
    v0 = rsqrt(one + ux*ux + uy*uy + uz*uz);
    ux *= cdt_dx;
    uy *= cdt_dy;
    uz *= cdt_dz;
    ux *= v0;
    uy *= v0;
    uz *= v0;     // ux,uy,uz are normalized displ (relative to cell size)
    v0 = dx + ux;
    v1 = dy + uy;
    v2 = dz + uz; // New particle midpoint
    v3 = v0 + ux;
    v4 = v1 + uy;
    v5 = v2 + uz; // New particle position
    outbnd = (v3>one) | (v3<neg_one) |
             (v4>one) | (v4<neg_one) |
             (v5>one) | (v5<neg_one);
    cmov(outbnd,dx,v3);  // Do not update outbnd particles
    cmov(outbnd,dy,v4);
    cmov(outbnd,dz,v5);
    deswizzle(v3,v4,v5,ii,&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx);
    
    // Accumulate current of inbnd particles
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step
    czero(outbnd,q);               // Do not accumulate outbnd particles
    dx = v0;                       // Streak midpoint (valid for inbnd only)
    dy = v1;
    dz = v2;
    v5 = q*ux*uy*uz*one_third;     // Charge conservation correction
    vp0 = (v4float *)(a0 + ii(0)); // Accumulator pointers
    vp1 = (v4float *)(a0 + ii(1));
    vp2 = (v4float *)(a0 + ii(2));
    vp3 = (v4float *)(a0 + ii(3));
#   define accumulate_j(X,Y,Z)                                      \
    v4  = q*u##X;   /* v4 = q ux                            */      \
    v1  = v4*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v4-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v4;       /* v1 = q ux (1+dy)                     */      \
    v4  = one+d##Z; /* v4 = 1+dz                            */      \
    v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v4  = one-d##Z; /* v4 = 1-dz                            */      \
    v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */      \
    transpose(v0,v1,v2,v3);                                         \
    load(vp0,v4); v4 += v0; store(v4,vp0); vp0++;                   \
    load(vp1,v4); v4 += v1; store(v4,vp1); vp1++;                   \
    load(vp2,v4); v4 += v2; store(v4,vp2); vp2++;                   \
    load(vp3,v4); v4 += v3; store(v4,vp3); vp3++
    accumulate_j(x,y,z);
    accumulate_j(y,z,x);
    accumulate_j(z,x,y);
#   undef accumulate_j

    // Update position and accumulate outbnd

#   define move_outbnd(N)                       \
    if( outbnd(N) && nm>0 ) {                   \
      pm->dispx = ux(N);                        \
      pm->dispy = uy(N);                        \
      pm->dispz = uz(N);                        \
      pm->i     = (p - p0) + N;                 \
      if( move_p( p0, pm, a0, g ) ) pm++, nm--; \
    }
    move_outbnd(0);
    move_outbnd(1);
    move_outbnd(2);
    move_outbnd(3);
#   undef move_outbnd

  }

  // Accumulate the final incomplete bundle
  // Note: advance_p_no_v4 will give a warning if ran out of movers
  return advance_p_no_v4( p0, i0, i1, q_m, pm, nm, a0, f0, g );
}

#endif

/* Returns the number of particle movers in use */
int advance_p( particle_t * RESTRICT ALIGNED p,
               const int n,
               const float q_m,
               particle_mover_t * RESTRICT ALIGNED pm,
               int nm,       
               accumulator_t * RESTRICT ALIGNED a,
               const interpolator_t * RESTRICT ALIGNED f,
               const grid_t * RESTRICT g ) {

  if( p==NULL  ) { ERROR(("Bad particle array"));      return 0; }
  if( n<0      ) { ERROR(("Bad number of particles")); return 0; }
  if( pm==NULL ) { ERROR(("Bad particle mover"));      return 0; }
  if( nm<0     ) { ERROR(("Bad number of movers"));    return 0; }
  if( a==NULL  ) { ERROR(("Bad accumulator"));         return 0; }
  if( f==NULL  ) { ERROR(("Bad interpolator"));        return 0; }
  if( g==NULL  ) { ERROR(("Bad grid"));                return 0; }

#ifdef V4VERSION  
  return advance_p_v4(    p, 0, n-1, q_m, pm, nm, a, f, g ) - pm;
#else
  return advance_p_no_v4( p, 0, n-1, q_m, pm, nm, a, f, g ) - pm;
#endif
}
