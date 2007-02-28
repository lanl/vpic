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
#include <species.h>

/* center_p half pushes the particle momentum to time center it
   All particles must be inbounds */

static void center_p_no_v4( particle_t * RESTRICT ALIGNED p,
                            int n,
                            const float q_m,
                            const interpolator_t * RESTRICT ALIGNED f0,
                            const grid_t * RESTRICT g ) {
  const interpolator_t *f;
  float ux, uy, uz;
  float hax, hay, haz;
  float cbx, cby, cbz;
  float v0, v1, v2, v3, v4;
  float qdt_2mc = 0.5*q_m*g->dt/g->cvac;
  
  for(;n;n--,p++) {
    /* Load the particle position */
    v0 = p->dx;
    v1 = p->dy;
    v2 = p->dz;
    f  = f0 + p->i;
    /* Interpolate the fields */
    hax = qdt_2mc*((f->ex+v1*f->dexdy) + v2*(f->dexdz+v1*f->d2exdydz));
    hay = qdt_2mc*((f->ey+v2*f->deydz) + v0*(f->deydx+v2*f->d2eydzdx));
    haz = qdt_2mc*((f->ez+v0*f->dezdx) + v1*(f->dezdy+v0*f->d2ezdxdy));
    cbx = f->cbx + v0*f->dcbxdx;
    cby = f->cby + v1*f->dcbydy;
    cbz = f->cbz + v2*f->dcbzdz;
    /* Load the particle momentum */
    ux = p->ux;
    uy = p->uy;
    uz = p->uz;
    /* Half E advance */
    ux += hax;
    uy += hay;
    uz += haz;
    /* Boris rotation - curl scalars (0.5 in v0 for half rotate) */
    v0 = 0.5*qdt_2mc/sqrt(1 + ux*ux + uy*uy + uz*uz);
    v1 = cbx*cbx + cby*cby + cbz*cbz;
    v2 = v0*v0*v1;
    v3 = v0*(1+(1./3.)*v2*(1+0.4*v2));
    v4 = v3/(1 + v1*v3*v3); v4 += v4;
    /* Boris rotation - uprime */
    v0 = ux + v3*( uy*cbz - uz*cby );
    v1 = uy + v3*( uz*cbx - ux*cbz );
    v2 = uz + v3*( ux*cby - uy*cbx );
    /* Boris rotation - u */
    ux += v4*( v1*cbz - v2*cby );
    uy += v4*( v2*cbx - v0*cbz );
    uz += v4*( v0*cby - v1*cbx );
    /* Store updated momentum */
    p->ux = ux;
    p->uy = uy;
    p->uz = uz;
  }
}

#ifdef V4VERSION
#include CONCAT3(<,V4VERSION,>)
using namespace v4;

static void center_p_v4( particle_t * RESTRICT ALIGNED p,
                         const int n,
                         const float q_m,
                         const interpolator_t * RESTRICT ALIGNED f0,
                         const grid_t * RESTRICT g ) {
  v4float dx, dy, dz; v4int ii;
  v4float ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v0, v1, v2, v3, v4, v5;
  v4float qdt_2mc(0.5*q_m*g->dt/g->cvac);
  v4float qdt_4mc(0.25*q_m*g->dt/g->cvac);
  v4float one(1.);
  v4float one_third(1./3.);
  v4float two_fifths(0.4);
  v4float *vp0, *vp1, *vp2, *vp3;
  particle_t *p_stop;

  for( p_stop=p+(n&(~3)); p<p_stop; p+=4 ) {

    // Interpolate the fields
    swizzle(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);
    vp0 = (v4float *)(f0 + ii(0));
    vp1 = (v4float *)(f0 + ii(1));
    vp2 = (v4float *)(f0 + ii(2));
    vp3 = (v4float *)(f0 + ii(3));
    swizzle(vp0++,vp1++,vp2++,vp3++,hax,v0,v1,v2); hax = qdt_2mc*((hax+dy*v0)+dz*(v1+dy*v2));
    swizzle(vp0++,vp1++,vp2++,vp3++,hay,v3,v4,v5); hay = qdt_2mc*((hay+dz*v3)+dx*(v4+dz*v5));
    swizzle(vp0++,vp1++,vp2++,vp3++,haz,v0,v1,v2); haz = qdt_2mc*((haz+dx*v0)+dy*(v1+dx*v2));
    swizzle(vp0++,vp1++,vp2++,vp3++,cbx,v3,cby,v4); cbx += dx*v3;
    /**/                                            cby += dy*v4;
    half_swizzle(vp0,vp1,vp2,vp3,cbz,v0);           cbz += dz*v0;

    // Update the momentum
    swizzle(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,q);
    ux += hax;
    uy += hay;
    uz += haz;
    v0 = qdt_4mc*rsqrt(one + ux*ux + uy*uy + uz*uz);
    v1 = cbx*cbx + cby*cby + cbz*cbz;
    v2 = v0*v0*v1;
    v3 = v0*(one+one_third*v2*(one+two_fifths*v2));
    v4 = v3*rcp(one+v1*v3*v3);
    v4 += v4;
    v0 = ux + v3*( uy*cbz - uz*cby );
    v1 = uy + v3*( uz*cbx - ux*cbz );
    v2 = uz + v3*( ux*cby - uy*cbx );
    ux += v4*( v1*cbz - v2*cby );
    uy += v4*( v2*cbx - v0*cbz );
    uz += v4*( v0*cby - v1*cbx );
    deswizzle(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux);
  }
  center_p_no_v4( p, n&3, q_m, f0, g );
}

#endif

void center_p( particle_t * RESTRICT ALIGNED p,
               const int n,
               const float q_m,
               const interpolator_t * RESTRICT ALIGNED f,
               const grid_t * RESTRICT g ) {

  if( p==NULL ) { ERROR(("Bad particle array"));      return; }
  if( n<0     ) { ERROR(("Bad number of particles")); return; }
  if( f==NULL ) { ERROR(("Bad interpolator"));        return; }
  if( g==NULL ) { ERROR(("Bad grid"));                return; }

#ifdef V4VERSION
  center_p_v4( p, n, q_m, f, g );
#else
  center_p_no_v4( p, n, q_m, f, g );
#endif
}
