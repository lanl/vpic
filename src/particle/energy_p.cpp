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

// energy_p computes the kinetic energy (excluding rest mass) stored in the
// local particle list for the given species. The energy is computed at the
// same time level as the particle position. Since Boris rotation leaves the
// energy unchanged, only have to half push E to get the timecentered energy.
//
// Calculation is based on:
//    gamma^2 = 1 + |u|^2
// -> gamma^2-1 = |u|^2
// -> (gamma-1)(gamma+1) = |u|^2
// -> (gamma-1) = |u|^2 / (gamma+1)
// This is a numerically safe way to compute gamma-1. It is accurate
// both the non-relativistic and extreme relativistic limits.

static double
local_energy_p_no_v4( const particle_t * RESTRICT ALIGNED p,
                      int n,
                      float q_m,
                      const interpolator_t * RESTRICT ALIGNED f0,
                      const grid_t * RESTRICT g ) {
  double en;
  float v0, v1, v2;
  const float qdt_2mc = 0.5*q_m*g->dt/g->cvac;
  const interpolator_t *f;
  
  en = 0;
  for(;n;n--,p++) {
    f  = f0 + p->i;
    v0 = p->ux + qdt_2mc*((f->ex+p->dy*f->dexdy) + p->dz*(f->dexdz+p->dy*f->d2exdydz));
    v1 = p->uy + qdt_2mc*((f->ey+p->dz*f->deydz) + p->dx*(f->deydx+p->dz*f->d2eydzdx));
    v2 = p->uz + qdt_2mc*((f->ez+p->dx*f->dezdx) + p->dy*(f->dezdy+p->dx*f->d2ezdxdy));
    v0 = v0*v0 + v1*v1 + v2*v2;
    v0 /= sqrt(1+v0)+1;
    en += v0*p->q;
  }
  return g->cvac*g->cvac*en/q_m;
}

#ifdef V4VERSION
#include CONCAT3(<,V4VERSION,>)
using namespace v4;
                       
static double local_energy_p_v4( const particle_t * RESTRICT ALIGNED p,
                                 int n,
                                 float q_m,
                                 const interpolator_t * RESTRICT ALIGNED f0,
                                 const grid_t * RESTRICT g ) {
  double en0, en1, en2, en3;
  v4int ii;
  v4float dx, dy, dz, ux, uy, uz, w;
  v4float hax, hay, haz;
  v4float v0, v1, v2;
  v4float qdt_2mc(0.5*q_m*g->dt/g->cvac);
  v4float one(1.);
  v4float *vp0, *vp1, *vp2, *vp3;
  int nextra = n&3;

  en0 = en1 = en2 = en3 = 0;
  for(n>>=2;n;n--,p+=4) {
      // Interpolate the fields
      swizzle(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);
      vp0 = (v4float *)(f0 + ii(0));
      vp1 = (v4float *)(f0 + ii(1));
      vp2 = (v4float *)(f0 + ii(2));
      vp3 = (v4float *)(f0 + ii(3));
      swizzle(vp0++,vp1++,vp2++,vp3++,hax,v0,v1,v2); hax = qdt_2mc*((hax+dy*v0)+dz*(v1+dy*v2));
      swizzle(vp0++,vp1++,vp2++,vp3++,hay,v0,v1,v2); hay = qdt_2mc*((hay+dz*v0)+dx*(v1+dz*v2));
      swizzle(vp0++,vp1++,vp2++,vp3++,haz,v0,v1,v2); haz = qdt_2mc*((haz+dx*v0)+dy*(v1+dx*v2));
      
      // Compute gamma-1 for the particles
      swizzle(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,w);
      ux += hax;
      uy += hay;
      uz += haz;
      v0 = ux*ux + uy*uy + uz*uz;
      v0 /= sqrt(one+v0)+one;
      v0 *= w;
      
      // Accumulate the energies (done in double precision to reduce round-off)
      en0 += v0(0);
      en1 += v0(1);
      en2 += v0(2);
      en3 += v0(3);
    }
    return g->cvac*g->cvac*(en0 + en1 + en2 + en3)/q_m +
      local_energy_p_no_v4(p,nextra,q_m,f0,g);
}

#endif

double energy_p( const particle_t * RESTRICT ALIGNED p,
                 int n,
                 float q_m,
                 const interpolator_t * RESTRICT ALIGNED f,
                 const grid_t * RESTRICT g ) {
  double local, global;

  if( p==NULL ) { ERROR(("Bad particle array"));      return 0; }
  if( n<0     ) { ERROR(("Bad number of particles")); return 0; }
  if( f==NULL ) { ERROR(("Bad interpolator"));        return 0; }
  if( g==NULL ) { ERROR(("Bad grid"));                return 0; }

#ifdef V4VERSION
  local = local_energy_p_v4( p, n, q_m, f, g );
#else
  local = local_energy_p_no_v4( p, n, q_m, f, g );
#endif

  mp_allsum_d( &local, &global, 1, g->mp );
  return global;
}
