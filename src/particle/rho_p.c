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

/* accumulate_rho_p adds the charge density associated with the supplied
   particle array to the rhof of the fields. Trilinear interpolation is used.
   rhof is known at the nodes at the same time as particle positions. No
   effort is made to fix up edges of the computational domain. All particles
   on the list must be inbounds. */

#define f(x,y,z) f0[INDEX_FORTRAN_3(x,y,z,0,g->nx+1,0,g->ny+1,0,g->nz+1)]

void accumulate_rho_p( field_t * ALIGNED f0,
                       const particle_t * ALIGNED p,
                       int n,
                       const grid_t * g ) {
  float w0, w1, w2, w3, w4, w5, w6, w7, t, r8V, *rho;
  int stride_10, stride_21, stride_43;

  if( f0==NULL ) { ERROR(("Bad field"));                return; }
  if( p==NULL  ) { ERROR(("Bad particle array"));       return; }
  if( n<0      ) { ERROR(("Bad number of particles"));  return; }
  if( g==NULL  ) { ERROR(("Bad grid"));                 return; }

  r8V = 0.125/(g->dx*g->dy*g->dz);

  stride_10 = &f(1,0,0).rhof - &f(0,0,0).rhof;
  stride_21 = &f(0,1,0).rhof - &f(1,0,0).rhof;
  stride_43 = &f(0,0,1).rhof - &f(1,1,0).rhof;

  for(;n;n--,p++) {
    /* Compute the trilinear coefficients */
    t   = p->dx;    /* t  = x                          */
    w0  = r8V*p->q; /* w0 = w/8                        */
    t  *= w0;       /* t  = wx/8                       */
    w1  = w0+t;     /* w1 = w/8 + wx/8 = (w/8)(1+x)    */
    w0 -= t;        /* w0 = w/8 - wx/8 = (w/8)(1-x)    */
    t   = p->dy;    /* t  = y                          */
    w3  = 1+t;      /* w3 = 1+y                        */
    w2  = w0*w3;    /* w2 = (w/8)(1-x)(1+y)            */
    w3 *= w1;       /* w3 = (w/8)(1+x)(1+y)            */
    t   = 1-t;      /* t  = 1-y                        */
    w0 *= t;        /* w0 = (w/8)(1-x)(1-y)            */
    w1 *= t;        /* w1 = (w/8)(1+x)(1-y)            */
    t   = p->dz;    /* t  = z                          */
    w7  = 1+t;      /* w7 = 1+z                        */
    w4  = w0*w7;    /* w4 = (w/8)(1-x)(1-y)(1+z) *Done */
    w5  = w1*w7;    /* w5 = (w/8)(1+x)(1-y)(1+z) *Done */
    w6  = w2*w7;    /* w6 = (w/8)(1-x)(1+y)(1+z) *Done */
    w7 *= w3;       /* w7 = (w/8)(1+x)(1+y)(1+z) *Done */
    t   = 1-t;      /* t  = 1-z                        */
    w0 *= t;        /* w0 = (w/8)(1-x)(1-y)(1-z) *Done */
    w1 *= t;        /* w1 = (w/8)(1+x)(1-y)(1-z) *Done */
    w2 *= t;        /* w2 = (w/8)(1-x)(1+y)(1-z) *Done */
    w3 *= t;        /* w3 = (w/8)(1+x)(1+y)(1-z) *Done */

    /* Accumulate the particle charge */
    rho = &f0[p->i].rhof; *rho += w0;
    rho += stride_10;     *rho += w1;
    rho += stride_21;     *rho += w2;
    rho += stride_10;     *rho += w3;
    rho += stride_43;     *rho += w4;
    rho += stride_10;     *rho += w5;
    rho += stride_21;     *rho += w6;
    rho += stride_10;     *rho += w7;
  }
}
