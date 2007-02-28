/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised and extended from earlier V4PIC versions
 *
 */

#include <particle.h>

/* move_p moves the particle m->p by m->dispx, m->dispy, m->dispz depositing
   particle current as it goes. If the particle was moved sucessfully
   (particle mover is no longer in use) returns 0. If the particle interacted
   with something this routine could not handle, this routine returns 1
   (particle mover is still in use). On a successful move, the particle
   position is updated and m->dispx, m->dispy and m->dispz are zerod. On a
   partial move, the particle position is updated to the point where the
   particle interacted and m->dispx, m->dispy, m->dispz contains the remaining
   particle displacement. The displacements are the physical displacments
   normalized current cell cell size.

   Because move_p is internal use only and frequently called, it does not
   check its input arguments. Higher level routines are responsible for
   insuring valid arguments. */

int move_p( particle_t       * RESTRICT ALIGNED p,
            particle_mover_t * RESTRICT ALIGNED pm,
            accumulator_t * RESTRICT ALIGNED a0,
            const grid_t * RESTRICT g ) {
  float s_midx, s_midy, s_midz;
  float s_dispx, s_dispy, s_dispz;
  float s_dir[3];
  float v0, v1, v2, v3, v4, v5;
  int type;
  INT64_TYPE neighbor;
  float *a;
  p += pm->i;

  for(;;) {
    s_midx = p->dx;
    s_midy = p->dy;
    s_midz = p->dz;

    s_dispx = pm->dispx;
    s_dispy = pm->dispy;
    s_dispz = pm->dispz;

    s_dir[0] = (s_dispx>0) ? 1 : -1;
    s_dir[1] = (s_dispy>0) ? 1 : -1;
    s_dir[2] = (s_dispz>0) ? 1 : -1;
    
    /* Compute the twice the fractional distance to each potential
       streak/cell face intersection. */
    v0 = (s_dispx==0) ? 3.4e38 : (s_dir[0]-s_midx)/s_dispx;
    v1 = (s_dispy==0) ? 3.4e38 : (s_dir[1]-s_midy)/s_dispy;
    v2 = (s_dispz==0) ? 3.4e38 : (s_dir[2]-s_midz)/s_dispz;

    /* Determine the fractional length and type of current streak. The streak
       ends on either the first face intersected by the particle track or at
       the end of the particle track.
         type 0,1 or 2 ... streak ends on a x,y or z-face respectively
         type 3        ... streak ends at end of the particle track */
    /**/      v3=2,  type=3;
    if(v0<v3) v3=v0, type=0;
    if(v1<v3) v3=v1, type=1;
    if(v2<v3) v3=v2, type=2;
    v3 *= 0.5;

    /* Compute the midpoint and the normalized displacement of the streak */
    s_dispx *= v3;
    s_dispy *= v3;
    s_dispz *= v3;
    s_midx += s_dispx;
    s_midy += s_dispy;
    s_midz += s_dispz;

    /* Accumulate the streak
       Note: accumulator values are 4 times the total physical charge that
       passed through the appropriate current quadrant in a time-step */
    v5 = p->q*s_dispx*s_dispy*s_dispz*(1./3.);
    a = (float *)(a0 + p->i);
#   define accumulate_j(X,Y,Z)                                        \
    v4  = p->q*s_disp##X; /* v2 = q ux                            */  \
    v1  = v4*s_mid##Y;    /* v1 = q ux dy                         */  \
    v0  = v4-v1;          /* v0 = q ux (1-dy)                     */  \
    v1 += v4;             /* v1 = q ux (1+dy)                     */  \
    v4  = 1+s_mid##Z;     /* v4 = 1+dz                            */  \
    v2  = v0*v4;          /* v2 = q ux (1-dy)(1+dz)               */  \
    v3  = v1*v4;          /* v3 = q ux (1+dy)(1+dz)               */  \
    v4  = 1-s_mid##Z;     /* v4 = 1-dz                            */  \
    v0 *= v4;             /* v0 = q ux (1-dy)(1-dz)               */  \
    v1 *= v4;             /* v1 = q ux (1+dy)(1-dz)               */  \
    v0 += v5;             /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */  \
    v1 -= v5;             /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */  \
    v2 -= v5;             /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */  \
    v3 += v5;             /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */  \
    a[0] += v0;                                                       \
    a[1] += v1;                                                       \
    a[2] += v2;                                                       \
    a[3] += v3
    accumulate_j(x,y,z); a += 4;
    accumulate_j(y,z,x); a += 4;
    accumulate_j(z,x,y);
#   undef accumulate_j

    /* Compute the remaining particle displacment */
    pm->dispx -= s_dispx;
    pm->dispy -= s_dispy;
    pm->dispz -= s_dispz;

    /* Compute the new particle offset */
    p->dx += s_dispx+s_dispx;
    p->dy += s_dispy+s_dispy;
    p->dz += s_dispz+s_dispz;

    /* If an end streak, return success (should be ~50% of the time) */
    if( type==3 ) return 0;

    /* Determine if the cell crossed into a local cell or if it hit a boundary
       Convert the coordinate system accordingly. Note: Crossing into a local
       cell should happen ~50% of the time; hitting a boundary is usually a
       rare event. Note: the entry / exit coordinate for the particle is
       guaranteed to be +/-1 _exactly_ for the particle. */
    v0 = s_dir[type];
    neighbor = g->neighbor[ 6*p->i + ((v0>0)?3:0) + type ];
    if( neighbor<g->rangel || neighbor>g->rangeh ) { /* Hit a boundary */
      (&(p->dx))[type] = v0;                         /* Put on boundary */
      if( neighbor!=reflect_particles ) return 1;    /* Cannot handle it */
      (&(p->ux))[type] = -(&(p->ux))[type];
      (&(pm->dispx))[type] = -(&(pm->dispx))[type];
    } else {
      p->i = neighbor - g->rangel; /* Compute local index of neighbor */
                                   /* Note: neighbor - g->rangel < 2^31 / 6 */
      (&(p->dx))[type] = -v0;      /* Convert coordinate system */
    }
  }
  return 0; /* Never get here ... avoid compiler warning */
}

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,g->nx+1,0,g->ny+1,0,g->nz+1)]

/* Removes a particle from the particle list. The removed particle is
   accumulated to rhob and swapped with the last particle on the particle
   last (backfilling). This routine returns 0 if the particle could not be
   removed and 1 if removed successfully.

   Even though remove_p is internal use only, because it is seldom called,
   remove_p lightly checks its input arguments to verify removability.

   Note: for particles injected from a neighboring domain, the removal of the
   particle from the neighboring domain will accumulate the particle rhob at
   the same position on this processor where the opposite charged particle in
   accumulated (see inject_p below). When rhob is synchronized, these charges
   cancel. */

int remove_p( particle_t * RESTRICT ALIGNED r,
              particle_t * RESTRICT ALIGNED p,
              int np,
              field_t * RESTRICT ALIGNED f,
              const grid_t * RESTRICT g ) {
  float w0, w1, w2, w3, w4, w5, w6, w7, t;
  int i, j, k;
  float *rhob;

  if( r==NULL || p==NULL || f==NULL || g==NULL || (r-p)<0 || (r-p)>=np )
    return 0;

  /* Backfill the particle to remove */
  np--;
  p += np;
  t = p->dx; p->dx = r->dx; r->dx = t;
  t = p->dy; p->dy = r->dy; r->dy = t;
  t = p->dz; p->dz = r->dz; r->dz = t;
  i = p->i;  p->i  = r->i;  r->i  = i;
  t = p->ux; p->ux = r->ux; r->ux = t;
  t = p->uy; p->uy = r->uy; r->uy = t;
  t = p->uz; p->uz = r->uz; r->uz = t;
  t = p->q;  p->q  = r->q;  r->q  = t;

  /* Compute the trilinear weights */
  t   = p->dx;                      /* t  = x                          */
  w0  = p->q/(8*g->dx*g->dy*g->dz); /* w0 = w/8                        */
  t  *= w0;                         /* t  = wx/8                       */
  w1  = w0+t;                       /* w1 = w/8 + wx/8 = (w/8)(1+x)    */
  w0 -= t;                          /* w0 = w/8 - wx/8 = (w/8)(1-x)    */
  t   = p->dy;                      /* t  = y                          */
  w3  = 1+t;                        /* w3 = 1+y                        */
  w2  = w0*w3;                      /* w2 = (w/8)(1-x)(1+y)            */
  w3 *= w1;                         /* w3 = (w/8)(1+x)(1+y)            */
  t   = 1-t;                        /* t  = 1-y                        */
  w0 *= t;                          /* w0 = (w/8)(1-x)(1-y)            */
  w1 *= t;                          /* w1 = (w/8)(1+x)(1-y)            */
  t   = p->dz;                      /* t  = z                          */
  w7  = 1+t;                        /* w7 = 1+z                        */
  w4  = w0*w7;                      /* w4 = (w/8)(1-x)(1-y)(1+z) *Done */
  w5  = w1*w7;                      /* w5 = (w/8)(1+x)(1-y)(1+z) *Done */
  w6  = w2*w7;                      /* w6 = (w/8)(1-x)(1+y)(1+z) *Done */
  w7 *= w3;                         /* w7 = (w/8)(1+x)(1+y)(1+z) *Done */
  t   = 1-t;                        /* t  = 1-z                        */
  w0 *= t;                          /* w0 = (w/8)(1-x)(1-y)(1-z) *Done */
  w1 *= t;                          /* w1 = (w/8)(1+x)(1-y)(1-z) *Done */
  w2 *= t;                          /* w2 = (w/8)(1-x)(1+y)(1-z) *Done */
  w3 *= t;                          /* w3 = (w/8)(1+x)(1+y)(1-z) *Done */
  
  /* Adjust the weights for a corrected local accumulation of rhob */
  i  = p->i;        /* i = INDEX_FORTRAN_3(ix,iy,iz,0,nx+1,0,ny+1,0,nz+1)
                         = ix + (nx+2)*( iy + (ny+2)*iz ) */
  j  = i/(g->nx+2); /* j = iy + (ny+2)*iz */
  i -= j*(g->nx+2); /* i = ix */
  k  = j/(g->ny+2); /* k = iz */
  j -= k*(g->ny+2); /* j = iy */
  if( i==1     ) w0 += w0, w2 += w2, w4 += w4, w6 += w6;
  if( i==g->nx ) w1 += w1, w3 += w3, w5 += w5, w7 += w7;
  if( j==1     ) w0 += w0, w1 += w1, w4 += w4, w5 += w5;
  if( j==g->ny ) w2 += w2, w3 += w3, w6 += w6, w7 += w7;
  if( k==1     ) w0 += w0, w1 += w1, w2 += w2, w3 += w3;
  if( k==g->nz ) w4 += w4, w5 += w5, w6 += w6, w7 += w7;
  
  /* Update rhob */
  i = &f(1,0,0).rhob - &f(0,0,0).rhob;
  j = &f(0,1,0).rhob - &f(1,0,0).rhob;
  k = &f(0,0,1).rhob - &f(1,1,0).rhob;
  rhob = &f[p->i].rhob; *rhob += w0;
  rhob += i;            *rhob += w1;
  rhob += j;            *rhob += w2;
  rhob += i;            *rhob += w3;
  rhob += k;            *rhob += w4;
  rhob += i;            *rhob += w5;
  rhob += j;            *rhob += w6;
  rhob += i;            *rhob += w7;

  return 1;
}

/* Inject a particle into the simulation. Inject first loads the particle
   position and mover according to the injector. Then it accumulates an equal
   and opposite charge to rhob (corrected local accumulation) to keep the
   charge conservation books balanced. The particle is then moved into its
   final position. Returns 0 if the supplied mover is no longer in use.
   Returns 1 if the supplied mover mover is in use still in use.

   Arguments are lightly tested for validity. If bad arguments are passed,
   an error message is printed and 0 (mover not in use) is returned.

   Note: for particles injected from a neighboring domain, the removal of the
   particle from the neighboring domain will accumulate the particle rhob at
   the same position on this processor where the opposite charged particle in
   accumulated (see remove_p above). When rhob is synchronized, these charges
   cancel. */

int inject_p( particle_t * RESTRICT ALIGNED p0,       /* Array to inject into */
              int n,                                  /* Where to inject */
              particle_mover_t * RESTRICT ALIGNED pm, /* Particle mover */
              field_t * RESTRICT ALIGNED f,
              accumulator_t * RESTRICT ALIGNED a,
              const particle_injector_t * RESTRICT pi,
              const grid_t * RESTRICT g ) {
  float w0, w1, w2, w3, w4, w5, w6, w7, t;
  int i, j, k;
  float *rhob;
  particle_t * p = p0 + n;

  if( p==NULL  ) { ERROR(("Bad particle"));    return -1; }
  if( pm==NULL ) { ERROR(("Bad mover"));       return -1; }
  if( f==NULL  ) { ERROR(("Bad field"));       return -1; }
  if( a==NULL  ) { ERROR(("Bad accumulator")); return -1; }
  if( pi==NULL ) { ERROR(("Bad injector"));    return -1; }
  if( g==NULL  ) { ERROR(("Bad grid"));        return -1; }

  /* Load the particle and particle mover from the injector */
  p->dx     = pi->dx;
  p->dy     = pi->dy;
  p->dz     = pi->dz;
  p->i      = pi->i;
  p->ux     = pi->ux;
  p->uy     = pi->uy;
  p->uz     = pi->uz;
  p->q      = pi->q;
  pm->dispx = pi->dispx;
  pm->dispy = pi->dispy;
  pm->dispz = pi->dispz;
  pm->i     = n;

  /* Compute the trilinear weights */
  t   = p->dx;                       /* t  = x                          */
  w0  = -p->q/(8*g->dx*g->dy*g->dz); /* w0 = w/8                        */
  t  *= w0;                          /* t  = wx/8                       */
  w1  = w0+t;                        /* w1 = w/8 + wx/8 = (w/8)(1+x)    */
  w0 -= t;                           /* w0 = w/8 - wx/8 = (w/8)(1-x)    */
  t   = p->dy;                       /* t  = y                          */
  w3  = 1+t;                         /* w3 = 1+y                        */
  w2  = w0*w3;                       /* w2 = (w/8)(1-x)(1+y)            */
  w3 *= w1;                          /* w3 = (w/8)(1+x)(1+y)            */
  t   = 1-t;                         /* t  = 1-y                        */
  w0 *= t;                           /* w0 = (w/8)(1-x)(1-y)            */
  w1 *= t;                           /* w1 = (w/8)(1+x)(1-y)            */
  t   = p->dz;                       /* t  = z                          */
  w7  = 1+t;                         /* w7 = 1+z                        */
  w4  = w0*w7;                       /* w4 = (w/8)(1-x)(1-y)(1+z) *Done */
  w5  = w1*w7;                       /* w5 = (w/8)(1+x)(1-y)(1+z) *Done */
  w6  = w2*w7;                       /* w6 = (w/8)(1-x)(1+y)(1+z) *Done */
  w7 *= w3;                          /* w7 = (w/8)(1+x)(1+y)(1+z) *Done */
  t   = 1-t;                         /* t  = 1-z                        */
  w0 *= t;                           /* w0 = (w/8)(1-x)(1-y)(1-z) *Done */
  w1 *= t;                           /* w1 = (w/8)(1+x)(1-y)(1-z) *Done */
  w2 *= t;                           /* w2 = (w/8)(1-x)(1+y)(1-z) *Done */
  w3 *= t;                           /* w3 = (w/8)(1+x)(1+y)(1-z) *Done */
  
  /* Adjust the weights for a corrected local accumulation of rhob */
  i  = p->i;        /* i = INDEX_FORTRAN_3(ix,iy,iz,0,nx+1,0,ny+1,0,nz+1)
                         = ix + (nx+2)*( iy + (ny+2)*iz ) */
  j  = i/(g->nx+2); /* j = iy + (ny+2)*iz */
  i -= j*(g->nx+2); /* i = ix */
  k  = j/(g->ny+2); /* k = iz */
  j -= k*(g->ny+2); /* j = iy */
  if( i==1     ) w0 += w0, w2 += w2, w4 += w4, w6 += w6;
  if( i==g->nx ) w1 += w1, w3 += w3, w5 += w5, w7 += w7;
  if( j==1     ) w0 += w0, w1 += w1, w4 += w4, w5 += w5;
  if( j==g->ny ) w2 += w2, w3 += w3, w6 += w6, w7 += w7;
  if( k==1     ) w0 += w0, w1 += w1, w2 += w2, w3 += w3;
  if( k==g->nz ) w4 += w4, w5 += w5, w6 += w6, w7 += w7;
  
  /* Update rhob */
  i = &f(1,0,0).rhob - &f(0,0,0).rhob;
  j = &f(0,1,0).rhob - &f(1,0,0).rhob;
  k = &f(0,0,1).rhob - &f(1,1,0).rhob;
  rhob = &f[p->i].rhob; *rhob += w0;
  rhob += i;            *rhob += w1;
  rhob += j;            *rhob += w2;
  rhob += i;            *rhob += w3;
  rhob += k;            *rhob += w4;
  rhob += i;            *rhob += w5;
  rhob += j;            *rhob += w6;
  rhob += i;            *rhob += w7;

  return move_p(p0,pm,a,g);
}
