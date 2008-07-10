#define IN_spa
#include "spa_private.h"

// move_p moves the particle m->p by m->dispx, m->dispy, m->dispz
// depositing particle current as it goes. If the particle was moved
// sucessfully (particle mover is no longer in use) returns 0. If the
// particle interacted with something this routine could not handle,
// this routine returns 1 (particle mover is still in use). On a
// successful move, the particle position is updated and m->dispx,
// m->dispy and m->dispz are zerod. On a partial move, the particle
// position is updated to the point where the particle interacted and
// m->dispx, m->dispy, m->dispz contains the remaining particle
// displacement. The displacements are the physical displacments
// normalized current cell size.
//
// Because move_p is internal use only and frequently called, it does
// not check its input arguments. Higher level routines are
// responsible for insuring valid arguments.

int
move_p( particle_t       * ALIGNED(128) p0,
        particle_mover_t * ALIGNED(16)  pm,
        accumulator_t    * ALIGNED(128) a0,
        const grid_t     *              g ) {
  float s_midx, s_midy, s_midz;
  float s_dispx, s_dispy, s_dispz;
  float s_dir[3];
  float v0, v1, v2, v3, v4, v5;
  int type;
  int64_t neighbor;
  float *a;
  particle_t * ALIGNED(32) p = p0 + pm->i;

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
    
    // Compute the twice the fractional distance to each potential
    // streak/cell face intersection.
    v0 = (s_dispx==0) ? 3.4e38 : (s_dir[0]-s_midx)/s_dispx;
    v1 = (s_dispy==0) ? 3.4e38 : (s_dir[1]-s_midy)/s_dispy;
    v2 = (s_dispz==0) ? 3.4e38 : (s_dir[2]-s_midz)/s_dispz;

    // Determine the fractional length and type of current streak. The
    // streak ends on either the first face intersected by the
    // particle track or at the end of the particle track.
    // 
    //   type 0,1 or 2 ... streak ends on a x,y or z-face respectively
    //   type 3        ... streak ends at end of the particle track
    /**/      v3=2,  type=3;
    if(v0<v3) v3=v0, type=0;
    if(v1<v3) v3=v1, type=1;
    if(v2<v3) v3=v2, type=2;
    v3 *= 0.5;

    // Compute the midpoint and the normalized displacement of the streak
    s_dispx *= v3;
    s_dispy *= v3;
    s_dispz *= v3;
    s_midx += s_dispx;
    s_midy += s_dispy;
    s_midz += s_dispz;

    // Accumulate the streak.  Note: accumulator values are 4 times
    // the total physical charge that passed through the appropriate
    // current quadrant in a time-step
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

    // Compute the remaining particle displacment
    pm->dispx -= s_dispx;
    pm->dispy -= s_dispy;
    pm->dispz -= s_dispz;

    // Compute the new particle offset
    p->dx += s_dispx+s_dispx;
    p->dy += s_dispy+s_dispy;
    p->dz += s_dispz+s_dispz;

    // If an end streak, return success (should be ~50% of the time)
    if( type==3 ) return 0;

    // Determine if the cell crossed into a local cell or if it hit a
    // boundary Convert the coordinate system accordingly.  Note:
    // Crossing into a local cell should happen ~50% of the time;
    // hitting a boundary is usually a rare event.  Note: the entry /
    // exit coordinate for the particle is guaranteed to be +/-1
    // _exactly_ for the particle.

    v0 = s_dir[type];
    neighbor = g->neighbor[ 6*p->i + ((v0>0)?3:0) + type ];
    if( neighbor<g->rangel || neighbor>g->rangeh ) { // Hit a boundary
      (&(p->dx))[type] = v0;                         // Put on boundary
      if( neighbor!=reflect_particles ) return 1;    // Cannot handle it
      (&(p->ux))[type] = -(&(p->ux))[type];
      (&(pm->dispx))[type] = -(&(pm->dispx))[type];
    } else {
      p->i = neighbor - g->rangel; // Compute local index of neighbor
      /**/                         // Note: neighbor - g->rangel < 2^31 / 6
      (&(p->dx))[type] = -v0;      // Convert coordinate system
    }
  }
  return 0; // Never get here ... avoid compiler warning
}
