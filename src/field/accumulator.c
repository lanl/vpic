/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version (uses algorithms heavily adapted from
 *                    earlier V4PIC versions)
 *
 */

#include <field.h>

/*****************************************************************************
 * Going into unload_accumulator, the accumulator contains 4 times the net
 * amount of charge that crossed the quarter face associated with each
 * accumulator component (has units of physical charge, i.e. C).
 * unload_accumulator computes the physical current density (A/m^2 in MKS
 * units) associated with all local quarter faces and accumulates the local
 * quarter faces to jf.
 *****************************************************************************/

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
#define a(x,y,z) a[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void unload_accumulator( field_t * RESTRICT ALIGNED f, 
                         const accumulator_t * RESTRICT ALIGNED a,
                         const grid_t * RESTRICT g ) {
  float cx, cy, cz;
  int x, y, z, nx, ny, nz;
  const float *pa;
  field_t *f0, *fx, *fy, *fz, *fyz, *fzx, *fxy;
  
  if( f==NULL ) { ERROR(("Bad field"));       return; }
  if( a==NULL ) { ERROR(("Bad accumulator")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));        return; }

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;

  cx = 0.25 / (g->dt*g->dy*g->dz);
  cy = 0.25 / (g->dt*g->dz*g->dx);
  cz = 0.25 / (g->dt*g->dx*g->dy);
  
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {

      pa  = &a(1,y,z).jx[0];
      f0  = &f(1,y,  z  );
      fx  = &f(2,y,  z  );
      fy  = &f(1,y+1,z  );
      fz  = &f(1,y,  z+1);
      fyz = &f(1,y+1,z+1);
      fzx = &f(2,y,  z+1);
      fxy = &f(2,y+1,z  );

      for( x=1; x<=nx; x++ ) {

        f0->jfx  += cx*pa[0];  /* f(x,y,  z  ).jfx += a(x,y,z).jx[0] */
        fy->jfx  += cx*pa[1];  /* f(x,y+1,z  ).jfx += a(x,y,z).jx[1] */
        fz->jfx  += cx*pa[2];  /* f(x,y,  z+1).jfx += a(x,y,z).jx[2] */
        fyz->jfx += cx*pa[3];  /* f(x,y+1,z+1).jfx += a(x,y,z).jx[3] */

        f0->jfy  += cy*pa[4];  /* f(x,  y,z  ).jfy += a(x,y,z).jy[0] */
        fz->jfy  += cy*pa[5];  /* f(x,  y,z+1).jfy += a(x,y,z).jy[1] */
        fx->jfy  += cy*pa[6];  /* f(x+1,y,z  ).jfy += a(x,y,z).jy[2] */
        fzx->jfy += cy*pa[7];  /* f(x+1,y,z+1).jfy += a(x,y,z).jy[3] */

        f0->jfz  += cz*pa[8];  /* f(x,  y,  z).jfz += a(x,y,z).jz[0] */
        fx->jfz  += cz*pa[9];  /* f(x+1,y,  z).jfz += a(x,y,z).jz[1] */
        fy->jfz  += cz*pa[10]; /* f(x,  y+1,z).jfz += a(x,y,z).jz[2] */
        fxy->jfz += cz*pa[11]; /* f(x+1,y+1,z).jfz += a(x,y,z).jz[3] */

        pa += 12;
        f0++;
        fx++;
        fy++;
        fz++;
        fyz++;
        fzx++;
        fxy++;
      }
    }
  }
}
