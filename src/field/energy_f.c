/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version (uses algorithms similar to earlier
 *                    V4PIC versions)
 *
 */

#include <field.h>

#define field(x,y,z) f0[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void energy_f( double * RESTRICT global,
               const field_t * RESTRICT ALIGNED f0,
               const material_coefficient_t * RESTRICT ALIGNED m,
               const grid_t * RESTRICT ALIGNED g ) {
  double v0, v1, local[6];
  int x, y, z, nx, ny, nz;
  const field_t *f;

  if( global==NULL ) { ERROR(("Bad energy"));                return; }
  if( f0==NULL )     { ERROR(("Bad field"));                 return; }
  if( m==NULL )      { ERROR(("Bad material coefficients")); return; }
  if( g==NULL )      { ERROR(("Bad grid"));                  return; }

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  
  /* Energy in the ex field */
  local[0] = 0;
  for( z=1; z<=nz+1; z++ ) {
    v0 = ((z==1 || z==nz+1) ? 0.5 : 1.);
    for( y=1; y<=ny+1; y++ ) {
      v1 = ((y==1 || y==ny+1) ? 0.5 : 1.)*v0;
      f = &field(1,y,z);
      for( x=1; x<=nx; x++ ) {
        local[0] += m[f->ematx].epsx*f->ex*f->ex*v1;
        f++;
      }
    }
  }

  /* Energy in the ey field */
  local[1] = 0;
  for( z=1; z<=nz+1; z++ ) {
    v0 = ((z==1 || z==nz+1) ? 0.5 : 1.);
    for( y=1; y<=ny; y++ ) {
      f = &field(1,y,z);
      for( x=1; x<=nx+1; x++ ) {
        v1 = ((x==1 || x==nx+1) ? 0.5 : 1.)*v0;
        local[1] += m[f->ematy].epsy*f->ey*f->ey*v1;
        f++;
      }
    }
  }

  /* Energy in the ez field */
  local[2] = 0;
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny+1; y++ ) {
      v0 = ((y==1 || y==ny+1) ? 0.5 : 1.);
      f = &field(1,y,z);
      for( x=1; x<=nx+1; x++ ) {
        v1 = ((x==1 || x==nx+1) ? 0.5 : 1.)*v0;
        local[2] += m[f->ematz].epsz*f->ez*f->ez*v1;
        f++;
      }
    }
  }

  /* Energy in the bx field */
  local[3] = 0;
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f = &field(1,y,z);
      for( x=1; x<=nx+1; x++ ) {
        v0 = ((x==1 || x==nx+1) ? 0.5 : 1.);
        local[3] += m[f->fmatx].rmux*f->cbx*f->cbx*v0;
        f++;
      }
    }
  }

  /* Energy in the by field */
  local[4] = 0;
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny+1; y++ ) {
      v0 = ((y==1 || y==ny+1) ? 0.5 : 1.);
      f = &field(1,y,z);
      for( x=1; x<=nx; x++ ) {
        local[4] += m[f->fmaty].rmuy*f->cby*f->cby*v0;
        f++;
      }
    }
  }

  /* Energy in the bz field */
  local[5] = 0;
  for( z=1; z<=nz+1; z++ ) {
    v0 = ((z==1 || z==nz+1) ? 0.5 : 1.);
    for( y=1; y<=ny; y++ ) {
      f = &field(1,y,z);
      for( x=1; x<=nx; x++ ) {
        local[5] += m[f->fmatz].rmuz*f->cbz*f->cbz*v0;
        f++;
      }
    }
  }

  /* Convert to physical units */
  v0 = 0.5*g->eps0*g->dx*g->dy*g->dz;
  local[0] *= v0;
  local[1] *= v0;
  local[2] *= v0;
  local[3] *= v0;
  local[4] *= v0;
  local[5] *= v0;

  mp_allsum_d( local, global, 6, g->mp );
}
