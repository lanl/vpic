/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version (uses similar algorithm to earlier
 *                    V4PIC versions)
 */

/* See notes in div_e.c */

#include <math.h> /* For sqrt */
#include <field.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
	     
/*****************************************************************************
 * compute_div_b_err applies the following difference equation:
 *   div_b_err = div cB
 *****************************************************************************/

void compute_div_b_err( field_t * RESTRICT ALIGNED f,
                        const grid_t * RESTRICT g ) {
  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) { ERROR(("Bad field")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));  return; }

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;

  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(2,y,z);
      fy = &f(1,y+1,z);
      fz = &f(1,y,z+1);
      for( x=1; x<=nx; x++ ) {
	f0->div_b_err = px*( fx->cbx - f0->cbx ) +
	                py*( fy->cby - f0->cby ) +
                        pz*( fz->cbz - f0->cbz );
	f0++;
	fx++;
	fy++;
	fz++;
      }
    }
  }
}

/*****************************************************************************
 * compute_rms_div_b_err returns
 *   eps0 sqrt( Integral |div_b_err|^2 / Volume )
 * This has units of (electric charge / volume) ... yes electric charge. The
 * integrals are done over all the domains. The volume is the total volume
 * of all domains. Every processor gets the same value.
 *****************************************************************************/

double compute_rms_div_b_err( field_t * RESTRICT ALIGNED f,
                              const grid_t * RESTRICT g ) {
  double err, local[2], global[2];
  field_t *f0;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) { ERROR(("Bad field")); return -1; }
  if( g==NULL ) { ERROR(("Bad grid"));  return -1; }

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;

  err = 0;
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      for( x=1; x<=nx; x++ ) {
        err += f0->div_b_err*f0->div_b_err;
        f0++;
      }
    }
  }

  local[0] = err*g->dx*g->dy*g->dz;
  local[1] = g->nx*g->ny*g->nz*g->dx*g->dy*g->dz;
  mp_allsum_d( local, global, 2, g->mp );
  return g->eps0*sqrt(global[0]/global[1]);
}
  
/*****************************************************************************
 * clean_div_b applies the following difference equation:
 *   cB_new = cB_old + alpha dt grad div_b_err
 * alpha is picked to rapidly reduce the rms_div_b_err
 *****************************************************************************/

void clean_div_b( field_t * RESTRICT ALIGNED f,
		  const grid_t * RESTRICT g ) {
  float alphadt, px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) { ERROR(("Bad field")); return; }
  if( g==NULL ) { ERROR(("Bad grid"));  return; }

# define MARDER_CBX() f0->cbx += px*( f0->div_b_err - fx->div_b_err )
# define MARDER_CBY() f0->cby += py*( f0->div_b_err - fy->div_b_err )
# define MARDER_CBZ() f0->cbz += pz*( f0->div_b_err - fz->div_b_err )

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;
  alphadt = 0.3888889/( px*px + py*py + pz*pz );
  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;

  /* Begin setting derr ghosts */
  begin_remote_ghost_div_b( f, g );
  local_ghost_div_b( f, g);

  /* Do Marder pass in interior */

  /* Do majority in a single pass */
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
	MARDER_CBX();
	MARDER_CBY();
	MARDER_CBZ();
	f0++;
	fx++;
	fy++;
	fz++;
      }
    }
  }

  /* Do left over bx */
  for( y=1; y<=ny; y++ ) {
    f0 = &f(2,y,1);
    fx = &f(1,y,1);
    for( x=2; x<=nx; x++ ) {
      MARDER_CBX();
      f0++;
      fx++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(2,1,z);
    fx = &f(1,1,z);
    for( x=2; x<=nx; x++ ) {
      MARDER_CBX();
      f0++;
      fx++;
    }
  }

  /* Left over by */
  for( z=1; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fy = &f(1,y-1,z);
      MARDER_CBY();
    }
  }
  for( y=2; y<=ny; y++ ) {
    f0 = &f(2,y,  1);
    fy = &f(2,y-1,1);
    for( x=2; x<=nx; x++ ) {
      MARDER_CBY();
      f0++;
      fy++;
    }
  }

  /* Left over bz */
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fz = &f(1,1,z-1);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBZ();
      f0++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fz = &f(1,y,z-1);
      MARDER_CBZ();
    }
  }

  /* Finish setting derr ghosts */
  end_remote_ghost_div_b( f, g );

  /* Do Marder pass in exterior */

  /* Exterior bx */
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(0,y,z);
      MARDER_CBX();
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,z);
      fx = &f(nx,  y,z);
      MARDER_CBX();
    }
  }

  /* Exterior by */
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fy = &f(1,0,z);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBY();
      f0++;
      fy++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fy = &f(1,ny,  z);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBY();
      f0++;
      fy++;
    }
  }

  /* Exterior bz */
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,1);
    fz = &f(1,y,0);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBZ();
      f0++;
      fz++;
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,nz+1);
    fz = &f(1,y,nz);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBZ();
      f0++;
      fz++;
    }
  }

  local_adjust_norm_b(f,g);
}
