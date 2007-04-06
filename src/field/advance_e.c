#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define UPDATE_EX()						            \
  f0->tcax = ( py*(f0->cbz*m[f0->fmatz].rmuz-fy->cbz*m[fy->fmatz].rmuz) -   \
               pz*(f0->cby*m[f0->fmaty].rmuy-fz->cby*m[fz->fmaty].rmuy) ) - \
    damp*f0->tcax;                                                          \
  f0->ex = m[f0->ematx].decayx*f0->ex +                                     \
           m[f0->ematx].drivex*( f0->tcax - cj*f0->jfx )
#define UPDATE_EY()						            \
  f0->tcay = ( pz*(f0->cbx*m[f0->fmatx].rmux-fz->cbx*m[fz->fmatx].rmux) -   \
               px*(f0->cbz*m[f0->fmatz].rmuz-fx->cbz*m[fx->fmatz].rmuz) ) - \
    damp*f0->tcay;                                                          \
  f0->ey = m[f0->ematy].decayy*f0->ey +                                     \
           m[f0->ematy].drivey*( f0->tcay - cj*f0->jfy )
#define UPDATE_EZ()						            \
  f0->tcaz = ( px*(f0->cby*m[f0->fmaty].rmuy-fx->cby*m[fx->fmaty].rmuy) -   \
               py*(f0->cbx*m[f0->fmatx].rmux-fy->cbx*m[fy->fmatx].rmux) ) - \
    damp*f0->tcaz;                                                          \
  f0->ez = m[f0->ematz].decayz*f0->ez +                                     \
           m[f0->ematz].drivez*( f0->tcaz - cj*f0->jfz )

void
advance_e( field_t * ALIGNED f,
           const material_coefficient_t * ALIGNED m,
           const grid_t * g ) {
  advance_e_pipeline_args_t args[1];
  
  float damp, px, py, pz, cj;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field"));
  if( m==NULL ) ERROR(("Bad material coefficients"));
  if( g==NULL ) ERROR(("Bad grid"));

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  damp = g->damp;
  px = (nx>1) ? (1+damp)*g->cvac*g->dt/g->dx : 0;
  py = (ny>1) ? (1+damp)*g->cvac*g->dt/g->dy : 0;
  pz = (nz>1) ? (1+damp)*g->cvac*g->dt/g->dz : 0;
  cj = g->dt/g->eps0;

  /***************************************************************************
   * Begin tangential B ghost setup
   ***************************************************************************/
  
  begin_remote_ghost_tang_b( f, g );
  local_ghost_tang_b( f, g );

  /***************************************************************************
   * Update interior fields
   * Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
   * Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
   * Note: ez all (1:nx+1,1:ny+1,1:nz  ) interior (1:nx,1:ny,2:nz)
   ***************************************************************************/

  /* Do majority interior e in a single pass.     
     FIXME: CHECK IF IT IS SAFE TO DISPATCH THE PIPELINES EVEN EARLIER AND
     COMPLETE THEM EVEN LATER.  I DON'T THINK IT IS SAFE TO DISPATCH
     THEM EARLIER UNDER ABSORBING BOUNDARY CONDITIONS BUT I AM NOT SURE.
     I THINK IT IS PROBABLY SAFE TO FINISH THEM LATER. */
     
  args->f = f;
  args->m = m;
  args->g = g;

  /* FIXME: FOR REASONS NOT YET UNDERSTOOD, THIS V4 PIPELINE IS ACTING
     UP ON OPTERON / 64-bit / GCC-4.1.2 BUT WORKS FINE on PENTIUM4 / 32-bit
     GCC-4.1.1 FOR OTHERWISE IDENTICAL COMPILE FLAGS. */
  dispatch_pipelines( advance_e_pipeline, args, 0 );
  
  /* Do left over interior ex */
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fy = &f(1,y-1,z);
      fz = &f(1,y,  z-1);
      UPDATE_EX();
    }
  }

  /* Do left over interior ey */
  for( z=2; z<=nz; z++ ) {
    f0 = &f(2,1,z);
    fx = &f(1,1,z);
    fz = &f(2,1,z-1);
    for( x=2; x<=nx; x++ ) {
      UPDATE_EY();
      f0++;
      fx++;
      fz++;
    }
  }

  /* Do left over interior ez */
  for( y=2; y<=ny; y++ ) {
    f0 = &f(2,y,  1);
    fx = &f(1,y,  1);
    fy = &f(2,y-1,1);
    for( x=2; x<=nx; x++ ) {
      UPDATE_EZ();
      f0++;
      fx++;
      fy++;
    }
  }

  wait_for_pipelines();
  
  /***************************************************************************
   * Finish tangential B ghost setup
   ***************************************************************************/

  end_remote_ghost_tang_b( f, g );

  /***************************************************************************
   * Update exterior fields
   ***************************************************************************/

  /* Do exterior ex */
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  1);
    fy = &f(1,y-1,1);
    fz = &f(1,y,  0);
    for( x=1; x<=nx; x++ ) {
      UPDATE_EX();
      f0++;
      fy++;
      fz++;
    }
  }
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  nz+1);
    fy = &f(1,y-1,nz+1);
    fz = &f(1,y,  nz);
    for( x=1; x<=nx; x++ ) {
      UPDATE_EX();
      f0++;
      fy++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fy = &f(1,0,z);
    fz = &f(1,1,z-1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_EX();
      f0++;
      fy++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fy = &f(1,ny,  z);
    fz = &f(1,ny+1,z-1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_EX();
      f0++;
      fy++;
      fz++;
    }
  }

  /* Do exterior ey */
  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(0,y,z);
      fz = &f(1,y,z-1);
      UPDATE_EY();
    }
  }
  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,z);
      fx = &f(nx,  y,z);
      fz = &f(nx+1,y,z-1);
      UPDATE_EY();
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(2,y,1);
    fx = &f(1,y,1);
    fz = &f(2,y,0);
    for( x=2; x<=nx; x++ ) {
      UPDATE_EY();
      f0++;
      fx++;
      fz++;
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(2,y,nz+1);
    fx = &f(1,y,nz+1);
    fz = &f(2,y,nz  );
    for( x=2; x<=nx; x++ ) {
      UPDATE_EY();
      f0++;
      fx++;
      fz++;
    }
  }

  /* Do exterior ez */
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fx = &f(0,1,z);
    fy = &f(1,0,z);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_EZ();
      f0++;
      fx++;
      fy++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(0,ny+1,z);
    fy = &f(1,ny,  z);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_EZ();
      f0++;
      fx++;
      fy++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fx = &f(0,y,  z);
      fy = &f(1,y-1,z);
      UPDATE_EZ();
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fx = &f(nx,  y,  z);
      fy = &f(nx+1,y-1,z);
      UPDATE_EZ();
    }
  }

  local_adjust_tang_e( f, g );
}
