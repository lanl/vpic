#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
 
/* WTF!  Under -ffast-math, gcc-4.1.1 thinks it is okay to treat the
   below as
     f0->cbx = ( f0->cbx + py*( blah ) ) - pz*( blah )
   even with explicit parenthesis are in there!  Oh my ...
   -ffast-math must not be used */

#define UPDATE_CBX() f0->cbx -= ( py*( fy->ez-f0->ez ) - pz*( fz->ey-f0->ey ) )
#define UPDATE_CBY() f0->cby -= ( pz*( fz->ex-f0->ex ) - px*( fx->ez-f0->ez ) )
#define UPDATE_CBZ() f0->cbz -= ( px*( fx->ey-f0->ey ) - py*( fy->ex-f0->ex ) )

void
advance_b( field_t * ALIGNED f,
           const grid_t * g,
           float frac ) {
  advance_b_pipeline_args_t args[1];
  
  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field")); 
  if( g==NULL ) ERROR(("Bad grid"));

  /* Do the bulk of the magnetic fields in the pipelines */
  
  args->f = f;
  args->g = g;
  args->frac = frac;
  dispatch_pipelines( advance_b_pipeline, args, 0 );
  
  /* While the pipelines are busy, do surface fields */
  
  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? frac*g->cvac*g->dt/g->dx : 0;
  py = (ny>1) ? frac*g->cvac*g->dt/g->dy : 0;
  pz = (nz>1) ? frac*g->cvac*g->dt/g->dz : 0;

  /* Do left over bx */
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fy = &f(nx+1,y+1,z);
      fz = &f(nx+1,y,  z+1);
      UPDATE_CBX();
    }
  }

  /* Do left over by */
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(2,ny+1,z);
    fz = &f(1,ny+1,z+1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_CBY();
      f0++;
      fx++;
      fz++;
    }
  }

  /* Do left over bz */
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,  nz+1);
    fx = &f(2,y,  nz+1);
    fy = &f(1,y+1,nz+1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_CBZ();
      f0++;
      fx++;
      fy++;
    }
  }

  local_adjust_norm_b(f,g);
  
  wait_for_pipelines();
}

