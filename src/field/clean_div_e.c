#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define MARDER_EX() \
    f0->ex += m[f0->ematx].drivex*px*(fx->div_e_err-f0->div_e_err)

#define MARDER_EY() \
    f0->ey += m[f0->ematy].drivey*py*(fy->div_e_err-f0->div_e_err)

#define MARDER_EZ() \
    f0->ez += m[f0->ematz].drivez*pz*(fz->div_e_err-f0->div_e_err)

void clean_div_e( field_t * ALIGNED f,
		  const material_coefficient_t * ALIGNED m,
		  const grid_t * g ) {
  clean_div_e_pipeline_args_t args[1];

  float alphadt, px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field"));
  if( m==NULL ) ERROR(("Bad material coefficients"));
  if( g==NULL ) ERROR(("Bad grid"));

  /* Do majority of field components in single pass on the pipelines */

  args->f = f;
  args->m = m;
  args->g = g;

  dispatch_pipelines( clean_div_e_pipeline_v4, args, 0 );
  
  /* Do left over field components on the host */

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
  
  /* Do left over ex */
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,nz+1);
    fx = &f(2,y,nz+1);
    for( x=1; x<=nx; x++ ) {
      MARDER_EX();
      f0++; fx++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(2,ny+1,z);
    for( x=1; x<=nx; x++ ) {
      MARDER_EX();
      f0++; fx++;
    }
  }

  /* Do left over ey */
  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fy = &f(nx+1,y+1,z);
      MARDER_EY();
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,  nz+1);
    fy = &f(1,y+1,nz+1);
    for( x=1; x<=nx; x++ ) {
      MARDER_EY();
      f0++; fy++;
    }
  }

  /* Do left over ez */
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fz = &f(1,ny+1,z+1);
    for( x=1; x<=nx+1; x++ ) {
      MARDER_EZ();
      f0++; fz++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,z);
      fz = &f(nx+1,y,z+1);
      MARDER_EZ();
    }
  }

  wait_for_pipelines(); /* FIXME: FINSIH EVEN LATER?? */

  local_adjust_tang_e(f,g);
}

