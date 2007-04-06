#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

double compute_rms_div_e_err( field_t * ALIGNED f,
                              const grid_t * g ) {
  compute_rms_div_e_err_pipeline_args_t args[1];
  int p;

  double err, local[2], global[2];
  int x, y, z, nx, ny, nz;
  field_t *f0;

  if( f==NULL ) ERROR(("Bad field"));
  if( g==NULL ) ERROR(("Bad grid"));
  
  /* Have the pipelines accumulate the interior of the local domain */

  args->f = f;
  args->g = g;
  dispatch_pipelines( compute_rms_div_e_err_pipeline, args, 0 );

  /* Have the host accumulator the exterior of the local domain */

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;

  err = 0;

  /* Do exterior faces */

  for( y=2; y<=ny; y++ ) {
    for( z=2; z<=nz; z++ ) {
      f0 = &f(   1, y, z); err += 0.5*(double)f0->div_e_err*(double)f0->div_e_err;
      f0 = &f(nx+1, y, z); err += 0.5*(double)f0->div_e_err*(double)f0->div_e_err;
    }
  }

  for( z=2; z<=nz; z++ ) {
    for( x=2; x<=nx; x++ ) {
      f0 = &f( x,   1, z); err += 0.5*(double)f0->div_e_err*(double)f0->div_e_err;
      f0 = &f( x,ny+1, z); err += 0.5*(double)f0->div_e_err*(double)f0->div_e_err;
    }
  }

  for( x=2; x<=nx; x++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(   x,   y,   1); err += 0.5*(double)f0->div_e_err*(double)f0->div_e_err;
      f0 = &f(   x,   y,nz+1); err += 0.5*(double)f0->div_e_err*(double)f0->div_e_err;
    }
  }

  /* Do exterior edges */

  for( x=2; x<=nx; x++ ) {
    f0 = &f(   x,   1,   1); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
    f0 = &f(   x,ny+1,   1); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
    f0 = &f(   x,   1,nz+1); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
    f0 = &f(   x,ny+1,nz+1); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
  }

  for( y=2; y<=ny; y++ ) {
    f0 = &f(   1,   y,   1); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
    f0 = &f(   1,   y,nz+1); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
    f0 = &f(nx+1,   y,   1); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
    f0 = &f(nx+1,   y,nz+1); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
  }

  for( z=2; z<=nz; z++ ) {
    f0 = &f(   1,   1,   z); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
    f0 = &f(nx+1,   1,   z); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
    f0 = &f(   1,nz+1,   z); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
    f0 = &f(nx+1,nz+1,   z); err += 0.25*(double)f0->div_e_err*(double)f0->div_e_err;
  }

  /* Do exterior corners */

  f0 = &f(   1,   1,   1); err += 0.125*(double)f0->div_e_err*(double)f0->div_e_err;
  f0 = &f(nx+1,   1,   1); err += 0.125*(double)f0->div_e_err*(double)f0->div_e_err;
  f0 = &f(   1,ny+1,   1); err += 0.125*(double)f0->div_e_err*(double)f0->div_e_err;
  f0 = &f(nx+1,ny+1,   1); err += 0.125*(double)f0->div_e_err*(double)f0->div_e_err;
  f0 = &f(   1,   1,nz+1); err += 0.125*(double)f0->div_e_err*(double)f0->div_e_err;
  f0 = &f(nx+1,   1,nz+1); err += 0.125*(double)f0->div_e_err*(double)f0->div_e_err;
  f0 = &f(   1,ny+1,nz+1); err += 0.125*(double)f0->div_e_err*(double)f0->div_e_err;
  f0 = &f(nx+1,ny+1,nz+1); err += 0.125*(double)f0->div_e_err*(double)f0->div_e_err;
  
  /* Reduce the results from the host and pipelines */

  wait_for_pipelines();

  for( p=0; p<n_pipeline; p++ ) err += args->err[p];

  /* Reduce the results from all nodes */

  local[0] = err*g->dx*g->dy*g->dz;
  local[1] = g->nx*g->ny*g->nz*g->dx*g->dy*g->dz;
  mp_allsum_d( local, global, 2, g->mp );
  return g->eps0*sqrt(global[0]/global[1]);
}

