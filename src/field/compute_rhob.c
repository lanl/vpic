#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define UPDATE_RHO_B() f0->rhob = m[f0->nmat].nonconductive*      \
   ( px*( m[f0->ematx].epsx*f0->ex - m[fx->ematx].epsx*fx->ex ) + \
     py*( m[f0->ematy].epsy*f0->ey - m[fy->ematy].epsy*fy->ey ) + \
     pz*( m[f0->ematz].epsz*f0->ez - m[fz->ematz].epsz*fz->ez ) - \
     f0->rhof )

void
compute_rhob( field_t * ALIGNED f,
              const material_coefficient_t * ALIGNED m,
              const grid_t * g ) {
  compute_rhob_pipeline_args_t args[1];

  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field"));
  if( m==NULL ) ERROR(("Bad material coefficients"));
  if( g==NULL ) ERROR(("Bad grid"));
  
  /* Have the pipelines work on the interior of the local domain */
  /* FIXME: CHECK IF THIS CAN BE STARTED THIS EARLY */

  args->f = f;
  args->m = m;
  args->g = g;
  
  dispatch_pipelines( compute_rhob_pipeline, args, 0 );

  /* Have the host work on the exterior of the local domain */

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? g->eps0/g->dx : 0;
  py = (ny>1) ? g->eps0/g->dy : 0;
  pz = (nz>1) ? g->eps0/g->dz : 0;

  /* Begin setting normal e ghosts */

  begin_remote_ghost_norm_e( f, g );

  local_ghost_norm_e( f, g );

  /* Finish setting normal E ghosts */

  end_remote_ghost_norm_e( f, g );

  /* Compute divergence error in exterior */

  /* z faces, x edges, y edges and all corners */
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  1);
    fx = &f(0,y,  1);
    fy = &f(1,y-1,1);
    fz = &f(1,y,  0);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_RHO_B();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  nz+1);
    fx = &f(0,y,  nz+1);
    fy = &f(1,y-1,nz+1);
    fz = &f(1,y,  nz);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_RHO_B();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
 
  /* y faces, z edges */
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fx = &f(0,1,z);
    fy = &f(1,0,z);
    fz = &f(1,1,z-1);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_RHO_B();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(0,ny+1,z);
    fy = &f(1,ny,  z);
    fz = &f(1,ny+1,z-1);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_RHO_B();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  /* x faces */
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fx = &f(0,y,  z);
      fy = &f(1,y-1,z);
      fz = &f(1,y,  z-1);
      UPDATE_RHO_B();
      f0 = &f(nx+1,y,  z);
      fx = &f(nx,  y,  z);
      fy = &f(nx+1,y-1,z);
      fz = &f(nx+1,y,  z-1);
      UPDATE_RHO_B();
    }
  }

  /* Finish up setting the interior */
  /* FIXME: CHECK EXACTLY HOW LATE THIS CAN BE DONE */

  wait_for_pipelines();

  local_adjust_rhob(f,g);
}

