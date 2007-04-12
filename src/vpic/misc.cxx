/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <vpic.hxx>

void vpic_simulation::IUO_allocate_fields(void) {
  int nx1 = grid->nx + 1, ny1 = grid->ny+1, nz1 = grid->nz+1;
  field = new_field(grid);
  interpolator = new_interpolator(grid);
  accumulator = new_accumulators(grid);
  hydro = new_hydro(grid);
  // Pre-size communications buffers
  // This is done to get 99.9% of memory allocation over with before
  // the simulation starts running
  mp_size_recv_buffer(BOUNDARY(-1, 0, 0),ny1*nz1*sizeof(hydro_t),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 1, 0, 0),ny1*nz1*sizeof(hydro_t),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 0,-1, 0),nz1*nx1*sizeof(hydro_t),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 0, 1, 0),nz1*nx1*sizeof(hydro_t),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 0, 0,-1),nx1*ny1*sizeof(hydro_t),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 0, 0, 1),nx1*ny1*sizeof(hydro_t),grid->mp);
  mp_size_send_buffer(BOUNDARY(-1, 0, 0),ny1*nz1*sizeof(hydro_t),grid->mp);
  mp_size_send_buffer(BOUNDARY( 1, 0, 0),ny1*nz1*sizeof(hydro_t),grid->mp);
  mp_size_send_buffer(BOUNDARY( 0,-1, 0),nz1*nx1*sizeof(hydro_t),grid->mp);
  mp_size_send_buffer(BOUNDARY( 0, 1, 0),nz1*nx1*sizeof(hydro_t),grid->mp);
  mp_size_send_buffer(BOUNDARY( 0, 0,-1),nx1*ny1*sizeof(hydro_t),grid->mp);
  mp_size_send_buffer(BOUNDARY( 0, 0, 1),nx1*ny1*sizeof(hydro_t),grid->mp);
}

error_code vpic_simulation::inject_particle( species_id id,
                                             double x,  double y,  double z,
                                             double ux, double uy, double uz,
                                             double q,  double age ) {
  species_t *sp;
  int ix, iy, iz;
  particle_injector_t pi;

  // Check input parameters
  if( age<0 || age>1 ) return ERROR_CODE("Invalid age");
  if( grid==NULL ) return ERROR_CODE("NULL grid");
  if( accumulator==NULL ) return ERROR_CODE("NULL accumulator");
  if( species_list==NULL ) return ERROR_CODE("No species defined");
  if( id<0 || id>species_list->id ) return ERROR_CODE("Invalid species ID");
  if( species_lookup==NULL ) return ERROR_CODE("Species list not finalized");
  sp = species_lookup[id];
  if( sp->p==NULL  || sp->np>=sp->max_np ||
      sp->pm==NULL || sp->nm>=sp->max_nm )
    return ERROR_CODE("No room to inject");
  
  // Compute the injection cell and coordinate in cell coordinate system
  // BJA:  Note the use of double precision here for accurate particle 
  //       placement on large meshes. 
  x -= (double)grid->x0;
  y -= (double)grid->y0;
  z -= (double)grid->z0;
  x /= (double)grid->dx;
  y /= (double)grid->dy;
  z /= (double)grid->dz;
  ix = (int)floor(x);
  iy = (int)floor(y);
  iz = (int)floor(z); 
  x -= ix;
  y -= iy;
  z -= iz;
  x += x-1;
  y += y-1;
  z += z-1;
  ix++;
  iy++;
  iz++;
    
  // Allow injection on the far walls of the local computational domain
  if( x<-1+1e-6 && ix==grid->nx+1 ) x=1, ix--;
  if( y<-1+1e-6 && iy==grid->ny+1 ) y=1, iy--;
  if( z<-1+1e-6 && iz==grid->nz+1 ) z=1, iz--;

  // Check if requested injection point is inbounds
  if( ix<1 || ix>grid->nx || iy<1 || iy>grid->ny || iz<1 || iz>grid->nz )
    return ERROR_CODE("Injection point is outside local domain");
  
  // Load the particle injector
  pi.dx = x;
  pi.dy = y;
  pi.dz = z;
  pi.i  = INDEX_FORTRAN_3(ix,iy,iz,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
  pi.ux = ux;
  pi.uy = uy;
  pi.uz = uz;
  pi.q  = q;

  age *= grid->cvac*grid->dt/sqrt( ux*ux + uy*uy + uz*uz + 1 );
  pi.dispx = ux*age/grid->dx;
  pi.dispy = uy*age/grid->dy;
  pi.dispz = uz*age/grid->dz;

  // Inject the particle
  sp->nm += inject_p( sp->p, sp->np, sp->pm+sp->nm, field, accumulator,
                      &pi, grid );
  sp->np++;
  // FIXME: if inject_p fails, sp->np should not be incremented

  return NO_ERROR;
}
