#include "boundary.h"
#include <stdio.h> // for debugging output

// Special absorbing boundary condition that generates a "link" output
// comprising data for particles absorbed at boundary.
// 
// Output is ASCII written to files link.XX where XX is the processor
// id.  Files contain the following:
// 
// species_id, x, y, z, ux, uy, uz, q
// 
// where x, y, z are (physical) position of the particle, ux, uy, uz
// are particle normalized momenta, and q is the particle charge.
// 
// Written by: Brian J. Albright, X-1, LANL January, 2006

int
link_boundary( link_boundary_t     * lb,
               particle_t          * r, 
               particle_mover_t    * pm,
               field_t             * f,
               accumulator_t       * a, 
               const grid_t        * g,
               species_t           * sp, 
               particle_injector_t * pi,
               rng_t               * rng,
               int                   face ) {
  static FILE *fp = NULL; 
  int ix, iy, iz;
  double x, y, z;
  char fname[512];

  if( !fp ) { 
    // FIXME: TIMESTAMP THE FILE ... TIMESTAMP SHOULD BE SHARED BY ALL
    // NODES
    //
    // char timestamp[16];
    // make_timestamp( timestamp );
    // sprintf( fname, "link.%s.%d", timestamp_string, world_rank );

    sprintf( fname, "%s.%d", lb->fbase, world_rank ); 
    fp = fopen( fname, "r" );
    if( fp )
      ERROR(( "File %s already exists (probably from a earlier run or "
              "recovery) ... Please move it or change link filename", fname ));
    fp = fopen( fname, "w" );
    if( !fp ) ERROR(( "Could not open file %s", fname ));
    
    fprintf( fp, "%% %17.0f\n", lb->n_out );
  }

  // FIXME: THIS CALCULATION SHOULD BE MADE MORE RIGOROUS!
  // FIXME: USE GRID_H HELPER HERE ... THEY ARE FASTER.

  iz = r->i/((g->nx+2)*(g->ny+2));
  iy = (r->i-iz*(g->nx+2)*(g->ny+2))/(g->nx+2);
  ix = r->i-(g->nx+2)*(iy+(g->ny+2)*iz);

  x = g->x0 + ((ix-1)+(r->dx+1)*0.5)*g->dx;
  y = g->y0 + ((iy-1)+(r->dy+1)*0.5)*g->dy;
  z = g->z0 + ((iz-1)+(r->dz+1)*0.5)*g->dz;

  // ASCII write of particle data as per Tom's request
  // FIXME: perform I/O without calls to fprintf and need for stdio.h
  // Fields: species id, x, y, z, ux, uy, uz, q

  fprintf( fp, "%d %e %e %e %e %e %e %e\n",
           sp->id, x, y, z, r->ux, r->uy, r->uz, sp->q*r->w );
  lb->n_out++;

  accumulate_rhob( f, r, g, sp->q );

  return 0;
}
