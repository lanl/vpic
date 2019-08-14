// Written by:
//   Kevin J. Bowers, Ph.D.
//   Plasma Physics Group (X-1)
//   Applied Physics Division
//   Los Alamos National Lab
// March/April 2004 - Original version

#include "vpic.h"

// FIXME: MOVE THIS INTO VPIC.HXX TO BE TRULY INLINE

void
vpic_simulation::inject_particle( species_t * sp,
                                  double x,  double y,  double z,
                                  double ux, double uy, double uz,
                                  double w,  double age,
                                  int update_rhob )
{
  int ix, iy, iz;

  // Check input parameters
  if( !accumulator_array ) ERROR(( "Accumulator not setup yet" ));
  if( !sp                ) ERROR(( "Invalid species" ));
  if( w < 0              ) ERROR(( "inject_particle: w < 0" ));

  const double x0 = (double)grid->x0, y0 = (double)grid->y0, z0 = (double)grid->z0;
  const double x1 = (double)grid->x1, y1 = (double)grid->y1, z1 = (double)grid->z1;
  const int    nx = grid->nx,         ny = grid->ny,         nz = grid->nz;

  // Do not inject if the particle is strictly outside the local domain
  // or if a far wall of local domain shared with a neighbor
  // FIXME: DO THIS THE PHASE-3 WAY WITH GRID->NEIGHBOR
  // NOT THE PHASE-2 WAY WITH GRID->BC

  if( (x<x0) | (x>x1) | ( (x==x1) & (grid->bc[BOUNDARY(1,0,0)]>=0 ) ) ) return;
  if( (y<y0) | (y>y1) | ( (y==y1) & (grid->bc[BOUNDARY(0,1,0)]>=0 ) ) ) return;
  if( (z<z0) | (z>z1) | ( (z==z1) & (grid->bc[BOUNDARY(0,0,1)]>=0 ) ) ) return;

  // This node should inject the particle

  if( sp->np>=sp->max_np ) ERROR(( "No room to inject particle" ));

  // Compute the injection cell and coordinate in cell coordinate system
  // BJA:  Note the use of double precision here for accurate particle
  //       placement on large meshes.

  // The ifs allow for injection on the far walls of the local computational
  // domain when necessary

  x  = ((double)nx)*((x-x0)/(x1-x0)); // x is rigorously on [0,nx]
  ix = (int)x;                        // ix is rigorously on [0,nx]
  x -= (double)ix;                    // x is rigorously on [0,1)
  x  = (x+x)-1;                       // x is rigorously on [-1,1)
  if( ix==nx ) x = 1;                 // On far wall ... conditional move
  if( ix==nx ) ix = nx-1;             // On far wall ... conditional move
  ix++;                               // Adjust for mesh indexing

  y  = ((double)ny)*((y-y0)/(y1-y0)); // y is rigorously on [0,ny]
  iy = (int)y;                        // iy is rigorously on [0,ny]
  y -= (double)iy;                    // y is rigorously on [0,1)
  y  = (y+y)-1;                       // y is rigorously on [-1,1)
  if( iy==ny ) y = 1;                 // On far wall ... conditional move
  if( iy==ny ) iy = ny-1;             // On far wall ... conditional move
  iy++;                               // Adjust for mesh indexing

  z  = ((double)nz)*((z-z0)/(z1-z0)); // z is rigorously on [0,nz]
  iz = (int)z;                        // iz is rigorously on [0,nz]
  z -= (double)iz;                    // z is rigorously on [0,1)
  z  = (z+z)-1;                       // z is rigorously on [-1,1)
  if( iz==nz ) z = 1;                 // On far wall ... conditional move
  if( iz==nz ) iz = nz-1;             // On far wall ... conditional move
  iz++;                               // Adjust for mesh indexing

  particle_t * p = sp->p + (sp->np++);
  p->dx = (float)x; // Note: Might be rounded to be on [-1,1]
  p->dy = (float)y; // Note: Might be rounded to be on [-1,1]
  p->dz = (float)z; // Note: Might be rounded to be on [-1,1]
  p->i  = VOXEL(ix,iy,iz, nx,ny,nz);
  p->ux = (float)ux;
  p->uy = (float)uy;
  p->uz = (float)uz;
  p->w  = w;

  if( update_rhob ) accumulate_rhob( field_array->f, p, grid, -sp->q );

  if( age!=0 ) {
    if( sp->nm>=sp->max_nm )
      WARNING(( "No movers available to age injected  particle" ));
    particle_mover_t * pm = sp->pm + sp->nm;
    age *= grid->cvac*grid->dt/sqrt( ux*ux + uy*uy + uz*uz + 1 );
    pm->dispx = ux*age*grid->rdx;
    pm->dispy = uy*age*grid->rdy;
    pm->dispz = uz*age*grid->rdz;
    pm->i     = sp->np-1;
    sp->nm += move_p( sp->p, pm, accumulator_array->a, grid, sp->q );
  }
}

// Add capability to modify certain fields "on the fly" so that one
// can, e.g., extend a run, change a quota, or modify a dump interval
// without having to rerun from the start.
//
// File is specified in arg 3 slot in command line inputs.  File is in
// ASCII format with each field in the form: field val [newline].
//
// Allowable values of field variables are: num_steps, quota,
// checkpt_interval, hydro_interval, field_interval, particle_interval
// ndfld, ndhyd, ndpar, ndhis, ndgrd, head_option,
// istride, jstride, kstride, stride_option, pstride
//
// [x]_interval sets interval value for dump type [x].  Set interval
// to zero to turn off dump type.
//
// FIXME-KJB: STRIP_CMDLINE ALLOWS SOMETHING CLEANER AND MORE POWERFUL

#define SETIVAR( V, A, S ) do {                                          \
    V=(A);                                                               \
    if ( rank()==0 ) log_printf( "*** Modifying %s to value %d", S, A ); \
  } while(0)

#define SETDVAR( V, A, S ) do {                                           \
    V=(A);                                                                \
    if ( rank()==0 ) log_printf( "*** Modifying %s to value %le", S, A ); \
  } while(0)

#define ITEST( V, N, A ) \
  if( sscanf( line, N " %d",  &iarg )==1 ) SETIVAR( V, A, N )

#define DTEST( V, N, A ) \
  if( sscanf( line, N " %le", &darg )==1 ) SETDVAR( V, A, N )

void
vpic_simulation::modify( const char *fname ) {
  FILE *handle=NULL;
  char line[128];
  int iarg=0;
  double darg=0;

  // Open the modfile
  handle = fopen( fname, "r" );
  if( !handle ) ERROR(( "Modfile read failed" ));
  // Parse modfile
  while( fgets( line, 127, handle ) ) {
    DTEST( quota,             "quota",             darg );
    ITEST( num_step,          "num_step",          iarg );
    ITEST( checkpt_interval,  "checkpt_interval",  (iarg<0 ? 0 : iarg) );
    ITEST( hydro_interval,    "hydro_interval",    (iarg<0 ? 0 : iarg) );
    ITEST( field_interval,    "field_interval",    (iarg<0 ? 0 : iarg) );
    ITEST( particle_interval, "particle_interval", (iarg<0 ? 0 : iarg) );

    // RFB: No existing decks I know of use these quantities. They are for
    // striding output and are a legacy hangover from 407
#ifdef EXTEND_MODIFY
    ITEST( ndfld, "ndfld", (iarg<0 ? 0 : iarg) );
    ITEST( ndhyd, "ndhyd", (iarg<0 ? 0 : iarg) );
    ITEST( ndpar, "ndpar", (iarg<0 ? 0 : iarg) );
    ITEST( ndhis, "ndhis", (iarg<0 ? 0 : iarg) );
    ITEST( ndgrd, "ndgrd", (iarg<0 ? 0 : iarg) );
    ITEST( head_option, "head_option", (iarg<0 ? 0 : iarg) );
    ITEST( istride, "istride", (iarg<1 ? 1 : iarg) );
    ITEST( jstride, "jstride", (iarg<1 ? 1 : iarg) );
    ITEST( kstride, "kstride", (iarg<1 ? 1 : iarg) );
    ITEST( stride_option, "stride_option", (iarg<1 ? 1 : iarg) );
    ITEST( pstride, "pstride", (iarg<1 ? 1 : iarg) );
    ITEST( stepdigit, "stepdigit", (iarg<0 ? 0 : iarg) );
    ITEST( rankdigit, "rankdigit", (iarg<0 ? 0 : iarg) );
#endif

  }
}

#undef SETIVAR
#undef SETDVAR
#undef ITEST
#undef DTEST

#if defined(ENABLE_OPENSSL)
#include "../util/checksum.h"

void vpic_simulation::checksum_fields(CheckSum & cs) {
  checkSumBuffer<field_array_t>(field_array, grid->nv, cs, "sha1");

  if(nproc() > 1) {
    const unsigned int csels = cs.length*nproc();
    unsigned char * sums(NULL);

    if(rank() == 0) {
      sums = new unsigned char[csels];
    } // if

    // gather sums from all ranks
    mp_gather_uc(cs.value, sums, cs.length);

    if(rank() == 0) {
      checkSumBuffer<unsigned char>(sums, csels, cs, "sha1");
      delete[] sums;
    } // if
  } // if
} // vpic_simulation::output_checksum_fields

void vpic_simulation::output_checksum_fields() {
  CheckSum cs;
  checkSumBuffer<field_array_t>(field_array, grid->nv, cs, "sha1");

  if(nproc() > 1) {
    const unsigned int csels = cs.length*nproc();
    unsigned char * sums(NULL);

    if( rank() == 0) {
      sums = new unsigned char[csels];
    } // if

    // gather sums from all ranks
    mp_gather_uc(cs.value, sums, cs.length);

    if( rank() == 0) {
      checkSumBuffer<unsigned char>(sums, csels, cs, "sha1");
      MESSAGE(("FIELDS SHA1CHECKSUM: %s", cs.strvalue));
      delete[] sums;
    } // if
  }
  else {
    MESSAGE(("FIELDS SHA1CHECKSUM: %s", cs.strvalue));
  } // if
} // vpic_simulation::output_checksum_fields

void vpic_simulation::checksum_species(const char * species, CheckSum & cs) {
  species_t * sp = find_species_name(species, species_list);
  if(sp == NULL) {
    ERROR(("Invalid species name \"%s\".", species));
  } // if

  checkSumBuffer<particle_t>(sp->p, sp->np, cs, "sha1");

  if(nproc() > 1) {
    const unsigned int csels = cs.length*nproc();
    unsigned char * sums(NULL);

	if(rank() == 0) {
    	sums = new unsigned char[csels];
	} // if

	// gather sums from all ranks
	mp_gather_uc(cs.value, sums, cs.length);

	if(rank() == 0) {
		checkSumBuffer<unsigned char>(sums, csels, cs, "sha1");
		MESSAGE(("SPECIES \"%s\" SHA1CHECKSUM: %s", species, cs.strvalue));
		delete[] sums;
    } // if
  } // if
} // vpic_simulation::checksum_species

void vpic_simulation::output_checksum_species(const char * species) {
  species_t * sp = find_species_name(species, species_list);
  if(sp == NULL) {
    ERROR(("Invalid species name \"%s\".", species));
  } // if

  CheckSum cs;
  checkSumBuffer<particle_t>(sp->p, sp->np, cs, "sha1");

  if(nproc() > 1) {
    const unsigned int csels = cs.length*nproc();
    unsigned char * sums(NULL);

	if( rank() == 0) {
    	sums = new unsigned char[csels];
	} // if

	// gather sums from all ranks
	mp_gather_uc(cs.value, sums, cs.length);

	if( rank() == 0) {
		checkSumBuffer<unsigned char>(sums, csels, cs, "sha1");
		MESSAGE(("SPECIES \"%s\" SHA1CHECKSUM: %s", species, cs.strvalue));
		delete[] sums;
    } // if
  }
  else {
    MESSAGE(("SPECIES \"%s\" SHA1CHECKSUM: %s", species, cs.strvalue));
  } // if
} // vpic_simulation::output_checksum_species

#endif // ENABLE_OPENSSL
