// Beam-plasma interaction
//
// Based on geometry and physics parameters provided by A. Sgro (LANL)
//
// This input deck was written by:
//   Kevin J Bowers, Ph.D.
//   Plasma Physics Group (X-1)
//   Applied Physics Division
//   Los Alamos National Lab
// April 2004 - revised from earlier V4PIC versions
//
// Updated Feb 2009 to work with new deck formatting - BJA

begin_globals {
  double energies_interval;
  double fields_interval;
  double bhydro_interval;
  double ehydro_interval;
  double ihydro_interval;
  double bparticle_interval;
  double eparticle_interval;
  double iparticle_interval;
  double restart_interval;

  double r_beam;
  double I_beam;
  double gamma_beam;
  double uth_beam;
  double ninj_beam;
};

begin_initialization {
  // Units (rationalized-MKSA)
  double cvac = 299792458;       // Speed of light (in m/s)
  double eps0 = 8.854187817e-12; // Permittivity (in F/m)
  double ec   = 1.602176462e-19; // Fundamental charge (in C)
  double me   = 9.10938188e-31;  // Electron mass (in kg)
  double mi   = 1.67262158e-27;  // Ion mas (in kg) ... Protons

  // Grid specification
  double lx        = 1e-2;       // x-size (in m)
  double ly        = 1e-2;       // y-size (in m)
  double lz        = 10e-2;      // z-size (in m)
  double nx        = 10;         // x-resolution
  double ny        = 10;         // y-resolution
  double nz        = 500;        // z-resolution
  double cfl_req   = 0.97;       // How close to Courant should we run
  double wpedt_max = 0.2;        // Max timestep if not Courant limited
  double tau       = 10*lz/cvac; // How long to run (in s)
  double damp      = 0.01;       // Radiation damping

  // Beam
  double r_beam    = 1e-3;     // Beam radius (in m)
  double I_beam    = 4000;     // Beam current (in A)
  double E_beam    = 5.5e6*ec; // Beam energy, excludes rest energy (in J)
  double kT_beam   = 0.026*ec; // Temperature in the beam frame (in J)
  double nipf_beam = 60;       // Particles to inject per cell face for beam

  // Plasma
  double z_plas    = 1e-2;     // Where the plasma begins (in m)
  double b_plas    = 1e-3;     // Isolation buffer (in m)
  double n_plas    = 1e19;     // Pair density (in m^-3)
  double kTi       = 0.026*ec; // Ion temperature (in J)
  double kTe       = 0.026*ec; // Electron temperature (in J)
  double nppc_plas = 100;      // Macro pairs per cell

  // Derived parameters
  double dx   = lx/nx;                           // x-cell size (in m)
  double dy   = ly/ny;                           // y-cell size (in m)
  double dz   = lz/nz;                           // z-cell size (in m)
  double dg   = courant_length(lx,ly,lz,nx,ny,nz); // Courant length (in m)
  double wpe  = sqrt(ec*ec*n_plas/(me*eps0));    // Plasma frequency (in rad/s)
  double dt   = cfl_req*dg/cvac;                 // Timestep (in s)
  if( wpe*dt>wpedt_max ) dt = wpedt_max/wpe;     // Override if wpe limited
  double vol  = ((nx==1) ? (lx-2*b_plas) : (lx-b_plas)) *
                ((ny==1) ? (ly-2*b_plas) : (ly-b_plas)) *
                (lz-b_plas-z_plas);              // Volume of plasma (in m^3)
  double np_phys = vol*n_plas;                   // Number of physical pairs
  double np   = round(nppc_plas*vol/(dx*dy*dz)); // Number of macro pairs
  double qp   = ec*np_phys/np;                   // Charge per macro particle
  double uti  = sqrt(kTi/mi)/cvac;               // H+ thermal norm momentum
  double ute  = sqrt(kTe/me)/cvac;               // e- thermal norm momentum
  double ninj = round(nipf_beam*0.25*M_PI*r_beam*r_beam/(dx*dy));
  if( nx==1 ) ninj *= 2;
  if( ny==1 ) ninj *= 2;
  
  // Setup high level parameters
  num_step             = int(tau/dt);
  status_interval      = 10;
  clean_div_e_interval = 50;
  clean_div_b_interval = 50;
  sync_shared_interval = 50;

  // Setup diagnostic dump intervals
  global->energies_interval  = 10;
  global->fields_interval    = 50;
  global->bhydro_interval    = 50;
  global->ehydro_interval    = 50;
  global->ihydro_interval    = 50;
  global->bparticle_interval = 250;
  global->eparticle_interval = 250;
  global->iparticle_interval = 250;
  global->restart_interval   = 0; // Do not dump

  // Setup custom beam injection parameters
  global->r_beam     = r_beam;
  global->I_beam     = I_beam;
  global->gamma_beam = 1 + E_beam/(me*cvac*cvac);
  global->uth_beam   = sqrt(kT_beam/me);
  global->ninj_beam  = ninj;

  // Setup the grid and domain decomposition
  define_units( cvac, eps0 );
  define_timestep( dt );
  define_periodic_grid( (nx==1) ? -0.5*lx : 0,  (ny==1) ? -0.5*ly : 0,  0,
                        (nx==1) ?  0.5*lx : lx, (ny==1) ?  0.5*ly : ly, lz,
                        nx, ny, nz,
                        1, 1, nproc() ); // Partition along z

  if( nx>1 ) { // Put in symmetry planes, metal if x-dimension is resolved
    set_domain_field_bc(    BOUNDARY(-1,0,0), symmetric_fields  );
    set_domain_particle_bc( BOUNDARY(-1,0,0), reflect_particles );
    set_domain_field_bc(    BOUNDARY( 1,0,0), anti_symmetric_fields  );
    set_domain_particle_bc( BOUNDARY( 1,0,0), absorb_particles );
  }

  if( ny>1 ) { // Put in symmetry planes, metal if y-dimension is resolved
    set_domain_field_bc(    BOUNDARY(0,-1,0), symmetric_fields  );
    set_domain_particle_bc( BOUNDARY(0,-1,0), reflect_particles );
    set_domain_field_bc(    BOUNDARY(0, 1,0), anti_symmetric_fields  );
    set_domain_particle_bc( BOUNDARY(0, 1,0), absorb_particles );
  }
  
  if( rank()==0 ) { // Put in beam launcher
    set_domain_field_bc(    BOUNDARY(0,0,-1), absorb_fields );
    set_domain_particle_bc( BOUNDARY(0,0,-1), absorb_particles );
  }

  if( rank()==nproc()-1 ) { // Put in Beam absorber
    set_domain_field_bc(    BOUNDARY(0,0,1), absorb_fields );
    set_domain_particle_bc( BOUNDARY(0,0,1), absorb_particles );
  }

  // Create the materials
  define_material( "vacuum", 1 );

  // Create the fields
  define_field_array( NULL, damp );

  // Create the species (20%-50% overhead in case of non-uniformity)
  /**/                   define_species( "beam", -ec, me, 1.5*ninj*lz/(cvac*dt*nproc()), -1, 20, 1 );
  species_t * electron = define_species( "e-",   -ec, me, 1.2*np/nproc(),                -1, 20, 1 );
  species_t * ion      = define_species( "H+",    ec, mi, 1.2*np/nproc(),                -1, 40, 1 );

  // Create the initial plasma
  // Note: This assumes ranks 0 and nproc()-1 have particles in them
  double xmin = (nx==1 ? -0.5*lx+b_plas : 0);
  double ymin = (ny==1 ? -0.5*ly+b_plas : 0);
  double zmin = (rank()==0 ? z_plas : grid->z0);
  double xmax = (nx==1 ? 0.5*lx-b_plas : lx-b_plas );
  double ymax = (ny==1 ? 0.5*ly-b_plas : ly-b_plas );
  double zmax = (rank()==nproc()-1 ? lz-b_plas : grid->z0+grid->dz*grid->nz);

  repeat( np*(zmax-zmin)/(lz-b_plas-z_plas) ) {
    double x = uniform( rng(0), xmin, xmax );
    double y = uniform( rng(0), ymin, ymax );
    double z = uniform( rng(0), zmin, zmax );
    inject_particle( electron, x, y, z,
                     normal( rng(0), 0, ute ),
                     normal( rng(0), 0, ute ),
                     normal( rng(0), 0, ute ), qp, 0, 0 );
    inject_particle( ion, x, y, z,
                     normal( rng(0), 0, uti ),
                     normal( rng(0), 0, uti ),
                     normal( rng(0), 0, uti ), qp, 0, 0 );
  }
  
  // Output diagnostics
  sim_log( "" );
  sim_log( "Numerical parameters" );
  sim_log( "Box size[cm] = " << lx*100 << " x " << ly*100 << " x " << lz*100 );
  sim_log( "Resolution = " << nx << " x " << ny << " x " << nz );
  sim_log( "Cell size[mm] = " << dx*1e3 << " x " << dy*1e3 << " x " << dz*1e3 );
  sim_log( "Duration[ns] = " << num_step*dt*1e9 );
  sim_log( "Number of step = " << num_step );
  sim_log( "Time step[ps] = " << dt*1e12 );
  sim_log( "damping = " << damp );
  sim_log( "courant = " << cvac*dt/dg );
  sim_log( "" );
  sim_log( "Beam parameters" );
  sim_log( "E_beam[MeV] = " << E_beam/(ec*1e6) );
  sim_log( "I_beam[kA]  = " << I_beam/1e3 );
  sim_log( "r_beam[mm]  = " << r_beam*1e3 );
  sim_log( "T_beam[eV]  = " << kT_beam/ec );
  sim_log( "ninj_beam   = " << ninj );
  sim_log( "" );
  sim_log( "Plasma Parameters" );
  sim_log( "ne[cm^-3] = "  << n_plas/1e6 << ", kTe[eV] = " << kTe/ec );
  sim_log( "ni[cm^-3] = "  << n_plas/1e6 << ", kTi[eV] = " << kTi/ec );
  sim_log( "mi/me = "      << mi/me );
  sim_log( "z_plas[cm] = " << z_plas*100 );
  sim_log( "b_plas[cm] = " << b_plas*100 );
  sim_log( "nppc_plas = "  << nppc_plas );
  sim_log( "" );
  sim_log( "Miscellaneous" );
  sim_log( "num_proc = " << nproc() );
  sim_log( "np = " << np );
  sim_log( "" );
}

begin_diagnostics {

  if( step()==0 ) {
    dump_grid("grid");
    dump_materials("materials");
    dump_species("species");
    dump_energies("energies",0);
    dump_fields("fields");
    dump_hydro("beam","bhydro");
    dump_hydro("e-","ehydro");
    dump_hydro("H+","ihydro");
  }

# define should_dump(x) ( step()>0 &&                               \
                          global->x##_interval>0 &&                 \
                          remainder(step(),global->x##_interval)==0 )
  if( should_dump(energies)  ) dump_energies("energies");
  if( should_dump(fields)    ) dump_fields("fields");
  if( should_dump(bhydro)    ) dump_hydro("beam","bhydro");
  if( should_dump(ehydro)    ) dump_hydro("e-",  "ehydro");
  if( should_dump(ihydro)    ) dump_hydro("H+",  "ihydro");
  if( should_dump(bparticle) ) dump_particles("beam","bparticle");
  if( should_dump(eparticle) ) dump_particles("e-",  "eparticle");
  if( should_dump(iparticle) ) dump_particles("H+",  "iparticle");
  if( should_dump(restart)   ) checkpt("restart",step()); 
# undef should_dump

}

begin_particle_injection {
  return;
  if( rank()!=0 ) return;

  species_t * beam = find_species_name( "beam", species_list );
  double ub = sqrt( global->gamma_beam*global->gamma_beam - 1 );
  double qb = 0.25*global->I_beam*grid->dt/global->ninj_beam;
  if( grid->nx==1 ) qb *= 2;
  if( grid->ny==1 ) qb *= 2;

  repeat(global->ninj_beam) {
    double x, y, z, ux, uy, uz, age;

    do {
      x = uniform( rng(0), -1, 1 );
      y = uniform( rng(0), -1, 1 );
    } while( x*x + y*y > 1 );
    if( grid->nx>1 && x<0 ) x=-x;
    if( grid->ny>1 && y<0 ) y=-y;
    x *= global->r_beam;
    y *= global->r_beam;
    z = 0;

    // FIXME: THIS DISTRIBUTION IS SLIGHTLY WRONG (BOTH RELATIVISTIC AND
    // NUMERICALLY FOR THE UZ COMPONENT (SEE MAXWELLIAN_RAND).
    ux = normal( rng(0), 0, global->uth_beam );
    uy = normal( rng(0), 0, global->uth_beam );
    uz = normal( rng(0), 0, global->uth_beam );
    uz = global->gamma_beam*uz + sqrt(1+ux*ux+uy*uy+uz*uz)*ub;

    age = uniform( rng(0), 0, 1 );

    inject_particle( beam, x, y, z, ux, uy, uz, qb, age, 1 );
  }
}

begin_current_injection {
}

begin_field_injection {
}


begin_particle_collisions { 
  // No collisions for this simulation 
}

