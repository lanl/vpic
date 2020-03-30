// Demonstrates how to use the Takizuka-Abe collision model to simulate
// an electron beam slowing down in a background plasma.

double
mu_prime(double x) {
  return 2*sqrt(x/M_PI)*exp(-x);
}  

double
mu(double x) {
  // Computes the definite integral
  // mu(x) = 2/sqrt(pi) int_0^x sqrt(y) exp(-y) dy

  double dx, x0 = 0, z = 0;
  if(x > 10) {
    x0 = 10;
    z = mu(x0);
  }

  dx = (x-x0)/1000;
  for(int i=0 ; i < 1000 ; ++i) {
    z  += 0.5*dx*mu_prime(x0);
    x0 += dx;
    z  += 0.5*dx*mu_prime(x0);
  }
  return z;
}

begin_globals {
  int ncoll;
};

begin_initialization {

  // Simulation parameters

  double mime        = 100;    // Ion to electron mass ratio
  double vthe        = 0.1;    // Electron thermal speed in c
  double Ebeam       = 0.5;    // Beam energy in Te
  double nbeam       = 0.01;   // Beam density relative to background
  double nu0         = 1e-2;   // Base collision frequency
  double beam_weight = 10;     // Stat. weight of the beam realtive to background
  int nppc           = 100000; // Particles/species/cell in the background
  int ncoll          = 1;      // How often to collide particles in steps
  int nstep          = 10;     // Number of steps to run for

  // Setup
  double dt = 0.98 * courant_length(1, 1, 1, 3, 3, 3);

  define_units( 1, 1 );
  define_timestep( dt );
  define_periodic_grid( 0, 0, 0,   // Grid low corner
                        1, 1, 1,   // Grid high corner
                        3, 3, 3,   // Grid resolution
                        1, 1, 1 ); // Processor configuration
  define_material("vacuum",1.0,1.0,0.0);
  define_field_array();
  set_region_field(everywhere, 0, 0, 0, 0, 0, 0);

  // Define the particle species.        name         q    m   Nmax  Nm Sort In-place
  int Nback = grid->nx*grid->ny*grid->nz * nppc;
  int Nbeam = nbeam*Nback/beam_weight;
  species_t * electron = define_species( "electron", -1,    1, Nback, -1, ncoll, 0 ),
            * ion      = define_species( "ion",       1, mime, Nback, -1, ncoll, 0 ),
            * beam     = define_species( "beam",     -1,    1, Nbeam, -1, ncoll, 0 );


  // Load the particles.
  double weight = 1/(double) Nback;
  double vthi   = vthe / sqrt(mime);
  double vbeam  = sqrt(2*Ebeam)*vthe;

  repeat(Nback){

    double x = uniform( rng(0), 0, 1 );
    double y = uniform( rng(0), 0, 1 );
    double z = uniform( rng(0), 0, 1 );

    inject_particle( electron, x, y, z,
                     normal(  rng(0), 0, vthe ),
                     normal(  rng(0), 0, vthe ),
                     normal(  rng(0), 0, vthe ),
                     weight, 0, 0 );

    inject_particle( ion, x, y, z,
                     normal(  rng(0), 0, vthi ),
                     normal(  rng(0), 0, vthi ),
                     normal(  rng(0), 0, vthi ),
                     weight, 0, 0 );

  }

  repeat(Nbeam){

    double x = uniform( rng(0), 0, 1 );
    double y = uniform( rng(0), 0, 1 );
    double z = uniform( rng(0), 0, 1 );
    inject_particle( beam, x, y, z, vbeam, 0, 0, weight*beam_weight, 0, 0 );

  }

  // Create the collision operators
  // 
  // On input, nu0 refers to the base collision frequency of the thermal
  // plasma, but for the collision operator we need to normalize this
  // by mc2 instead of T. This amounts to rescaling nu0 by (vthe/c)^3

  double nu0_c3 = nu0 * pow(vthe, 3);

  define_collision_op(takizuka_abe("ee_bulk", electron, electron, entropy, nu0_c3, ncoll));
  define_collision_op(takizuka_abe("ei_bulk", electron, ion     , entropy, nu0_c3, ncoll));
  define_collision_op(takizuka_abe("ee_beam", electron, beam    , entropy, nu0_c3, ncoll));
  define_collision_op(takizuka_abe("ei_beam", ion,      beam    , entropy, nu0_c3, ncoll));
  // Ignore beam-beam collisions for benchmarking.

  // Finalize the simulation paramters
  num_step             = nstep;
  status_interval      = 200;
  sync_shared_interval = -status_interval/2;
  clean_div_e_interval = -status_interval/2;
  clean_div_b_interval = -status_interval/2;

  // For informational purposes, compute the transport rates for the beam.
  // See e.g., the Plasma Formulary for these rates.
  double x;
  double nu_beam = nu0 * pow(Ebeam, -1.5);

  // Rates due to collisions with electrons.
  x = Ebeam;
  double nu_se = -2 * mu(x) * nu_beam;
  double nu_de = 2 * ( (1-0.5/x)*mu(x) + mu_prime(x) ) * nu_beam;
  double nu_te = -2 * ( mu(x) - mu_prime(x) ) * nu_beam;

  // Rates due to collisions with ions
  x = Ebeam * mime;
  double nu_si = -(1+1/mime)*mu(x) * nu_beam;
  double nu_di = 2 * ( (1-0.5/x)*mu(x) + mu_prime(x) ) * nu_beam;
  double nu_ti = -2 * ( mu(x)/mime - mu_prime(x) ) * nu_beam;

  // Log the paramters.

  sim_log("Electron Beam Parameters");
  sim_log("-----------------------------------------------");
  sim_log("Ion / electron mass ratio    = " << mime         );
  sim_log("Beam energy / temperature    = " << Ebeam        );
  sim_log("Beam density / background    = " << nbeam        );
  sim_log("Beam weights / backgorund    = " << beam_weight  );
  sim_log("Background particles         = " << Nback        );
  sim_log("Beam particles               = " << Nbeam        );
  sim_log("-----------------------------------------------");
  sim_log("Theoretical slow down  rate  = " << nu_se + nu_si);
  sim_log("Theoretical deflection rate  = " << nu_de + nu_di);
  sim_log("Theoretical energy loss rate = " << nu_te + nu_ti);
  sim_log("-----------------------------------------------");

  // Set globals
  global->ncoll = ncoll;

}

begin_diagnostics {

  if( step() % global->ncoll == 0) {

    // Find the beam and collect the hydro moments
    species_t * beam = find_species("beam");
    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, beam, interpolator_array );
    synchronize_hydro_array( hydro_array );

    // Average over the volume.
    int i, j, k;
    hydro_t avg = {0};
    double ncells = grid->nx * grid->ny * grid->nz; 
    for( k=1; k <= grid->nz ; ++k)
      for( j=1; j <= grid->ny ; ++j)
        for( i=1; i <= grid->nx ; ++i) {
	    	  
	  hydro_t cell = hydro_array->h[ voxel(i,j,k) ];
	  avg.jx  += cell.jx  / ncells;
	  avg.jy  += cell.jy  / ncells;
	  avg.jz  += cell.jz  / ncells;
	  avg.ke  += cell.ke  / ncells;
	  avg.txx += cell.txx / ncells;
	  avg.tyy += cell.tyy / ncells;
	  avg.tzz += cell.tzz / ncells;

	}

    // Now log the slow-down rate.
    static double jx = avg.jx;
    double nu_s = (2*(avg.jx-jx)/(avg.jx+jx))/(grid->dt*global->ncoll);
    jx = avg.jx;

    // Log the deflection rate and the energy loss rate.
    static double Tperp = avg.tyy + avg.tzz;
    static double Tbeam = avg.txx + avg.tyy + avg.tzz;
    
    double Tp = avg.tyy + avg.tzz;
    double T  = Tp + avg.txx;
    double nu_d = (2*(Tp-Tperp)/(T+Tbeam))/(grid->dt*global->ncoll);
    double nu_t = (2*(T -Tbeam)/(T+Tbeam))/(grid->dt*global->ncoll);
    Tperp = Tp;
    Tbeam = T;

    sim_log("Beam current = " << jx << ", slow down rate   = " << nu_s);
    sim_log("Beam Tperp   = " << Tp << ", deflection rate  = " << nu_d);
    sim_log("Beam energy  = " << T  << ", energy loss rate = " << nu_t);

  } 
		  
}	     
	     	  
begin_particle_injection {
}

begin_current_injection {
}

begin_field_injection {
  // Force field to 0 to prevent grid instabilities since dx may be >> Debye
  set_region_field(everywhere, 0, 0, 0, 0, 0, 0);
}

begin_particle_collisions {

}
