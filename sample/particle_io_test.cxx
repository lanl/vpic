begin_globals {
};

begin_initialization {


  num_step = 1;                      // Run for one step so we can do output once
  status_interval = 100;             // Basically don't print status
  sync_shared_interval = status_interval;
  clean_div_e_interval = status_interval;
  clean_div_b_interval = status_interval;

  define_units(1, 1);                // speed of light and eps0
  define_timestep( 0.95/sqrt(3.0) ); // This implies dx=dy=dz = 1

  int gnx = 8;
  int gny = 4;
  int gnz = 4;

  int topox = nproc();
  int topoy = 1;
  int topoz = 1;

  define_periodic_grid(0,0,0,              // Low corner
                       gnx,gny,gnz,        // High corner
                       gnx,gny,gnz,        // Resolution
                       topox,topoy,topoz); // Topology

  // Space is by default filled with first material defined
  define_material("vacuum",1.0);

  // Create the field array
  define_field_array(NULL, 0.00);

  // Create species
  species_t* elec = define_species("electron",  // name
                                          -1.,  // charge
                                           1.,  // mass
                                            2,  // maximum number of local particles
                                           -1,  // automatic numbner of particle movers
                                            0,  // sort never
                                            1); // sort mode out-of-place

  species_t* fill = define_species("fill",  // name
                                      -1.,  // charge
                                       1.,  // mass
                  100*gnx*gny*gnz/nproc(),  // maximum number of local particles
                                       -1,  // automatic numbner of particle movers
                                        0,  // sort never
                                        1); // sort mode out-of-place

  species_t* prot = define_species("proton",  // name
                                         1.,  // charge
                                         1.,  // mass
                                          1,  // maximum number of local particles
                                         -1,  // automatic numbner of particle movers
                                          0,  // sort never
                                          1); // sort mode out-of-place


  // Set field values
  set_region_field(everywhere, 0.,0.,0., 0.,0.,0.);

  // Generate particles
  seed_entropy(rank()); // different random numbers on different ranks

  if(rank() == 0) {
    const float x = 3.9;
    const float y = 2.1;
    const float z = 2.2;

    const float ux = 0.999; // nearly c
    const float uy = 0.;
    const float uz = 0.;

    inject_particle(elec, x,y,z,    ux,uy,uz, 1., 0, 0);
    inject_particle(elec, x,y-1.,z, ux,uy,uz, 1., 0, 0);
    inject_particle(prot, x,y,z,    ux,uy,uz, 1., 0, 0);
  }

  for(int n=0; n<fill->max_np; n++) {
    // Pick a uniform random location in the local domain
    const float x = uniform( rng(0), grid->x0, grid->x1 );
    const float y = uniform( rng(0), grid->y0, grid->y1 );
    const float z = uniform( rng(0), grid->z0, grid->z1 );

    // Pick random velocity from normal distribution to get thermal VDF
    const float vdriftx = 0.1;
    const float vdrifty = 0.;
    const float vdriftz = 0.;
    const float vthx = 0.1;
    const float vthy = 0.2;
    const float vthz = 0.3;
    const float ux = normal( rng(0), vdriftx, vthx );
    const float uy = normal( rng(0), vdrifty, vthy );
    const float uz = normal( rng(0), vdriftz, vthz );
    inject_particle(fill, x,y,z, ux,uy,uz, 1., 0, 0);
  }

#ifdef VPIC_GLOBAL_PARTICLE_ID
  // Add global IDs by moving to a tracer species;
  sim_log("Making all eletrons tracers");
  species_t * sp2 = make_tracers_by_percentage(elec, 100.0, Tracertype::Move);
#endif

  sim_log("Done with setup");
  fflush(NULL);
}

begin_diagnostics {
  dump_particles_hdf5("electron", "elec");
  dump_particles_hdf5("fill", "fill");
  dump_particles_hdf5("proton", "prot");
}

begin_particle_injection {
}

begin_current_injection {
}

begin_field_injection {
}

begin_particle_collisions {
}


