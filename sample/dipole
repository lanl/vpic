// Simulate Hertzian dipole radiation
//
// This input deck was written by:
//   Kevin J Bowers, Ph.D.
//   Plasma Physics Group (X-1)
//   Applied Physics Division
//   Los Alamos National Lab
// March/April 2004 - conversion of v4pic2 fdtd test code into an input deck

begin_globals {
};

begin_initialization {
  if( nproc()!=1 ) {
    sim_log( "This test case requires one processor" );
    abort(0);
  }

  num_step = 256;
  status_interval = 1;
  define_units( 1, 1 );
  define_timestep( 0.95/sqrt(3.0) );
  define_absorbing_grid( -32, -32, -32, // Low corner
                          32,  32,  32, // High corner
                          64, 64, 64,   // Resolution
                          1, 1, 1,      // Topology
                          absorb_particles ); 

  // Space is by default filled with first material defined
  define_material("vacuum",1.0,1.0,0.0);

  // Create the field array
  define_field_array( NULL, 0.01 );
}

begin_diagnostics {
  if( step()%8==0 ) dump_fields( "fields" );
}

begin_particle_injection {
}

begin_current_injection {
  if( rank()==0 ) field(33,33,32).jfz = sin((2.0*M_PI/16.0)*step());
  if( rank()==0 ) field(33,33,33).jfz = sin((2.0*M_PI/16.0)*step());
}

begin_field_injection {
}

begin_particle_collisions {
  // No collisions for this simulation
}


