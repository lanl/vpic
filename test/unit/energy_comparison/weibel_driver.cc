//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#define CATCH_CONFIG_RUNNER // We will provide a custom main
#include "catch.hpp"

// TODO: this import may ultimately be a bad idea, but it lets you paste an input deck in...

#include "deck/wrapper.h"

#include "src/species_advance/species_advance.h"
#include "src/vpic/vpic.h"

#include "compare_energies.h"

begin_globals {
  double energies_interval;
  double fields_interval;
  double ehydro_interval;
  double ihydro_interval;
  double eparticle_interval;
  double iparticle_interval;
  double restart_interval;
};

std::string energy_file_name = "./energies";

std::string energy_gold_file_name = EXPAND_AND_STRINGIFY( GOLD_ENERGY_FILE );

void vpic_simulation::user_diagnostics() {
    dump_energies(energy_file_name.c_str(), 1);
}

begin_initialization {
// AKA:
//void
//vpic_simulation::user_initialization( int num_cmdline_arguments,
                                      //char ** cmdline_argument )
//{
  // At this point, there is an empty grid and the random number generator is
  // seeded with the rank. The grid, materials, species need to be defined.
  // Then the initial non-zero fields need to be loaded at time level 0 and the
  // particles (position and momentum both) need to be loaded at time level 0.

  // Arguments can be passed from the command line to the input deck
  // if( num_cmdline_arguments!=3 ) {
  //   sim_log( "Usage: " << cmdline_argument[0] << " mass_ratio seed" );
  //   abort(0);
  // }
  seed_entropy(1); //seed_entropy( atoi( cmdline_argument[2] ) );

  // Diagnostic messages can be passed written (usually to stderr)
  sim_log( "Computing simulation parameters");

  // Define the system of units for this problem (natural units)
  //double L    = 1; // Length normalization (sheet thickness)
  double de   = 1; // Length normalization (electron inertial length)
  double ec   = 1; // Charge normalization
  double me   = 1; // Mass normalization
  double c    = 1; // Speed of light
  double eps0 = 1; // Permittivity of space

  // Physics parameters
  double mi_me   = 1836; //25; //atof(cmdline_argument[1]); // Ion mass / electron mass
  double vthe = 0.25/sqrt(2.0); //0.0424264068711;       //0.424264068711;       // Electron thermal velocity
  double vthi = 0.25/sqrt(2.0); //0.0424264068711;       //0.424264068711;       // Ion thermal velocity
  double vthex =0.05/sqrt(2.0); //0.0141421356237;      // 0.141421356237;      // Electron thermal velocity in x-direction.
  double vthix =0.05/sqrt(2.0); //0.0141421356237;      // 0.141421356237;Ion thermal velocity in x-direction.

  double n0      = 1.0;    //  Background plasma density
  double b0 = 0.0;         // In plane magnetic field.
  double tauwpe    = 200000;    // simulation wpe's to run

  // Numerical parameters
  double topology_x = nproc();  // Number of domains in x, y, and z
  double topology_y = 1;
  double topology_z = 1;  // For load balance, best to keep "1" or "2" for Harris sheet
  double Lx        = 2.09439510239320; //4.62*de; //6.7*de; //10.0*de;  // How big should the box be in the x direction
  double Ly        = 1; //0.0721875*de;  // How big should the box be in the y direction
  double Lz        = 1; //0.0721875*de;  // How big should the box be in the z direction
  double nx        = 16; //64; //64; //32;    // Global resolution in the x direction
  double ny        = 1;    // Global resolution in the y direction
  double nz        = 1; //32;     // Global resolution in the z direction
  double nppc      = 200; //800; //200; //2048; //1024; //128;    // Average number of macro particles per cell (both species combined!)
  double cfl_req   = 0.99f; //0.99;  // How close to Courant should we try to run
  double wpedt_max = 0.36;  // How big a timestep is allowed if Courant is not too restrictive
  double damp      = 0.0; // Level of radiation damping


  // Derived quantities
  double mi = me*mi_me;             // Ion mass
  double wpe  = c/de;               // electron plasma frequency
  double wpi  = wpe/sqrt(mi_me);    // ion plasma frequency
  double di   = c/wpi;              // ion inertial length

  double hx = Lx/nx;
  double hy = Ly/ny;
  double hz = Lz/nz;

  double Npe = n0*Ly*Lz*Lx;    // Number physical electrons.
  double Npi = Npe;            // Number of physical ions in box
  double Ne  = nppc*nx*ny*nz;  // total macro electrons in box

  Ne = trunc_granular(Ne,nproc());
  double Ni   = Ne;                                   // Total macro ions in box

  double we   = Npe/Ne;                               // Weight of a macro electron
  double wi   = Npi/Ni;                               // Weight of a macro ion


  // Determine the timestep
  double dg = courant_length(Lx,Ly,Lz,nx,ny,nz);      // Courant length
  double dt = cfl_req*dg/c;                           // Courant limited time step
  // printf("in harris.cxx: dt=%.7f\n",  dt);
  // exit(1);
  if( wpe*dt>wpedt_max ) dt=wpedt_max/wpe;            // Override time step if plasma frequency limited

  ////////////////////////////////////////
  // Setup high level simulation parmeters

  num_step             = 700; //4000; // int(tauwpe/(wpe*dt));
  status_interval      = 0; //2000;
  sync_shared_interval = 0; //status_interval;
  clean_div_e_interval = 0; //turn off cleaning (GY)//status_interval;
  clean_div_b_interval = 0; //status_interval; //(GY)

  global->energies_interval  = 1; //status_interval;
  global->fields_interval    = status_interval;
  global->ehydro_interval    = status_interval;
  global->ihydro_interval    = status_interval;
  global->eparticle_interval = status_interval; // Do not dump
  global->iparticle_interval = status_interval; // Do not dump
  global->restart_interval   = status_interval; // Do not dump

  ///////////////////////////
  // Setup the space and time

  // Setup basic grid parameters
  define_units( c, eps0 );
  define_timestep( dt );
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = c;
  //grid->damp = damp;

  // Parition a periodic box among the processors sliced uniformly along y
  // define_periodic_grid( -0.5*Lx, 0, 0,    // Low corner
  //                        0.5*Lx, Ly, Lz,  // High corner
  //                        nx, ny, nz,      // Resolution
  //                        1, nproc(), 1 ); // Topology
  define_periodic_grid(  0, -0.5*Ly, -0.5*Lz,    // Low corner
			  Lx, 0.5*Ly, 0.5*Lz,     // High corner
			  nx, ny, nz,             // Resolution
			  topology_x, topology_y, topology_z); // Topology

  //   printf("in harris.cxx: g->neighbor[6*265]=%jd\n",  grid->neighbor[6*265]);
  // Override some of the boundary conditions to put a particle reflecting
  // perfect electrical conductor on the -x and +x boundaries
  // set_domain_field_bc( BOUNDARY(-1,0,0), pec_fields );
  // set_domain_field_bc( BOUNDARY( 1,0,0), pec_fields );
  // set_domain_particle_bc( BOUNDARY(-1,0,0), reflect_particles );
  // set_domain_particle_bc( BOUNDARY( 1,0,0), reflect_particles );

  define_material( "vacuum", 1 );
  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.

#ifdef ADD_EXTRA_MATERIAL
  define_material( "not-vacuum", 1 );
#endif

  // If you pass NULL to define field array, the standard field array will
  // be used (if damp is not provided, no radiation damping will be used).
  define_field_array( NULL, damp );

  ////////////////////
  // Setup the species

  // Allow 50% more local_particles in case of non-uniformity
  // VPIC will pick the number of movers to use for each species
  // Both species use out-of-place sorting
  // species_t * ion      = define_species( "ion",       ec, mi, 1.5*Ni/nproc(), -1, 40, 1 );
  // species_t * electron = define_species( "electron", -ec, me, 1.5*Ne/nproc(), -1, 20, 1 );
  //species_t *electron = define_species("electron",-ec,me,2.4*Ne/nproc(),-1,25,0);
  //species_t *ion      = define_species("ion",      ec,mi,2.4*Ne/nproc(),-1,25,0);

  species_t *electron = define_species("electron",-ec,me,2.4*Ne/nproc(),-1,0,0); //turn off sorting (GY)
  species_t *ion      = define_species("ion",      ec,mi,2.4*Ne/nproc(),-1,0,0); //(GY)

  ///////////////////////////////////////////////////
  // Log diagnostic information about this simulation

  sim_log( "***********************************************" );
  sim_log ( "mi/me = " << mi_me );
  sim_log ( "tauwpe = " << tauwpe );
  sim_log ( "num_step = " << num_step );
  sim_log ( "Lx/di = " << Lx/di );
  sim_log ( "Lx/de = " << Lx/de );
  sim_log ( "Ly/di = " << Ly/di );
  sim_log ( "Ly/de = " << Ly/de );
  sim_log ( "Lz/di = " << Lz/di );
  sim_log ( "Lz/de = " << Lz/de );
  sim_log ( "nx = " << nx );
  sim_log ( "ny = " << ny );
  sim_log ( "nz = " << nz );
  sim_log ( "damp = " << damp );
  sim_log ( "courant = " << c*dt/dg );
  sim_log ( "nproc = " << nproc ()  );
  sim_log ( "nppc = " << nppc );
  sim_log ( " b0 = " << b0 );
  sim_log ( " di = " << di );
  sim_log ( " Ne = " << Ne );
  sim_log ( "total # of particles = " << 2*Ne );
  sim_log ( "dt*wpe = " << wpe*dt );
  sim_log ( "dx/de = " << Lx/(de*nx) );
  sim_log ( "dy/de = " << Ly/(de*ny) );
  sim_log ( "dz/de = " << Lz/(de*nz) );
  sim_log ( "dx/debye = " << (Lx/nx)/(vthe/wpe)  );
  sim_log ( "n0 = " << n0 );
  sim_log ( "vthi/c = " << vthi/c );
  sim_log ( "vthe/c = " << vthe/c );
  sim_log( "" );

  ////////////////////////////
  // Load fields and particles

  // sim_log( "Loading fields" );

  // set_region_field( everywhere, 0, 0, 0,                    // Electric field
  //                   0, -sn*b0*tanh(x/L), cs*b0*tanh(x/L) ); // Magnetic field
  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specied as logical equations (i.e. x>0 && x+y<2)

  sim_log( "Loading particles" );

  // Do a fast load of the particles
  //seed_rand( rng_seed*nproc() + rank() );  //Generators desynchronized
  double xmin = grid->x0 , xmax = grid->x0+(grid->dx)*(grid->nx);
  double ymin = grid->y0 , ymax = grid->y0+(grid->dy)*(grid->ny);
  double zmin = grid->z0 , zmax = grid->z0+(grid->dz)*(grid->nz);

  sim_log( "-> Uniform Bi-Maxwellian" );

  double n1,n2,n3;

  repeat ( Ne/nproc() ) {

      double x = uniform( rng(0), xmin, xmax );
      double y = uniform( rng(0), ymin, ymax );
      double z = uniform( rng(0), zmin, zmax );
      n1 = normal(rng(0),0,vthex);
      n2 = normal(rng(0),0,vthe );
      n3 = normal(rng(0),0,vthe );

      inject_particle( electron, x, y, z,
              n1,
              n2,
              n3,we, 0, 0);

      n1 = normal(rng(0),0,vthix);
      n2 = normal(rng(0),0,vthi );
      n3 = normal(rng(0),0,vthi );

      inject_particle( ion, x, y, z,
              n1,
              n2,
              n3,wi, 0 ,0 );

  }

  sim_log( "Finished loading particles" );

  //exit(1);

  // Upon completion of the initialization, the following occurs:
  // - The synchronization error (tang E, norm B) is computed between domains
  //   and tang E / norm B are synchronized by averaging where discrepancies
  //   are encountered.
  // - The initial divergence error of the magnetic field is computed and
  //   one pass of cleaning is done (for good measure)
  // - The bound charge density necessary to give the simulation an initially
  //   clean divergence e is computed.
  // - The particle momentum is uncentered from u_0 to u_{-1/2}
  // - The user diagnostics are called on the initial state
  // - The physics loop is started
  //
  // The physics loop consists of:
  // - Advance particles from x_0,u_{-1/2} to x_1,u_{1/2}
  // - User particle injection at x_{1-age}, u_{1/2} (use inject_particles)
  // - User current injection (adjust field(x,y,z).jfx, jfy, jfz)
  // - Advance B from B_0 to B_{1/2}
  // - Advance E from E_0 to E_1
  // - User field injection to E_1 (adjust field(x,y,z).ex,ey,ez,cbx,cby,cbz)
  // - Advance B from B_{1/2} to B_1
  // - (periodically) Divergence clean electric field
  // - (periodically) Divergence clean magnetic field
  // - (periodically) Synchronize shared tang e and norm b
  // - Increment the time step
  // - Call user diagnostics
  // - (periodically) Print a status message
}

TEST_CASE( "Check if Weibel gives correct energy (within tol)", "[energy]" )
{
    // Before we run this, we must make sure we remove the energy file
    std::ofstream ofs;
    ofs.open(energy_file_name, std::ofstream::out | std::ofstream::trunc);
    ofs.close();

    // Init and run sim
    vpic_simulation simulation = vpic_simulation();

    // TODO: We should do this in a safer manner
    simulation.initialize( 0, NULL );

    while( simulation.advance() );

    simulation.finalize();

    if( world_rank==0 ) log_printf( "normal exit\n" );

    std::cout << "Comparing " << energy_file_name << " to " <<
        energy_gold_file_name << std::endl;

    // Compare energies to make sure everything worked out OK (within 1%)
    const unsigned short e_mask = 0x000E; // = 0b0000001110; <- Non-standard C++11
    const unsigned short b_mask = 0x0070; // = 0b0001110000;    Requires C++14/GNU extensions
    const unsigned short particle_mask = 0x00C0; // = 0b011000000;

    // Test the sum of the e_field
    REQUIRE(
            test_utils::compare_energies(energy_file_name, energy_gold_file_name,
                0.3, e_mask, test_utils::FIELD_ENUM::Sum, 1, "Weibel.e.out")
           );

    // Test the sum of the b_field
    REQUIRE(
            test_utils::compare_energies(energy_file_name, energy_gold_file_name,
                0.03, b_mask, test_utils::FIELD_ENUM::Sum, 1, "Weibel.b.out")
           );


    // Test particle energies individually
    REQUIRE(
            test_utils::compare_energies(energy_file_name, energy_gold_file_name,
                0.01, particle_mask, test_utils::FIELD_ENUM::Sum, 1, "Weibel.p.out")
           );

}

begin_particle_injection {

  // No particle injection for this simulation

}

begin_current_injection {

  // No current injection for this simulation

}

begin_field_injection {

  // No field injection for this simulation

}

begin_particle_collisions{

  // No collisions for this simulation

}

// Manually implement catch main
int main( int argc, char* argv[] )
{

    // Setup
    boot_services( &argc, &argv );

    int result = Catch::Session().run( argc, argv );

    // clean-up...
    halt_services();

    return result;
}
