//==============================================================================
/* Grid heating test case for electros in 2D

   This deck is made for quickly testing code changes by testing the grid
   heating rate. It initializes a 2D box of electrons with periodic boundary
   conditions and lets it evolve for for 0.5 picoseconds.  The electrons are
   started very hot in an approximately linear heating regime.  Ions are
   excluded for speed.

   This is not a test of the correctness of the physics in the code, i.e., it
   does not compare to an analytic physical result and asess if the code gives
   an answer close to that.  Instead, it tests (more specifically,
   gridHeatingTestElec.py tests) if the code output has changed appreciably
   from a reference.

   This deck is very heavily modified from a short pulse deck origionally
   written by Brian J. Albright, 2005.

   Written by Scott V. Luedtke, XCP-6, August 15, 2019

*/
//==============================================================================

//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#define CATCH_CONFIG_RUNNER // We will provide a custom main
#include "catch.hpp"

// TODO: this import may ultimately be a bad idea, but it lets you paste an input deck in...
#include "deck/wrapper.h"

begin_globals {
    int energies_interval;        // how frequently to dump energies
};

begin_diagnostics {

    //if ( step()%1==0 ) sim_log("Time step: "<<step());

# define should_dump(x) \
    (global->x##_interval>0 && remainder(step(),global->x##_interval)==0)

    if ( step()==0 ) {
        if ( rank()==0 ) {
            dump_mkdir("heating_rundata");
        } // if
    }

    // energy in various fields/particles
    if( should_dump(energies) ) {
        dump_energies( "heating_rundata/energies", step() ==0 ? 0 : 1 );
    } //if

}

begin_initialization {

    // do not change parameters in this block:
    double elementary_charge  = 4.8032e-10;             // stat coulomb
    double elementary_charge2 = elementary_charge * elementary_charge;
    double speed_of_light     = 2.99792458e10;          // cm/sec
    double m_e                = 9.1094e-28;             // g
    double k_boltz            = 1.6022e-12;             // ergs/eV
    double mec2               = m_e*speed_of_light*speed_of_light/k_boltz;
    double eps0               = 1;

    double cfl_req            = 0.98;        // How close to Courant should we try to run
    double damp               = 0.0;         // Level of radiation damping

    double n_e_over_n_crit    = 1e2;       // n_e/n_crit in solid slab
    double vacuum_wavelength  = 1e3 *1e-7; // used for normalization
    double delta = (vacuum_wavelength/(2.0*M_PI))/sqrt(n_e_over_n_crit); // c/wpe

    double nx = 10;
    double ny = 10;
    double nz = 1;

    double box_size_x         = nx*0.5*delta;// microns
    double box_size_y         = ny*0.5*delta;// microns
    double box_size_z         = nz*0.5*delta;// microns

    int load_particles = 1;         // Flag to turn off particle load for testing wave launch. William Daughton.
    //????????????????????????????????????????????????????????????????????
    double nppc        = 64   ;          // Average number of macro particles/cell of each species

    double t_e = 20e4;                   // electron temp, eV
    double uthe    = sqrt(t_e/mec2);    // vthe/c

    double n_e = speed_of_light*speed_of_light*m_e/(4.0*M_PI*elementary_charge2*delta*delta);
    // n_e is used for emax only
    double debye = uthe*delta;

    //??????????????????????????????????????????????????????????????????????????
    int topology_x = 1;
    int topology_y = 1;
    int topology_z = 1;

    double cell_size_t = box_size_y/(debye*ny);  // in debye
    double cell_size_tz= box_size_z/(debye*nz);
    double cell_size_l = box_size_x/(debye*nx);

    double hy = debye*cell_size_t/delta;
    double hz = debye*cell_size_tz/delta;
    double hx = debye*cell_size_l/delta;

    double Lx = nx*hx;          // in c/wpe
    double Ly = ny*hy;
    double Lz = nz*hz;

    double particles_alloc = nppc*ny*nz*nx;

    double dt = cfl_req*courant_length(Lx, Ly, Lz, nx, ny, nz);
    double dt_courant = dt;

    double t_stop = 500. * 1e-15*speed_of_light/delta;// runtime in 1/omega_pe

    // Diagnostics intervals.
    int energies_interval = 100;
    global->energies_interval = energies_interval;

    double omega_0   = sqrt(1.0/n_e_over_n_crit);  // w0/wpe

    double delta_0 = delta*sqrt(n_e_over_n_crit); // c/w0

    double Ne    = nppc*nx*ny*nz;             // Number of macro electrons in box
    Ne = trunc_granular(Ne, nproc());         // Make Ne divisible by number of processors
    double Npe   = Lx*Ly*Lz;                  // Number of physical electrons in box, wpe = 1
    double qe    = -Npe/Ne;                   // Charge per macro electron

    // Print stuff that I need for plotters and such, and with enough sig figs!
    // Be very careful modifying this.  Plotters depend on explicit locations of
    // some of these numbers.  Generally speaking, add lines at the end only.
    if(rank() == 0){
        FILE * out;
        out = fopen("heating_params.txt", "w");
        fprintf(out, "# Parameter file used for plotters.\n");
        fprintf(out, "%.14e   Time step (dt), code units\n", dt);
        fprintf(out, "%.14e   Laser wavelength, SI\n", vacuum_wavelength*1e-2);
        fprintf(out, "%.14e   Ratio of electron to critical density\n", n_e_over_n_crit);
        fprintf(out, "%d   Number of cells in x\n", int(nx));
        fprintf(out, "%d   Number of cells in y\n", int(ny));
        fprintf(out, "%d   Number of cells in z\n", int(nz));
        fprintf(out, "%.14e   Box size x, microns\n", box_size_x*1e4);
        fprintf(out, "%.14e   Box size y, microns\n", box_size_y*1e4);
        fprintf(out, "%.14e   Box size z, microns\n", box_size_z*1e4);
        fclose(out);
    }

    // PRINT SIMULATION PARAMETERS

    sim_log("***** Simulation parameters *****");
    sim_log("* Processors:                    "<<nproc());
    sim_log("* dt_courant"<<dt_courant);
    sim_log("* delta/dx,delta/dz =  "<<1.0/hx<<" "<<1.0/hz);
    sim_log("* Time step, max time, nsteps =  "<<dt<<" "<<t_stop<<" "<<int(t_stop/(dt)));
    sim_log("* Debye length,cell size_l,cell size_t,delta,delta_0 = "<<debye<<" "<<cell_size_l<<" "<<cell_size_t<<" "<<delta<<" "<<delta_0);
    sim_log("* Lx, Ly, Lz =                   "<<Lx<<" "<<Ly<<" "<<Lz);
    sim_log("* nx, ny, nz =                   "<<nx<<" "<<ny<<" "<<nz);
    sim_log("* Charge/macro electron =        "<<qe);
    sim_log("* particles_alloc =              "<<particles_alloc);
    sim_log("* Average particles/processor:   "<<Ne/nproc());
    sim_log("* Average particles/cell:        "<<nppc);
    sim_log("* Omega_0, Omega_pe:             "<<(omega_0)<<" "<<1);
    sim_log("* Plasma density, ne/nc:         "<<n_e<<" "<<n_e_over_n_crit);
    sim_log("* T_e,m_e:  "<<t_e<<" "<<1);
    sim_log("* Radiation damping:             "<<damp);
    sim_log("* Fraction of courant limit:     "<<cfl_req);
    sim_log("* vthe/c:                        "<<uthe);
    sim_log("* energies_interval:             "<<energies_interval);
    sim_log("*********************************");

    // SETUP HIGH-LEVEL SIMULATION PARMETERS
    // FIXME : proper normalization in these units for: xfocus, ycenter, zcenter, waist
    sim_log("Setting up high-level simulation parameters. ");
    num_step             = int(t_stop/(dt));
    status_interval      = 20000;
    //????????????????????????????????????????????????????????????????????????????????
    sync_shared_interval = status_interval/10;
    clean_div_e_interval = status_interval/10;
    clean_div_b_interval = status_interval/10;
    verbose = 0;

    // SETUP THE GRID
    sim_log("Setting up computational grid.");
    grid->dx = hx;
    grid->dy = hy;
    grid->dz = hz;
    grid->dt = dt;
    grid->cvac = 1;
    grid->eps0 = eps0;

    // Partition a periodic box among the processors sliced uniformly in z:
    define_periodic_grid( 0,      -0.5*Ly,  -0.5*Lz,    // Low corner
            Lx,      0.5*Ly,   0.5*Lz,    // High corner
            nx,      ny,       nz,        // Resolution
            topology_x, topology_y, topology_z); // Topology

    sim_log("Setting up electrons. ");

    //???????????????????????????????????????????????????????????????????????
    // How oversized should the particle buffers be in case of non-uniform plasma?
    // The left and right ranks are half filled, and the middle will get to 1.5x
    // filled.  An over_alloc_fac between .5 and 1.5 should drop if not
    // dynamically resizing.
    double over_alloc_fac = 1.3;
    double max_local_np_e            = over_alloc_fac*particles_alloc/nproc();
    double max_local_nm_e            = max_local_np_e / 1.;
    species_t * electron = define_species("electron", -1, 1, max_local_np_e, max_local_nm_e, 0, 1);

    // SETUP THE MATERIALS
    sim_log("Setting up materials. ");
    define_material( "vacuum", 1 );
    define_field_array( NULL, damp );

    // LOAD PARTICLES

    // Load particles using rejection method (p. 290 Num. Recipes in C 2ed, Press et al.)

    if ( load_particles!=0 ) {
        sim_log( "Loading particles" );

        seed_entropy( 2995471 );  // Kevin said it should be this way
        // Fast load of particles
        double xmin = grid->x0;
        double xmax = (grid->x0+grid->nx*grid->dx);
        double ymin = grid->y0;
        double ymax = (grid->y0+grid->ny*grid->dy);
        double zmin = grid->z0;
        double zmax = (grid->z0+grid->nz*grid->dz);

        repeat( (Ne)/(topology_x*topology_y*topology_z) ) {
            double x = uniform( rng(0), xmin, xmax );
            double y = uniform( rng(0), ymin, ymax );
            double z = uniform( rng(0), zmin, zmax );

            // Rejection method, based on user-defined density function
            // Simplified for the uniform plasma used in this test
            if ( uniform( rng(0), 0, 1 ) < 1. ) {
                inject_particle( electron, x, y, z,
                        normal( rng(0), 0, uthe ),
                        normal( rng(0), 0, uthe ),
                        normal( rng(0), 0, uthe ), fabs(qe), 0, 0 );
            }
        }

    } // if load_particles

}

TEST_CASE( "Check if Weibel gives correct energy (within tol)", "[energy]" )
{
    // Init and run sim
    vpic_simulation simulation = vpic_simulation();

    // TODO: We should do this in a safer manner
    simulation.initialize( 0, NULL );

    while( simulation.advance() );

    simulation.finalize();
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
