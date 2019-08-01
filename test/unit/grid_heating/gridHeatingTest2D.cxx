//==============================================================================
/* Grid heating test case in 2D

   This deck is made for quickly testing code changes by testing the grid
   heating rate. It initializes a 2D box of plasma with periodic boundary
   conditions and lets it evolve for for 2 picoseconds.  With unchanged
   parameters, the electrons will enter a linear heating regime after about 0.7
   picoseconds.  The protons will not heat up enough to enter a linear regime.

   This is not a test of the correctness of the physics in the code, i.e., it
   does not compare to an analytic physical result and asess if the code gives
   an answer close to that.  Instead, it tests (more specifically,
   gridHeatingTest.py tests) if the code output has changed appreciably from a
   reference.

   This deck is very heavily modified from a short pulse deck origionally
   written by Brian J. Albright, 2005.

   Written by Scott V. Luedtke, XCP-6, July 31, 2019

*/
//==============================================================================
//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#define CATCH_CONFIG_RUNNER // We will provide a custom main
#include "catch.hpp"

// TODO: this import may ultimately be a bad idea, but it lets you paste an input deck in...

#include "deck/wrapper.h"

#include "src/species_advance/species_advance.h"
#include "src/vpic/vpic.h"


begin_globals {
  double omega_0;                // w0/wpe
  double vthe;                   // vthe/c   <- these are needed to make movie files
  double vthi_I2;                // vthi_I2/c
  int energies_interval;        // how frequently to dump energies
  int    field_interval;         // how frequently to dump field built-in diagnostic
  int    mobile_ions;	         // flag: 0 if ions are not to be pushed
  int    I1_present;             // flag nonzero when H ions are present. 
  int    I2_present;             // flag nonzero when He ions are present.  

  int    eparticle_interval;
  int    I2particle_interval;
  int    load_particles;         // Flag to turn off particle load for testing wave launch. 

  // Dump parameters for standard VPIC output formats
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hI2dParams;
  std::vector<DumpParameters *> outputParams;

};

// Density function helpers for defining how dense the plasma is as a function of position.  Uses SI units to avoid confusion.  The initialization function is responsible for feeding SI units.  Returns a value between 0 and 1 which corresponds to the percetnage of n_0 the density should be for position (x,y,z)
double density_func(double x, double y, double z, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax){
    // Uniform plasma to test grid heating
    return 1;
}

begin_initialization {

  // do not change parameters in this block:
  double elementary_charge  = 4.8032e-10;             // stat coulomb
  double elementary_charge2 = elementary_charge * elementary_charge;    
  double speed_of_light     = 2.99792458e10;          // cm/sec
  double m_e                = 9.1094e-28;             // g
  double k_boltz            = 1.6022e-12;             // ergs/eV
  double mec2               = m_e*speed_of_light*speed_of_light/k_boltz;
  double mpc2               = mec2*1836.0;
  double eps0               = 1;

  double cfl_req            = 0.98;        // How close to Courant should we try to run
  double damp               = 0.0;         // Level of radiation damping

  double n_e_over_n_crit    = 1e2;       // n_e/n_crit in solid slab
  double vacuum_wavelength  = 1000 * 1e-7; // 1 micron laser (not used)
  double delta = (vacuum_wavelength/(2.0*M_PI))/sqrt(n_e_over_n_crit); // c/wpe

  double nx = 10;
  double ny = 10;
  double nz = 1;

  double box_size_x         = nx*0.5*delta;//0.02 * 1e-4;  // microns
  double box_size_y         = ny*0.5*delta;//0.02 * 1e-4;  // microns (ignored if 1d or 2d out of plane)
  double box_size_z         = nz*0.5*delta;//0.02 * 1e-4;  // microns (ignored if 1d or 2d out of plane)

  double FWHM               = 3.0 * 1e-4;         // from Wilks et al.
  float dfrac               = 0.0;                // fraction of charge density for n_Al12/ne

  int load_particles = 1;         // Flag to turn off particle load for testing wave launch. William Daughton.
  //????????????????????????????????????????????????????????????????????
  double nppc        = 64   ;          // Average number of macro particles/cell of each species
  int mobile_ions         = 1;           // whether or not to push ions

  //????????????????????????????????????????????????????????????????????
  double quota = 3.9;             // Run quota in hours.  
  double quota_sec = quota*3600;  // Run quota in seconds. 

  double A_I2     = 1.0;             // proton
  double Z_I2     = 1;
  double mic2_I2  = mpc2*A_I2;
  double mime_I2  = mic2_I2/mec2;

  double t_e = 1e4;                   // electron temp, eV
  double t_i = 1e4;                      // ion temp, eV
  double uthe    = sqrt(t_e/mec2);    // vthe/c
  double uthi_I2 = sqrt(t_i/mic2_I2); // vthi/c

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

  double t_stop = 2000 * 1e-15*speed_of_light/delta;                // runtime in 1/omega_pe

  // Diagnostics intervals.  
  int energies_interval = 100;
  int field_interval    = 100000;
//?????????????????????????????????????????????????????????????????????????????
  int I2particle_interval = field_interval ;
  int eparticle_interval = field_interval;

//??????????????????????????????????????????????????????????????????????????????????????????

  double omega_0   = sqrt(1.0/n_e_over_n_crit);  // w0/wpe

  double delta_0 = delta*sqrt(n_e_over_n_crit); // c/w0

  double Ne    = nppc*nx*ny*nz;             // Number of macro electrons in box
  Ne = trunc_granular(Ne, nproc());         // Make Ne divisible by number of processors       
  double Ni    = Ne;
  double Npe   = Lx*Ly*Lz;                  // Number of physical electrons in box, wpe = 1
  double qe    = -Npe/Ne;                   // Charge per macro electron
  double qi_I1 = -dfrac*qe;                 // Charge per macro ion of type 1. Note that species
  double qi_I2 = -(1.0-dfrac)*qe;           // I2 and I1 are separate from one another in the loading.
  if ( load_particles==0 ) Ne=Ni=98.7654;   // A weird number to signal that load_particles turned off. 

  int I1_present=0;
  int I2_present=1;


  // Print stuff that I need for plotters and such, and with enough sig figs!
  // Be very careful modifying this.  Plotters depend on explicit locations of
  // some of these numbers.  Generally speaking, add lines at the end only.
  if(rank() == 0){
    FILE * out;
    out = fopen("params.txt", "w");
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
    fprintf(out, "%d   Field Interval\n", field_interval);
    fprintf(out, "%d   Tracer Interval\n", 0);
    fprintf(out, "%d   Number of tracers per species\n", 0);
    fprintf(out, "%d   Number of tracer species\n", 0);
    fprintf(out, "%d   Number of variables per tracer (possibly wrong)\n", 0);
    fprintf(out, "%d   Number of steps in the entire simulation\n", int(t_stop/(dt)));
    fprintf(out, "%d   Number of ranks\n", nproc());
    fprintf(out, "%.14e   Spec max for electron, code units (gamma-1)s\n", 0.);
    fprintf(out, "%.14e   Spec max for I2, code units (gamma-1)\n", 0.);
    fprintf(out, "%d   Number of bins in the spectra\n", 0);
    fprintf(out, "%d   Spectra Interval\n", 0);
    fprintf(out, "%d   This is my rank\n", rank());
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
  sim_log("* Charge/macro I2 =              "<<qi_I2);
  sim_log("* Charge/macro I1 =              "<<qi_I1);
  sim_log("* particles_alloc =              "<<particles_alloc);
  sim_log("* Average particles/processor:   "<<Ne/nproc());
  sim_log("* Average particles/cell:        "<<nppc);
  sim_log("* Do we have mobile ions?        "<<(mobile_ions ? "Yes" : "No"));
  sim_log("* Omega_0, Omega_pe:             "<<(omega_0)<<" "<<1);
  sim_log("* Plasma density, ne/nc:         "<<n_e<<" "<<n_e_over_n_crit);
  sim_log("* T_e, T_i, m_e, m_i_I1, m_i_I2:  "<<t_e<<" "<<t_i<<" "<<1<<" "<<1.<<" "<<mime_I2);
  sim_log("* Radiation damping:             "<<damp);
  sim_log("* Fraction of courant limit:     "<<cfl_req);
  sim_log("* vthe/c:                        "<<uthe);
  sim_log("* vthi_I1/c, vth_I2/c:            "<<0.<<" "<<uthi_I2);
  sim_log("* energies_interval:             "<<energies_interval);
  sim_log("* field_interval:                "<<field_interval);
  sim_log("*********************************");

  // SETUP HIGH-LEVEL SIMULATION PARMETERS
  // FIXME : proper normalization in these units for: xfocus, ycenter, zcenter, waist
  sim_log("Setting up high-level simulation parameters. "); 
  num_step             = int(t_stop/(dt)); 
  status_interval      = 200; 
//????????????????????????????????????????????????????????????????????????????????
//sync_shared_interval = status_interval/10;
//clean_div_e_interval = status_interval/10;
//clean_div_b_interval = status_interval/10;
  sync_shared_interval = status_interval;
  clean_div_e_interval = status_interval;
  clean_div_b_interval = status_interval;
  verbose = 0;
  global->energies_interval        = energies_interval;
  global->field_interval           = field_interval; 
  global->vthe                     = uthe;     // c=1
  global->vthi_I2                  = uthi_I2;  // c=1
  global->omega_0                  = omega_0;
  global->mobile_ions              = mobile_ions; 

  global->eparticle_interval          = eparticle_interval; 
  global->I2particle_interval           = I2particle_interval; 
  global->load_particles           = load_particles; 

  global->I1_present           = I1_present; 
  global->I2_present           = I2_present; 

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
  double max_local_np_i1            = max_local_np_e;
  double max_local_np_i2            = max_local_np_e;
  double max_local_nm_e            = max_local_np_e / 1.;
  double max_local_nm_i1            = max_local_nm_e;
  double max_local_nm_i2            = max_local_nm_e;
  species_t * electron = define_species("electron", -1, 1, max_local_np_e, max_local_nm_e, 20, 1);
  species_t *ion_I1, *ion_I2;
  if ( mobile_ions ) {
  sim_log("Setting up ions. ");
    if ( I2_present  ) ion_I2 = define_species("I2", Z_I2, mime_I2, max_local_np_i2, max_local_nm_i2, 80, 1);
  }


  // SETUP THE MATERIALS
  // Note: the semantics of Kevin's boundary handler for systems where the particles are absorbed
  //       requires that the field array be defined prior to setting up the custom boundary handlers

  sim_log("Setting up materials. "); 
  define_material( "vacuum", 1 );
  // define_material( "impermeable_vacuum", 1 );
//define_material( "impermeable_vacuum_xr", 1 );
  define_field_array( NULL, damp ); 

  // LOAD PARTICLES

  // Load particles using rejection method (p. 290 Num. Recipes in C 2ed, Press et al.)  

  if ( load_particles!=0 ) {
    sim_log( "Loading particles" );
  
    seed_entropy( 2392 );  // Kevin said it should be this way
    // Fast load of particles
    double xmin = grid->x0;
    double xmax = (grid->x0+grid->nx*grid->dx);
    double ymin = grid->y0;
    double ymax = (grid->y0+grid->ny*grid->dy);
    double zmin = grid->z0;
    double zmax = (grid->z0+grid->nz*grid->dz);
 
    // Convert to SI for the density function
    double lengthtoSI = box_size_x/(1e2*Lx);
    double SItolength = 1./lengthtoSI;
    // This must match how you set up the grid.
    double minx = 0.*lengthtoSI;
    double maxx = Lx*lengthtoSI;
    double miny = -0.5*Ly*lengthtoSI;
    double maxy = 0.5*Ly*lengthtoSI;
    double minz = -0.5*Lz*lengthtoSI;
    double maxz = 0.5*Lz*lengthtoSI;

    
    repeat( (Ne)/(topology_x*topology_y*topology_z) ) {
      double x = uniform( rng(0), xmin, xmax );
      double y = uniform( rng(0), ymin, ymax );   
      double z = uniform( rng(0), zmin, zmax );

      // Rejection method, based on user-defined density function
      if ( uniform( rng(0), 0, 1 ) < density_func(x*lengthtoSI,y*lengthtoSI,z*lengthtoSI,minx,maxx,miny,maxy,minz,maxz) ) {
          inject_particle( electron, x, y, z,
                         normal( rng(0), 0, uthe ),
                         normal( rng(0), 0, uthe ),
                         normal( rng(0), 0, uthe ), fabs(qe), 0, 0 );
          inject_particle( ion_I2, x, y, z,
                         normal( rng(0), 0, uthi_I2 ),
                         normal( rng(0), 0, uthi_I2 ),
                         normal( rng(0), 0, uthi_I2 ), fabs(qi_I2)/Z_I2, 0, 0 );
      }
    }
 
 } // if load_particles

 /*--------------------------------------------------------------------------
  * New dump definition
  *------------------------------------------------------------------------*/

 /*--------------------------------------------------------------------------
  * Set data output format
  * 
  * This option allows the user to specify the data format for an output
  * dump.  Legal settings are 'band' and 'band_interleave'.  Band-interleave
  * format is the native storage format for data in VPIC.  For field data,
  * this looks something like:
  * 
  *   ex0 ey0 ez0 div_e_err0 cbx0 ... ex1 ey1 ez1 div_e_err1 cbx1 ...
  *   
  * Banded data format stores all data of a particular state variable as a
  * contiguous array, and is easier for ParaView to process efficiently. 
  * Banded data looks like:
  * 
  *   ex0 ex1 ex2 ... exN ey0 ey1 ey2 ...
  *   
  *------------------------------------------------------------------------*/
  sim_log("Setting up hydro and field diagnostics.");

  global->fdParams.format = band;
  sim_log ( "Field output format          : band" );

  global->hedParams.format = band;
  sim_log ( "Electron hydro output format : band" );

  global->hI2dParams.format = band;
  sim_log ( "I2 hydro output format   : band" );


 /*--------------------------------------------------------------------------
  * Set stride
  * 
  * This option allows data down-sampling at output.  Data are down-sampled
  * in each dimension by the stride specified for that dimension.  For
  * example, to down-sample the x-dimension of the field data by a factor
  * of 2, i.e., half as many data will be output, select:
  * 
  *   global->fdParams.stride_x = 2;
  *
  * The following 2-D example shows down-sampling of a 7x7 grid (nx = 7,
  * ny = 7.  With ghost-cell padding the actual extents of the grid are 9x9.
  * Setting the strides in x and y to equal 2 results in an output grid of
  * nx = 4, ny = 4, with actual extents 6x6.
  *
  * G G G G G G G G G
  * G X X X X X X X G
  * G X X X X X X X G         G G G G G G
  * G X X X X X X X G         G X X X X G
  * G X X X X X X X G   ==>   G X X X X G
  * G X X X X X X X G         G X X X X G
  * G X X X X X X X G         G X X X X G
  * G X X X X X X X G         G G G G G G
  * G G G G G G G G G
  *
  * Note that grid extents in each dimension must be evenly divisible by
  * the stride for that dimension:
  *
  *   nx = 150;
  *   global->fdParams.stride_x = 10; // legal -> 150/10 = 15
  *
  *   global->fdParams.stride_x = 8; // illegal!!! -> 150/8 = 18.75
  *------------------------------------------------------------------------*/

  // Strides for field and hydro arrays.  Note that here we have defined them 
  // the same for fields and all hydro species; if desired, we could use different
  // strides for each.   Also note that strides must divide evenly into the number 
  // of cells in a given domain. 

  // Define strides and test that they evenly divide into grid->nx, ny, nz
//int stride_x = 2, stride_y = 4, stride_z = 4; 
//??????????????????????????????????????????????????????????????????????
//int stride_x = 4, stride_y = 4, stride_z = 4;
//int stride_x = 4, stride_y = 3, stride_z = 3;
  int stride_x = 1, stride_y = 1, stride_z = 1;
//int stride_x = 1, stride_y = 3, stride_z = 3;
  if ( int(grid->nx)%stride_x ) ERROR(("Stride doesn't evenly divide grid->nx."));
  if ( int(grid->ny)%stride_y ) ERROR(("Stride doesn't evenly divide grid->ny."));
  if ( int(grid->nz)%stride_z ) ERROR(("Stride doesn't evenly divide grid->nz."));

  //----------------------------------------------------------------------
  // Fields

  // relative path to fields data from global header
  sprintf(global->fdParams.baseDir, "field");

  // base file name for fields output
  sprintf(global->fdParams.baseFileName, "fields");

  // set field strides
  global->fdParams.stride_x = stride_x;
  global->fdParams.stride_y = stride_y;
  global->fdParams.stride_z = stride_z;
  sim_log ( "Fields x-stride " << global->fdParams.stride_x );
  sim_log ( "Fields y-stride " << global->fdParams.stride_y );
  sim_log ( "Fields z-stride " << global->fdParams.stride_z );

  // add field parameters to list
  global->outputParams.push_back(&global->fdParams);

  //----------------------------------------------------------------------
  // Electron hydro

  // relative path to electron species data from global header
  sprintf(global->hedParams.baseDir, "ehydro");

  // base file name for fields output
  sprintf(global->hedParams.baseFileName, "e_hydro");

  // set electron hydro strides
  global->hedParams.stride_x = stride_x;
  global->hedParams.stride_y = stride_y;
  global->hedParams.stride_z = stride_z;
  sim_log ( "Electron species x-stride " << global->hedParams.stride_x );
  sim_log ( "Electron species y-stride " << global->hedParams.stride_y );
  sim_log ( "Electron species z-stride " << global->hedParams.stride_z );

  // add electron hydro parameters to list
  global->outputParams.push_back(&global->hedParams);

  //----------------------------------------------------------------------
  // ion I2 hydro

  // relative path to electron species data from global header
  sprintf(global->hI2dParams.baseDir, "I2hydro");

  // base file name for fields output
  sprintf(global->hI2dParams.baseFileName, "I2_hydro");

  // set helium hydro strides
  global->hI2dParams.stride_x = stride_x;
  global->hI2dParams.stride_y = stride_y;
  global->hI2dParams.stride_z = stride_z;
  sim_log ( "Ion species x-stride " << global->hI2dParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hI2dParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hI2dParams.stride_z );

  // add helium hydro parameters to list
  global->outputParams.push_back(&global->hI2dParams);

 /*-----------------------------------------------------------------------
  * Set output fields
  *
  * It is now possible to select which state-variables are output on a
  * per-dump basis.  Variables are selected by passing an or-list of
  * state-variables by name.  For example, to only output the x-component
  * of the electric field and the y-component of the magnetic field, the
  * user would call output_variables like:
  *
  *   global->fdParams.output_variables( ex | cby );
  *
  * NOTE: OUTPUT VARIABLES ARE ONLY USED FOR THE BANDED FORMAT.  IF THE
  * FORMAT IS BAND-INTERLEAVE, ALL VARIABLES ARE OUTPUT AND CALLS TO
  * 'output_variables' WILL HAVE NO EFFECT.
  *
  * ALSO: DEFAULT OUTPUT IS NONE!  THIS IS DUE TO THE WAY THAT VPIC
  * HANDLES GLOBAL VARIABLES IN THE INPUT DECK AND IS UNAVOIDABLE.
  *
  * For convenience, the output variable 'all' is defined:
  *
  *   global->fdParams.output_variables( all );
  *------------------------------------------------------------------------*/
 /* CUT AND PASTE AS A STARTING POINT
  * REMEMBER TO ADD APPROPRIATE GLOBAL DUMPPARAMETERS VARIABLE

   output_variables( all );

   output_variables( electric | div_e_err | magnetic | div_b_err |
                     tca      | rhob      | current  | rhof |
                     emat     | nmat      | fmat     | cmat );

   output_variables( current_density  | charge_density |
                     momentum_density | ke_density     | stress_tensor );
  */

  //global->fdParams.output_variables( all );
  global->fdParams.output_variables( electric | magnetic );

  //global->hedParams.output_variables( all );
  //global->hedParams.output_variables( current_density | momentum_density );
  global->hedParams.output_variables(  current_density  | charge_density |
                                       momentum_density | ke_density |
                                       stress_tensor );
  global->hI2dParams.output_variables( current_density  | charge_density |
                                       momentum_density | ke_density |
                                       stress_tensor );

 /*--------------------------------------------------------------------------
  * Convenience functions for simlog output
  *------------------------------------------------------------------------*/
  char varlist[256];

  create_field_list(varlist, global->fdParams);
  sim_log ( "Fields variable list: " << varlist );

  create_hydro_list(varlist, global->hedParams);
  sim_log ( "Electron species variable list: " << varlist );

  create_hydro_list(varlist, global->hI2dParams);
  sim_log ( "I2 species variable list: " << varlist );

 /*------------------------------------------------------------------------*/


  sim_log("*** Finished with user-specified initialization ***"); 

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


begin_diagnostics {
//int mobile_ions=global->mobile_ions, 
//    I1_present=global->I1_present,
//    I2_present=global->I2_present;

  //if ( step()%1==0 ) sim_log("Time step: "<<step()); 

# define should_dump(x) \
  (global->x##_interval>0 && remainder(step(),global->x##_interval)==0)

  if ( step()==0 ) {
    // A grid dump contains all grid parameters, field boundary conditions,
    // particle boundary conditions and domain connectivity information. This
    // is stored in a binary format. Each rank makes a grid dump
    //dump_grid("grid");

    // A materials dump contains all the materials parameters. This is in a
    // text format. Only rank 0 makes the materials dump
    //dump_materials("materials");

    // A species dump contains the physics parameters of a species. This is in
    // a text format. Only rank 0 makes the species dump
    //dump_species("species");

    if ( rank()==0 ) {
    dump_mkdir("rundata_2d");
    dump_mkdir("field");
    dump_mkdir("ehydro");
    dump_mkdir("I2hydro");
    dump_mkdir("restart");


    global_header("global", global->outputParams);

    } // if 

  }


  // energy in various fields/particles 
  if( should_dump(energies) ) {
            dump_energies( "rundata_2d/energies", step() ==0 ? 0 : 1 );
  } //if

  if ( should_dump(field) ) {
    field_dump( global->fdParams );

    if ( global->load_particles ) {
      hydro_dump( "electron", global->hedParams );
      if ( global->mobile_ions ) {
        if ( global->I2_present ) hydro_dump( "I2", global->hI2dParams );
      }
    }
  }
}


begin_particle_injection {
  // No particle injection for this simulation
}


begin_current_injection {
  // No current injection for this simulation
}

begin_field_injection { 
} 

begin_particle_collisions {
  // No particle collisions for this simulation
}

TEST_CASE( "grid heating test in 2d", "[grid_heating]" )
{
    // Init and run sim
    vpic_simulation simulation = vpic_simulation();

    // TODO: We should do this in a safer manner
    simulation.initialize( 0, NULL );

    while( simulation.advance() );

    simulation.finalize();

    if( world_rank==0 ) log_printf( "normal exit\n" );

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
