////////////////////////////////////////////////////////////////////////
//
//  Reconnection Problem --> single Force-Free Current Sheet with conductive BC
//
///////////////////////////////////////////////////////////////////////
#include <math.h>

#include "VPICAdaptor.h"

// structure to hold the data for energy diagnostics
struct edata {
  species_id sp_id;         /* species id */
  double     vth;          /* thermal energy */
  char fname[256];        /* file to save data */
};

// naming convention for the hydro dump files
#define HYDRO_FILE_FORMAT "hydro/T.%d/%s.%d.%d"

// Vadim's in-line average
 #define ALLOCATE(A,LEN,TYPE)                                             \
   if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate."));

begin_globals {

  int restart_interval;
  int energies_interval;
  int fields_interval;
  int ehydro_interval;
  int Hhydro_interval;
  int eparticle_interval;
  int Hparticle_interval;
  int quota_check_interval;  //  How frequently to check if quote exceeded

  int rtoggle;             // enables save of last two restart dumps for safety
  double quota_sec;        // Run quota in seconds
  double b0;               // B0
  double bg;               // Guide field
  double topology_x;       // domain topology
  double topology_y;
  double topology_z;

//  Variables for new output format

  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;

  // Variables for the energy diagnostics

  edata ede;                        // parameters for electron species
  edata edi;                        // parameters for ion species
  double emax;                       // maximum energy (in units of me*c**2)
  int nex;                           // number of energy bins

  //Vadim:  modified restart machinary

  int write_restart ;                 // global flag for all to write restart files
  int write_end_restart ;             // global flag for all to write restart files

  //Vadim: E*j variables

  int nsp;                            // number of species
  int dis_nav;                         // number of steps to average over
  int dis_interval;                    // number of steps between outputs
  int dis_iter;                        // iteration count. 0 means we are not avraging at the moment  : initialized in dissipation.cxx
  int dis_begin_int;                  // the first time step of the interval  : initialized in dissipation.cxx
};


begin_initialization {

 // use natural PIC units

 double ec   = 1;         // Charge normalization
 double me   = 1;         // Mass normalization
 double c    = 1;         // Speed of light
 double de   = 1;         // Length normalization (electron inertial length)
 double eps0 = 1;         // Permittivity of space

  double cfl_req   = 0.7;  // How close to Courant should we try to run
  double wpedt_max = 0.2;  // How big a timestep is allowed if Courant is not too restrictive
  int rng_seed     = 1;     // Random number seed increment

  // Physics parameters

  double mi_me   = 1.0;    // Ion mass / electron mass
  double L_di    = 1.0/sqrt(mi_me);  // Sheet thickness / ion inertial length
  double Ti_Te   = 1.0;      // Ion temperature / electron temperature
  double vthe    = 0.2;     //  Electron thermal speed over c
  double wpe_wce = 2.0;      // electron plasma freq / electron cyclotron freq
  double bg = 0.0;           // electron plasma freq / electron cyclotron freq
  double theta   = 0.0;        // B0 = Bx
  double taui    = 5;      // simulation wci's to run

  double quota   = 14.4;   // run quota in hours
  double quota_sec = quota*3600;  // Run quota in seconds

  double pi = 3.1415927;
  double cs   = cos(theta/180.0*pi);
  double sn   = sin(theta/180.0*pi);

  //derived qunatities

  double mi = me*mi_me;                                   // Ion mass
  double vthi = vthe*sqrt(Ti_Te/mi_me);                   // Ion thermal velocity
  double wci  = 1.0/(mi_me*wpe_wce);                      // Ion cyclotron frequency
  double wce  = wci*mi_me;                                // Electron cyclotron freqeuncy
  double wpe  = wce*wpe_wce;                              // electron plasma frequency
  double wpi  = wpe/sqrt(mi_me);                          // ion plasma frequency
  double di   = c/wpi;                                    // ion inertial length
  double L    = L_di*di;                                  // Sheet thickness in c/wpe

  double ion_sort_interval = 25;        //  Injector moments are also updated at this internal
  double electron_sort_interval=25;    //  Injector moments are also updated at this internal

  // Numerical parameters

  double nppc  =  10; // Average number of macro particle per cell per species

  double Lx  = 30.0/sqrt(mi_me)*di; // size of box in x dimension
  double Ly  = 15.0/sqrt(mi_me)*di;     // size of box in y dimension
  double Lz  = 15.0/sqrt(mi_me)*di; // size of box in z dimension

  sim_log("ACB LX " << Lx << " LY " << Ly << " LZ " << Lz );

  double topology_x = 1;  // Number of domains in x, y, and z
  double topology_y = 2;
  double topology_z = 2;

  double nx = 64;
  double ny = 64;
  double nz = 64;

  double hx = Lx/nx;
  double hy = Ly/ny;
  double hz = Lz/nz;

  double b0  = me*c*wce/ec; // Asymptotic magnetic field strength
  double n0  = me*eps0*wpe*wpe/(ec*ec);  // Peak electron (ion) density
  double Ne  = nppc*nx*ny*nz;  // total macro electrons in box
  double Np  = n0*Lx*Ly*Lz;  //  total number of physical electrons
  Ne  = trunc_granular(Ne,nproc()); // Make it divisible by number of processors
  double weight = Np/Ne;
  double qe = -ec*Np/Ne;  // Charge per macro electron
  double qi =  ec*Np/Ne;  // Charge per macro ion
  double Lpert = Lx;   // wavelength of perturbation
  double dbz = 0.02*b0; //  Perturbation in Bz relative to Bo (Only change here)
  double dbx = -dbz*Lpert/(2.0*Lz); // Set Bx perturbation so that div(B) = 0

  // Determine the time step

  double dg = courant_length(Lx,Ly,Lz,nx,ny,nz);        // courant length
  double dt = cfl_req*dg/c;                             // courant limited time step
  if( wpe*dt>wpedt_max) dt=wpedt_max/wpe;               // override timestep if plasma frequency limited

  //  int restart_interval = int(10000.0/(wci*dt));
  int restart_interval = 8000;
  int energies_interval = 100;
  int interval = int(5.0/(wci*dt));
  int fields_interval = interval;
  int ehydro_interval = interval;
  int Hhydro_interval = interval;
  int eparticle_interval = 200000*interval;
  int Hparticle_interval = 200000*interval;
  int quota_check_interval     = 100;

  //  Determine which domains area along the boundaries - Use macro from grid/partition.c

# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {                   \
    int _ix, _iy, _iz;                                                    \
    _ix  = (rank);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(topology_x);   /* iy = iy+gpy*iz */                    \
    _ix -= _iy*int(topology_x);   /* ix = ix */                           \
    _iz  = _iy/int(topology_y);   /* iz = iz */                           \
    _iy -= _iz*int(topology_y);   /* iy = iy */ 	        	  \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \
  } END_PRIMITIVE

  int ix, iy, iz ;
  RANK_TO_INDEX( int(rank()), ix, iy, iz );


  ///////////////////////////////////////////////
  // Setup high level simulation parameters
  num_step             = int(taui/(wci*dt));

  std::cout << "numstep is " << num_step << std::endl;

  status_interval      = 200;
  sync_shared_interval = status_interval/2;
  clean_div_e_interval = status_interval/2;
  clean_div_b_interval = status_interval/2;

  global->restart_interval   = restart_interval;
  global->energies_interval  = energies_interval;
  global->fields_interval    = fields_interval;
  global->ehydro_interval    = ehydro_interval;
  global->Hhydro_interval    = Hhydro_interval;
  global->eparticle_interval = eparticle_interval;
  global->Hparticle_interval = Hparticle_interval;
  global->quota_check_interval     = quota_check_interval;
  global->quota_sec          = quota_sec;

  global->rtoggle            = 0;
  global->b0  = b0;
  global->bg  = bg;

  global->topology_x  = topology_x;
  global->topology_y  = topology_y;
  global->topology_z  = topology_z;


  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the grid

  // Setup basic grid parameters
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = c;
  grid->eps0 = eps0;

  // Define the grid

  define_periodic_grid(  0, -0.5*Ly, -0.5*Lz,    // Low corner
                          Lx, 0.5*Ly, 0.5*Lz,     // High corner
                          nx, ny, nz,             // Resolution
                          topology_x, topology_y, topology_z); // Topology

 // ***** Set Field Boundary Conditions *****

  // sim_log("Conducting fields on X & Z-boundaries");
  if ( iz==0 )            set_domain_field_bc( BOUNDARY(0,0,-1), pec_fields );
  if ( iz==topology_z-1 ) set_domain_field_bc( BOUNDARY( 0,0,1), pec_fields );

  if ( iz==0 )            set_domain_particle_bc( BOUNDARY(0,0,-1), reflect_particles );
  if ( iz==topology_z-1 ) set_domain_particle_bc( BOUNDARY(0,0,1), reflect_particles );

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup materials

  sim_log("Setting up materials. ");

  define_material( "vacuum", 1 );

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.

////////////////////////////////////////////////////////////////////////////////////////////
//  Finalize Field Advance

  sim_log("Finalizing Field Advance");

  define_field_array(NULL);

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the species

  sim_log("Setting up species. ");
  double nmax = 2.0*Ne/nproc();
  double nmovers = 0.1*nmax;
  double sort_method = 1;   //  0=in place and 1=out of place
  species_t *electron = define_species("electron",-ec, me, nmax, nmovers, electron_sort_interval, sort_method);
  species_t *ion = define_species("ion", ec, mi, nmax, nmovers, ion_sort_interval, sort_method);

  ///////////////////////////////////////////////////
  // Log diagnostic information about this simulation

  sim_log( "***********************************************" );
  sim_log("* Topology:                       "<<topology_x<<" "<<topology_y<<" "<<topology_z);
  sim_log ( "L_di   = " << L_di );
  sim_log ( "Ti/Te = " << Ti_Te ) ;
  sim_log ( "wpe/wce = " << wpe_wce );
  sim_log ( "mi/me = " << mi_me );
  sim_log ( "theta = " << theta );
  sim_log ( "taui = " << taui );
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
  sim_log ( "courant = " << c*dt/dg );
  sim_log ( "nproc = " << nproc()  );
  sim_log ( "nppc = " << nppc );
  sim_log ( " b0 = " << b0 );
  sim_log ( " di = " << di );
  sim_log ( " Ne = " << Ne );
  sim_log ( "total # of particles = " << 2*Ne );
  sim_log ( " qi = " << qi );
  sim_log ( " qe = " << qe );
  sim_log ( "dt*wpe = " << wpe*dt );
  sim_log ( "dt*wce = " << wce*dt );
  sim_log ( "dt*wci = " << wci*dt );
  sim_log ( " energies_interval: " << energies_interval );
  sim_log ( "dx/de = " << Lx/(de*nx) );
  sim_log ( "dy/de = " << Ly/(de*ny) );
  sim_log ( "dz/de = " << Lz/(de*nz) );
  sim_log ( "dx/rhoe = " << (Lx/nx)/(vthe/wce)  );
  sim_log ( "L/debye = " << L/(vthe/wpe)  );
  sim_log ( "dx/debye = " << (Lx/nx)/(vthe/wpe)  );
  sim_log ( "n0 = " << n0 );
  sim_log ( "vthi/c = " << vthi/c );
  sim_log ( "vthe/c = " << vthe/c );


  // Dump simulation information to file "info"
  if (rank() == 0 ) {

    FileIO fp_info;

    if ( ! (fp_info.open("info", io_write)==ok) ) ERROR(("Cannot open file."));

    fp_info.print( "           ***** Simulation parameters ***** \n");
    fp_info.print( " L/di=%e\n", L_di);
    fp_info.print( " L/de=%e\n", L/de);
    fp_info.print( " Ti/Te=%e\n", Ti_Te );
    fp_info.print( " wpe/wce = %e\n", wpe_wce );
    fp_info.print( " mi/me =%e\n", mi_me );
    fp_info.print( " theta =%e\n", theta );
    fp_info.print( " taui =%e\n", taui );
    fp_info.print( " num_step = %i\n", num_step );
    fp_info.print( " Lx/de = %e\n", Lx/de );
    fp_info.print( " Ly/de = %e\n", Ly/de );
    fp_info.print( " Lz/de =%e\n", Lz/de );
    fp_info.print( " Lx/di = %e\n", Lx/di );
    fp_info.print( " Ly/di = %e\n", Ly/di );
    fp_info.print( " Lz/di =%e\n", Lz/di );
    fp_info.print( " nx = %e\n", nx );
    fp_info.print( " ny = %e\n", ny );
    fp_info.print( " nz =%e\n", nz );
    fp_info.print( " courant = %e\n", c*dt/dg );
    fp_info.print( " nproc  = %e\n", nproc() );
    fp_info.print( " nppc = %e\n", nppc );
    fp_info.print( " b0 =%e\n", b0 );
    fp_info.print( " di = %e\n", di );
    fp_info.print( " Ne = %e\n", Ne );
    fp_info.print( " total # of particles = %e\n", 2*Ne );
    fp_info.print( " dt*wpe = %e\n", wpe*dt );
    fp_info.print( " dt*wce = %e\n", wce*dt );
    fp_info.print( " dt*wci = %e\n", wci*dt );
    fp_info.print( " energies_interval: %i\n", energies_interval);
    fp_info.print( " dx/de =%e\n", Lx/(de*nx) );
    fp_info.print( " dy/de =%e\n", Ly/(de*ny) );
    fp_info.print( " dz/de =%e\n", Lz/(de*nz) );
    fp_info.print( " L/debye =%e\n", L/(vthe/wpe) );
    fp_info.print( " dx/rhoi =%e\n", (Lx/nx)/(vthi/wci) );
    fp_info.print( " dx/rhoe = %e\n", (Lx/nx)/(vthe/wce) );
    fp_info.print( " dx/debye = %e\n", (Lx/nx)/(vthe/wpe) );
    fp_info.print( " n0 =            %e\n", n0 );
    fp_info.print( " vthi/c =%e\n", vthi/c );
    fp_info.print( " vthe/c =%e\n", vthe/c );
    fp_info.print( " ***************************\n");
    fp_info.close();


    // for the parallized translate.f90 written by Vadim
    // write binary info file

    if ( ! (fp_info.open("info.bin", io_write)==ok) ) ERROR(("Cannot open file."));

        fp_info.write(&topology_x, 1 );
        fp_info.write(&topology_y, 1 );
        fp_info.write(&topology_z, 1 );

        fp_info.write(&Lx, 1 );
        fp_info.write(&Ly, 1 );
        fp_info.write(&Lz, 1 );

        fp_info.write(&nx, 1 );
        fp_info.write(&ny, 1 );
        fp_info.write(&nz, 1 );

        fp_info.write(&dt, 1 );

        fp_info.write(&mi_me, 1 );
        fp_info.write(&wpe_wce, 1 );
        fp_info.write(&vthe, 1 );
        fp_info.write(&vthi, 1 );

        fp_info.close();

}

  ////////////////////////////
  // Load fields


  // Define some function to load profiles

# define BX b0*tanh(z/L)
# define BY sqrt(b0*b0 + bg*bg*b0*b0 - BX*BX)
# define VDY -0.5*(b0/L)/(cosh(z/L)*cosh(z/L))
# define VDX VDY*BX/BY
# define VD sqrt(VDX*VDX+VDY*VDY)
# define GVD 1./sqrt(1.-VD*VD/(c*c))

// modified perturbation
# define DBX dbx*cos(2.0*pi*(x-0.5*Lx)/Lpert)*sin(pi*z/Lz)
# define DBZ dbz*cos(pi*z/Lz)*sin(2.0*pi*(x-0.5*Lx)/Lpert)


//LO

  sim_log( "Loading fields" );
  //set_region_field( everywhere, 0, 0, 0, BX+DBX, BY, DBZ);
  set_region_field( everywhere, 0, 0, 0, (BX+DBX)*cs+BY*sn, -(BX+DBX)*sn+BY*cs, DBZ);

  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specied as logical equations (i.e. x>0 && x+y<2)

  // LOAD PARTICLES

  sim_log( "Loading particles" );

  // Do a fast load of the particles

  //  seed_entropy( rank() );  //Generators desynchronized
  seed_entropy( 101 );  // Kevin said it should be this way
  double xmin = grid->x0 , xmax = grid->x0+(grid->dx)*(grid->nx);
  double ymin = grid->y0 , ymax = grid->y0+(grid->dy)*(grid->ny);
  double zmin = grid->z0 , zmax = grid->z0+(grid->dz)*(grid->nz);

  // Load Harris population

  sim_log( "-> Force Free Sheet" );

  repeat ( Ne/nproc() ) {
    double x, y, z, ux, uy, uz, upa1, upe1, uz1, gu1;

    x = uniform( rng(0), xmin, xmax);
    y = uniform( rng(0), ymin, ymax);
    z = uniform( rng(0), zmin, zmax);

   // inject_particles() will return an error for particles no on this
   // node and will not inject particle locally

   //  Load electrons as drifting Maxwellian with velocity specified to be consistent with B field

    upa1 = normal( rng(0), 0, vthe);
    upe1 = normal( rng(0), 0 ,vthe);
    uz1 = normal( rng(0), 0, vthe);
    gu1 = sqrt(1.0+upa1*upa1+upe1*upe1+uz1*uz1);
    ux = (GVD*upa1*VDX/VD - upe1*VDY/VD) + GVD*VDX*gu1;
    uy = (GVD*upa1*VDY/VD + upe1*VDX/VD) + GVD*VDY*gu1;
    uz = uz1;

    //L.O. inject_particle(electron, x, y, z, ux, uy, uz, qe, 0, 0 );
    inject_particle(electron, x, y, z, ux*cs+uy*sn,-ux*sn+uy*cs, uz, weight, 0, 0 );

    //  Ions are spatially uniform Maxwellian with no drifts

    upa1 = normal(rng(0),0,vthi);
    upe1 = normal(rng(0),0,vthi);
    uz1 = normal(rng(0),0,vthi);
    gu1 = sqrt(1.0+upa1*upa1+upe1*upe1+uz1*uz1);
    ux = (-GVD*upa1*VDX/VD + upe1*VDY/VD) - GVD*VDX*gu1;
    uy = (-GVD*upa1*VDY/VD - upe1*VDX/VD) - GVD*VDY*gu1;
    uz = uz1;


    //L.O. inject_particle(ion, x, y, z, ux, uy, uz, qi, 0, 0 );
    inject_particle(ion, x, y, z, ux*cs+uy*sn,-ux*sn+uy*cs, uz,weight, 0, 0 );

  }

  sim_log( "Finished loading particles" );

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

	global->fdParams.format = band;

	sim_log ( "Fields output format = band" );

	global->hedParams.format = band;

	sim_log ( "Electron species output format = band" );

	global->hHdParams.format = band;

	sim_log ( "Ion species output format = band" );

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

	// relative path to fields data from global header
	sprintf(global->fdParams.baseDir, "fields");

	// base file name for fields output
	sprintf(global->fdParams.baseFileName, "fields");

	global->fdParams.stride_x = 1;
	global->fdParams.stride_y = 1;
	global->fdParams.stride_z = 1;

	// add field parameters to list
	global->outputParams.push_back(&global->fdParams);

	sim_log ( "Fields x-stride " << global->fdParams.stride_x );
	sim_log ( "Fields y-stride " << global->fdParams.stride_y );
	sim_log ( "Fields z-stride " << global->fdParams.stride_z );

	// relative path to electron species data from global header
	sprintf(global->hedParams.baseDir, "hydro");

	// base file name for fields output
	sprintf(global->hedParams.baseFileName, "ehydro");

	global->hedParams.stride_x = 1;
	global->hedParams.stride_y = 1;
	global->hedParams.stride_z = 1;

	// add electron species parameters to list
	global->outputParams.push_back(&global->hedParams);

	sim_log ( "Electron species x-stride " << global->hedParams.stride_x );
	sim_log ( "Electron species y-stride " << global->hedParams.stride_y );
	sim_log ( "Electron species z-stride " << global->hedParams.stride_z );

	// relative path to electron species data from global header
	sprintf(global->hHdParams.baseDir, "hydro");

	// base file name for fields output
	sprintf(global->hHdParams.baseFileName, "Hhydro");

	global->hHdParams.stride_x = 1;
	global->hHdParams.stride_y = 1;
	global->hHdParams.stride_z = 1;

	sim_log ( "Ion species x-stride " << global->hHdParams.stride_x );
	sim_log ( "Ion species y-stride " << global->hHdParams.stride_y );
	sim_log ( "Ion species z-stride " << global->hHdParams.stride_z );

	// add electron species parameters to list
	global->outputParams.push_back(&global->hHdParams);

    /*--------------------------------------------------------------------------
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

	global->fdParams.output_variables( electric | magnetic );
	global->hedParams.output_variables( current_density | charge_density | stress_tensor );
	global->hHdParams.output_variables( current_density | charge_density | stress_tensor );

	//global->fdParams.output_variables( all );
	//global->hedParams.output_variables( all );
	//global->hHdParams.output_variables( all );

	/*--------------------------------------------------------------------------
	 * Convenience functions for simlog output
	 *------------------------------------------------------------------------*/

	char varlist[512];
	create_field_list(varlist, global->fdParams);

	sim_log ( "Fields variable list: " << varlist );

	create_hydro_list(varlist, global->hedParams);

	sim_log ( "Electron species variable list: " << varlist );

	create_hydro_list(varlist, global->hHdParams);

	sim_log ( "Ion species variable list: " << varlist );




	/* ---------------------------------------------

	   now add parameters for the energy diagnostics

	 --------------------------------------------- */

	global->ede.sp_id = electron->id;
	global->ede.vth = sqrt(2.0)*vthe;
	sprintf(global->ede.fname,global->hedParams.baseFileName);

	global->edi.sp_id = ion->id;
	global->edi.vth = sqrt(2.0)*vthi;
	sprintf(global->edi.fname, global->hHdParams.baseFileName);

	global->nex  = 6;
	global->emax = 120;

        std::vector<std::string> pythonNames;
        pythonNames.push_back("scatterplot.py");
        pythonNames.push_back("histogram.py");
        pythonNames.push_back("particlewriter.py");

        coprocessorinitialize(pythonNames);

        sim_log ( "Finished coprocessor initialization ");

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

} //begin_initialization

#define should_dump(x) \
  (global->x##_interval>0 && remainder(step(), global->x##_interval) == 0)

begin_diagnostics {


	/*--------------------------------------------------------------------------
	 * NOTE: YOU CANNOT DIRECTLY USE C FILE DESCRIPTORS OR SYSTEM CALLS ANYMORE
	 *
	 * To create a new directory, use:
	 *
	 *   dump_mkdir("full-path-to-directory/directoryname")
	 *
	 * To open a file, use: FileIO class
	 *
	 * Example for file creation and use:
	 *
	 *   // declare file and open for writing
	 *   // possible modes are: io_write, io_read, io_append,
	 *   // io_read_write, io_write_read, io_append_read
	 *   FileIO fileIO;
	 *   FileIOStatus status;
	 *   status= fileIO.open("full-path-to-file/filename", io_write);
	 *
	 *   // formatted ASCII  output
	 *   fileIO.print("format string", varg1, varg2, ...);
	 *
	 *   // binary output
	 *   // Write n elements from array data to file.
	 *   // T is the type, e.g., if T=double
	 *   // fileIO.write(double * data, size_t n);
	 *   // All basic types are supported.
	 *   fileIO.write(T * data, size_t n);
	 *
	 *   // close file
	 *   fileIO.close();
     *------------------------------------------------------------------------*/

     /*--------------------------------------------------------------------------
	 * Data output directories
	 * WARNING: The directory list passed to "global_header" must be
	 * consistent with the actual directories where fields and species are
	 * output using "field_dump" and "hydro_dump".
	 *
	 * DIRECTORY PATHES SHOULD BE RELATIVE TO
	 * THE LOCATION OF THE GLOBAL HEADER!!!
     *------------------------------------------------------------------------*/


  /*--------------------------------------------------------------------------
   * Normal rundata dump
   *------------------------------------------------------------------------*/
  if(step()==0) {
		dump_mkdir("fields");
		dump_mkdir("hydro");
		dump_mkdir("rundata");
		dump_mkdir("restart0");
		dump_mkdir("restart1");  // 1st backup
		dump_mkdir("restart2");  // 2nd backup
		dump_mkdir("particle");

		dump_grid("rundata/grid");
		dump_materials("rundata/materials");
		dump_species("rundata/species");
		global_header("global", global->outputParams);
	} // if

	/*--------------------------------------------------------------------------
	 * Normal rundata energies dump
	 *------------------------------------------------------------------------*/
	if(should_dump(energies)) {
	  dump_energies("rundata/energies", step() == 0 ? 0 : 1);
	} // if

	/*--------------------------------------------------------------------------
	 * Field data output
	 *------------------------------------------------------------------------*/

	if(step() == 1 || should_dump(fields)) field_dump(global->fdParams);

	/*--------------------------------------------------------------------------
	 * Electron species output
	 *------------------------------------------------------------------------*/

	if(should_dump(ehydro)) hydro_dump("electron", global->hedParams);

	/*--------------------------------------------------------------------------
	 * Ion species output
	 *------------------------------------------------------------------------*/

	if(should_dump(Hhydro)) hydro_dump("ion", global->hHdParams);

	/*--------------------------------------------------------------------------
	 * Energy Spectrum Output
	 *------------------------------------------------------------------------*/

	//	#include "energy.cxx"   //  Subroutine to compute energy spectrum diagnostic

         //Vadim:
         //#include "dissipation.cxx"
         //#include "Ohms_exp_all_v2.cxx"

	/*--------------------------------------------------------------------------
	 * Restart dump
	 *------------------------------------------------------------------------*/

  //Vadim:
	if (step() && !(step()%global->restart_interval))
	  global->write_restart = 1; // set restart flag. the actual restart files are written during the next step
	else
	  if (global->write_restart) {

	    global->write_restart = 0; // reset restart flag

	    double dumpstart = uptime();

	    if(!global->rtoggle) {
	      global->rtoggle = 1;
	      checkpt("restart1/restart", 0);
	    }
	    else {
	      global->rtoggle = 0;
	      checkpt("restart2/restart", 0);
	    } // if

	    mp_barrier(  ); // Just to be safe

	    double dumpelapsed = uptime() - dumpstart;
	    sim_log("Restart duration "<<dumpelapsed);

      //Vadim
    if (rank()==0) {

        FileIO fp_restart_info;
        if ( ! (fp_restart_info.open("latest_restart", io_write)==ok) ) ERROR(("Cannot open file."));
        if(!global->rtoggle) {
          fp_restart_info.print("restart restart2/restart");
        } else
          fp_restart_info.print("restart restart1/restart");

        fp_restart_info.close();
       }

	} // if


  // Dump particle data

	char subdir[36];
	if ( should_dump(eparticle) && step() !=0 && step() > 0*(global->fields_interval)  ) {
	  sprintf(subdir,"particle/T.%d",step());
	  dump_mkdir(subdir);
	  sprintf(subdir,"particle/T.%d/eparticle",step());
	  dump_particles("electron",subdir);
	}

        //float mytime = step*(grid->dt);
        long long mystep = step();
        double mytime = grid->dt*step();
        int topology[3] = {static_cast<int>(global->topology_x),
                           static_cast<int>(global->topology_y),
                           static_cast<int>(global->topology_z)};
        coprocessorProcess(mystep, mytime, this, topology, global->outputParams);


  // Shut down simulation when wall clock time exceeds global->quota_sec.
  // Note that the mp_elapsed() is guaranteed to return the same value for all
  // processors (i.e., elapsed time on proc #0), and therefore the abort will
  // be synchronized across processors. Note that this is only checked every
  // few timesteps to eliminate the expensive mp_elapsed call from every
  // timestep. mp_elapsed has an ALL_REDUCE in it!

        //Vadim:
	if  (( step()>0 && global->quota_check_interval>0 && (step()&global->quota_check_interval)==0) || (global->write_end_restart) ) {
	   if ( (global->write_end_restart) ) {

		   global->write_end_restart = 0; // reset restart flag

             sim_log( "Allowed runtime exceeded for this job.  Terminating....\n");
	     double dumpstart = uptime();

	    checkpt("restart0/restart",0);

	    mp_barrier(  ); // Just to be safe
	    sim_log( "Restart dump restart completed." );
	    double dumpelapsed = uptime() - dumpstart;
	    sim_log("Restart duration "<<dumpelapsed);

            //Vadim:
            if (rank()==0) {
               FileIO fp_restart_info;
               if ( ! (fp_restart_info.open("latest_restart", io_write)==ok) ) ERROR(("Cannot open file."));
               fp_restart_info.print("restart restart0/restart");
               fp_restart_info.close();
             }

            coprocessorfinalize();

	    exit(0); // Exit or abort?
	  }
	   if( uptime( ) > global->quota_sec )   global->write_end_restart = 1;
	}

} // end diagnostics

// *******************  PARTICLE INJECTION  - OPEN BOUNDARY ***************************

begin_particle_injection {


}  // end particle injection


//   *******************  CURRENT INJECTION ***************************

begin_current_injection {

  // No current injection for this simulation

}

//   *******************  FIELD INJECTION ***************************

begin_field_injection {


}  // end field injection


begin_particle_collisions {

  // No particle collisions in this simulation


}
