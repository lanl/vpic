//========================================================================
//
// Interface dynamics for Brian's LDRD DR
//
// rm -rf rundata field ehydro Hhydro Hehydro particle fft movie poynting velocity restart
//
//========================================================================

// ASCII Logging of IO writes to disk.  Since wall clock elapsed time is written,
// and these require an MPI_Allreduce, don't include this macro inside code
// which only executes on a single processor, else deadlock will ensue.
// Turn off logging by setting "#if 0"

// Flag to put barriers in the begin_diagnostics{} stub in order to help
// debug I/0
#define DEBUG_SYNCHRONIZE_IO 0

// Employ turnstiles to partially serialize the high-volume file writes.
// In this case, the restart dumps.  Set NUM_TURNSTILES to be the desired
// number of simultaneous writes.
#define NUM_TURNSTILES 128

  // From grid/partition.c: used to determine which domains are on edge
# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {                   \
    int _ix, _iy, _iz;                                                    \
    _ix  = (rank);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(global->topology_x);   /* iy = iy+gpy*iz */            \
    _ix -= _iy*int(global->topology_x);   /* ix = ix */                   \
    _iz  = _iy/int(global->topology_y);   /* iz = iz */                   \
    _iy -= _iz*int(global->topology_y);   /* iy = iy */                   \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \
  } END_PRIMITIVE

begin_globals {
  float vthe; 			// vthe/c   <- these are needed to make movie files
  float vthi_I2;                // vthi_I2/c
  float vthi_I1;                // vthi_I1/c
  int field_interval;           // how frequently to dump field built-in diagnostic
  int energies_interval;
  int restart_interval; 	// how frequently to write restart file.
  int quota_check_interval;    // how often to check if quote exceeded
  int velocity_interval;        // how frequently to dump velocity space
  int fft_ex_interval;          // how frequently to save ex fft data
  int eparticle_interval;       // how frequently to dump particle data
  int I1particle_interval;       //
  int I2particle_interval;      //
  int rtoggle;                  // Enables save of last two restart dumps for safety
  double quota_sec;             // Run quota in seconds

  int do_collisions;            // Flag for whether to do collisions
  int self_collisions_only;     // Whether to turn on cross-species collisions
  int tstep_coll;               // How many time steps to sub-cycle the collision operator
  int load_particles;           // Whether to load particles

  // For collision operator

  int nppc_max;                 // Maximum number particles/cell possible
  double inv_nppc;              // 1 / nppc in the high density region - used in density scaling
  double cvar;                  // Collision variable
  double mime_I1;               // proton to electron mass ratio
  double mime_I2;               // alpha to electron mass ratio
  double Z_e;                   // charge (in electronic charge e) of electron
  double Z_I1;                  // charge (in electronic charge e) of ion species 1
  double Z_I2;                  // charge (in electronic charge e) of ion species 2

  double topology_x;            // domain topology needed to normalize Poynting diagnostic
  double topology_y;
  double topology_z;

  double wpe1ps;

  // Dump parameters for standard VPIC output formats
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hI1dParams;
  DumpParameters hI2dParams;
  std::vector<DumpParameters *> outputParams;
};

begin_initialization {

  sim_log( "*** Begin initialization. ***" );
  barrier(); // Barrier to ensure we started okay.
  sim_log( "*** Begin initialization2. ***" );

  float elementary_charge  = 4.8032e-10;       // stat coulomb
  float elementary_charge2 = elementary_charge * elementary_charge;
  float speed_of_light     = 2.99792458e10;    // cm/sec
  float m_e                = 9.1094e-28;  // g

  float mp_me              = 1836; // Proton to electron mass ratio - physical value is 1836.
                                   // All ions' masses are scaled to this ratio
  float k_boltz            = 1.6022e-12;       // ergs/eV
  float mec2               = m_e*speed_of_light*speed_of_light/k_boltz;
  float mpc2               = mec2 * mp_me;
  float eps0               = 1;

  double cfl_req           = 0.98;      // How close to Courant should we try to run
  double damp              = 0.0;         // Level of radiation damping
  double iv_thick          = 2;         // Thickness (in cells) of imperm. vacuum region

  float t_e                = 4000;               // Electron temp, eV
  float t_i                = t_e;              // Ion temp, eV
//???????????????????????????????????????????????????????????????????????
  float box_size_x         = 8.0 * 1e-4;  // in cm; 8 um
  float box_size_z         = 0.1     * 1e-4;  // in cm; 0.003 um

  int load_particles       = 1;         // Flag to turn off particle load for testing wave launch.
//??????????????????????????????????????????????????????????????????
  double nppc              = 4000; // for species I1 if qi_I1 > qi_I2; for species I2 if qi_I2 > qi_I1
                                  // electron and ion number in each region are the same

  // FIXME:  Put in the real values here rather than approximate values:
  float A_I1               = 27;       // Al
  float Z_I1               = 13;
  float A_I2               = 2;       // D
  float Z_I2               = 1;
  float mic2_I1            = mpc2*A_I1;
  float mic2_I2            = mpc2*A_I2;
  float mime_I1            = mic2_I1/mec2;
  float mime_I2            = mic2_I2/mec2;

  double uthe              = sqrt(t_e/mec2);     // vthe/c
  double uthi_I1           = sqrt(t_i/mic2_I1);  // vthi/c
  double uthi_I2           = sqrt(t_i/mic2_I2);  // vthi/c

  // NOTE: The following definitions make the assumption that Te = Ti = constant everywhere at t=0 and that
  //       total pressure balance is established across the interface at t=0
  double n_i_I2            = 2.9197e24;                              // ion density of D (1/cc)
  double n_i_I1            = n_i_I2 * (Z_I2 + 1.0) / (Z_I1 + 1.0) ;  // ion density of Al (assuming pressure balance)

  double rho0              = n_i_I1 * A_I1 * m_e * mp_me;            // ion mass density of Al in g/cc

  double n_e;              // Set n_e to be the electron density in region of highest e- density
  if ( n_i_I1 * Z_I1 > n_i_I2 * Z_I2 ) {  // region 1 has highest e- density

    n_e                    = Z_I1 * n_i_I1;                          // e- density in Al ( > e- density in D)

  } else {                                // region 2 has highest e- density

    sim_log("REGION 2 HAS HIGHER ELECTRON DENSITY - SHOULD NOT HAPPEN!!!" );
    exit(0);
  }

  // delta is skin depth in initial gold plasma
  double delta             = speed_of_light / sqrt( 4.0 * M_PI * n_e * elementary_charge2 / m_e ); // cm
  double debye             = uthe*delta;
  double wpe1ps            = 1e-12 * speed_of_light/delta;

//????????????????????????????????????????????????????????????????????????????????????
  double nx                = 80000;
  double ny                = 1;         // 2D problem in x-z plane
  double nz                = 1000;

  double hx                = box_size_x/(delta*nx);   // in c/wpe
  double hz                = box_size_z/(delta*nz);
  double hy                = hz;

  double cell_size_x       = delta*hx/debye;         // Cell size in Debye lengths
  double cell_size_z       = delta*hz/debye;         // Cell size in Debye lengths

  double Lx                = nx*hx;          // in c/wpe
  double Ly                = ny*hy;
  double Lz                = nz*hz;

  // Size of the perturbation in same system of units as Lx.
//????????????????????????????????????????????????????????????????????????????????????
  double perturb_ampl      = Lx * 0.125;

  double dt = cfl_req*courant_length(Lx, Ly, Lz, nx, ny, nz);

// BJA - change topology here to go with a 2D decomposition

//????????????????????????????????????????????????????????????????????????????????????
  double topology_x        = 2500;
  double topology_y        = 1;
  double topology_z        = 20;

  // DEBUG - smaller problem of 16 cores
  topology_x = 4;
  topology_z = 4;
  nx        /= 625;
  nz        /= 5;
  Lx        /= 625;
  Lz        /= 5;

  // work out stopping time in terms of ion sound speed crossing time across box in x for Au

  double cs_crossing_time_ps = 1e12 * box_size_x / sqrt( Z_I2 * t_e * k_boltz / (m_e * mp_me * A_I2) ) ; // ps
  float t_stop               = cs_crossing_time_ps * 100 * wpe1ps;  // Runtime (in 1/wpe) equal to 100 Cs crossing times

  int fft_ex_interval        = int(M_PI / (1.2 * dt));       // Num steps between writing Ex in fft slice
//????????????????????????????????????????????????????????????????????????????????????
  int energies_interval = 1000;
//int field_interval         = int( cs_crossing_time_ps * 0.001 * wpe1ps / dt );         // Num. steps between saving field, hydro data
  int field_interval         = 5000;         // Num. steps between saving field, hydro data
  int restart_interval       = 4000;  // DEBUG - test launch

  int quota_check_interval   = 20;
  int velocity_interval      = field_interval;     //  Num steps between writing poynting flux; not used in NIC

  int eparticle_interval     = 0;
  int I1particle_interval    = 0;
  int I2particle_interval    = 0;

//????????????????????????????????????????????????????????????????????????????????????
  int ele_sort_freq          = 160;
  int ion_sort_freq          = ele_sort_freq;

//????????????????????????????????????????????????????????????????????????????????????
  double quota              = 15.7;            // Run quota in hours.
  double quota_sec           = quota*3600;  // Run quota in seconds.

  double Ne                  = nppc*nx*ny*nz;                // Number of macro electrons in box in region with highest e- density
  Ne                         = trunc_granular(Ne, nproc());  // Make Ne divisible by number of processors
  double Npe                 = Lx*Ly*Lz;                     // Number of physical electrons in box in sim. units: density is 1, wpe = 1

  double qe_1, qe_2, qi_I1, qi_I2;

  // NOTE: The following definitions make the assumption that Te = Ti = constant everywhere at t=0

  // NOTE: RECOGNIZE THAT THE SEMANTICS ARE DIFFERENT FOR qi_I1 VS. qi_I2 DEPENDING
  //       ON WHICH REGION HAS THE HIGHEST ELECTRON DENSITY! THE LOWER DENSITY REGION
  //       qi_I1 OR qi_I2 ONLY IS USED FOR DEFINING THE AVERAGE NPPC IN THE PARTICLE
  //       LOAD; THIS IS DONE BECAUSE WE ARE FORCING UNIFORM WEIGHTING OF PARTICLES IN
  //       THE COLLISION OPERATOR.

  if ( n_i_I1 * Z_I1 > n_i_I2 * Z_I2 ) {  // region 1 has highest e- density

    qe_1                = -Npe/( Ne * Z_I1 * ( (Z_I2+1.0)/(Z_I1+1.0) ) );  // Charge per macro electron in region 1 - note: we load Z_I1 electrons per load.
    qi_I1               = -qe_1 * Z_I1;                                    // Charge per I1 macro ion =  Z_I1 * charge per macro electron
    qe_2                = qe_1;                                            // Charge per macro electron in region 2 - same as region 1
    qi_I2               = -qe_1 * Z_I2;                                    // Charge per I2 macro ion =  Z_I2 * charge per macro electron

  } else {                                // region 2 has highest e- density

    sim_log("REGION 2 HAS HIGHER ELECTRON DENSITY - SHOULD NOT HAPPEN!!!" );
    exit(0);

#   if 0
    qe_2                = -Npe/Ne;                      // Charge per macro electron in region 2
    qi_I2               = -qe_2;                        // Charge per I2 macro ion
    qi_I1               = -qe_2 * (n_i_I1 * Z_I1) / (n_i_I2 * Z_I2) ; //  (qi_I1 / qi_I2) * nppc defines the average number of ion macroparticles loaded in region 1
#   endif
  }

#if 0
  // old
  double qe_2                = -Npe/Ne;                      // Charge per macro electron in region 2
  double qi_I2               = -qe_2;                        // Charge per I2 macro ion
  double qe_1                = qe_2 * (n_i_I1 * Z_I1) / (n_i_I2 * Z_I2) ; // Charge per macro electron in region 1
  double qi_I1               = -qe_1;                        // Charge per H macro ion
#endif

  // Collision parameters
  // In CGS, variance of tan theta = 2 pi e^4 n_e dt_coll loglambda / (m_ab^2 c^3)
  sim_log("Setting up particle collision parameters. ");
  double nu_e;                           // Electron collision frequency nu_0^{e\e} in units of wpe
                                         // (See NRL Formulary, p. 32)
  // Define log-lambda (we could compute this, but let's set it to fixed value for now)
  double loglambda = 3;
  {
    // Compute nu_e from plasma parameters given.  Spitzer collisions assumed.
    // According to NRL plasma formuary, nu_e = nu_0^{e\e} (2 u_the^3)
    // where nu_0^{e\e} = 4 pi e^4 n_e lambda / (m_e^2 v_e^3)
    // We need to compute n_0^{e\e} for v_e = u_the*c

    // First, compute n_e_cgs
    double n_e_cgs = n_e;


    // Then, compute nu_e, assuming that v_e = uthe*c
    double nu_e_cgs = 4*M_PI*pow( 4.8032e-10,4 )*n_e_cgs*loglambda*pow( 9.1094e-28,-2 )*pow( uthe*3e10,-3 );

    // Finally, normalize to wpe
    double wpe_cgs = 5.64e4*sqrt(n_e_cgs);
    nu_e = nu_e_cgs / wpe_cgs;
    sim_log("* nu_e/wpe= "<<nu_e);
  }

  int do_collisions = 1;                 // Flag to turn on/off Takizuka and Abe collisions
  int self_collisions_only = 0;          // Flag for particles to scatter off same species only
  int tstep_coll = ion_sort_freq;        // How frequently to apply collision operator:
                                         // N.B. Needs to be a multiple of ion_sort_freq
  double dt_coll = tstep_coll*dt;        // Units of wpe

  // cvar = variance of tan theta in collision operator for electron-electron collisions
  //      = (8 pi e^4 n_e log-lambda) / (m_e^2 c^3)
  double cvar=2.0*nu_e*pow(uthe,3)*dt_coll;  // Coefficient for collision model
  double nppc_max = 40*nppc;             // Used internally in collision algorithm.  Should be sufficient?

  // Throw warning if collision interval is too large
  if ( nu_e*2*dt_coll>0.1 ) {
    sim_log( "*** Warning: scattering interval may be too large to sample electron scattering accurately." );
    sim_log( "*** typical tan theta step: "<<sqrt(nu_e*2*dt_coll) );
  }


  // Initialize with Kim's similarity solution
  double delta_x0, half_box_size_in_xi;
//delta_x0                = 0.5 * 1.0e-4;                          /* initial layer thickness in cm */
  delta_x0                = 0.05 * 1.0e-4;                          /* initial layer thickness in cm */
  half_box_size_in_xi     = (0.5 * box_size_x) / delta_x0;

  // Proton ion self diffusion coefficient for pure hydrogen plasma at mass density 1 g/cc at
  // ion temperature 1 keV (used in initializing diffusion velocities in layer)
  double D_proton         = 628.5 * ( 4.0 / loglambda );           // g/(cm s)

  // PRINT SIMULATION PARAMETERS
  sim_log("***** Simulation parameters *****");
  sim_log("* Processors:                    "<<nproc());
  sim_log("* nu_e*2*dt_coll=                "<<nu_e*2*dt_coll);
  sim_log("* Time step, max time, nsteps =  "<<dt<<" "<<t_stop<<" "<<int(t_stop/(dt)));
  sim_log("* wpe1ps =                       "<<wpe1ps);
  sim_log("* Debye length, delta, cell size in x, z = "<<debye<<" "<<delta<<" "<<cell_size_x<<" "<<cell_size_z);
  sim_log("* Lx, Ly, Lz =                   "<<Lx<<" "<<Ly<<" "<<Lz);
  sim_log("* perturb_ampl =                 "<<perturb_ampl);
  sim_log("* nx, ny, nz =                   "<<nx<<" "<<ny<<" "<<nz);
  sim_log("* Charge/macro electron 1      = "<<qe_1);
  sim_log("* Charge/macro electron 2      = "<<qe_2);
  sim_log("* Charge/macro ion spec 1      = "<<qi_I1);
  sim_log("* Charge/macro ion spec 2      = "<<qi_I2);
  sim_log("* Average particles/processor:   "<<Ne/nproc());
  sim_log("* Average particles/cell:        "<<nppc);
  sim_log("* Omega_pe:                      "<<1);
  sim_log("* Plasma density                 "<<n_e);
  sim_log("* T_e, T_i, m_e, m_i_I1, m_i_I2:  "<<t_e<<" "<<t_i<<" "<<1<<" "<<mime_I1<<" "<<mime_I2);
  sim_log("* Radiation damping:             "<<damp);
  sim_log("* Fraction of courant limit:     "<<cfl_req);
  sim_log("* vthe/c:                        "<<uthe);
  sim_log("* vthi_I1/c, vth_I2/c:           "<<uthi_I1<<" "<<uthi_I2);
  sim_log("* restart interval:              "<<restart_interval);
  sim_log("* energies_interval:             "<< energies_interval );
  sim_log("* quota_check_interval:          "<<quota_check_interval);
  sim_log("* velocity interval:             "<<velocity_interval);
  sim_log("* ex save interval:              "<<fft_ex_interval);
  sim_log("* quota (hours):                 "<<quota);
  sim_log("* load_particles:                "<<(load_particles ? "Yes" : "No"));
  sim_log("* do_collisions:                 "<<(do_collisions ? "Yes" : "No"));
  sim_log("* self_collisions_only:          "<<(self_collisions_only ? "Yes" : "No"));
  sim_log("* tstep_coll:                    "<<tstep_coll);
  sim_log("* nppc_max:                      "<<nppc_max);
  sim_log("* nu_e:                          "<<nu_e);
  sim_log("* cvar:                          "<<cvar);
  sim_log("* mime_I1:                       "<<mime_I1);
  sim_log("* mime_I2:                       "<<mime_I2);
  sim_log("* ele_sort_freq:                 "<<ele_sort_freq);
  sim_log("* ion_sort_freq:                 "<<ion_sort_freq);
  sim_log("* delta_x0:                      "<<delta_x0 * 1e4<<" microns");
  sim_log("* half_box_size_in_xi:           "<<half_box_size_in_xi);
  sim_log("* D_proton:                      "<<D_proton<<" g/(cm s)");
  sim_log("* rho0:                          "<<rho0<<" g/cc");
  sim_log("*********************************");


  // SETUP HIGH-LEVEL SIMULATION PARMETERS
  sim_log("Setting up high-level simulation parameters. ");
  num_step             = int(t_stop/(dt));
  status_interval      = 200;
//??????????????????????????????????????????????????????????????????????????????????????????
  sync_shared_interval = status_interval/1;
  clean_div_e_interval = status_interval/1;
  clean_div_b_interval = status_interval/10;

  // For maxwellian reinjection, we need more than the default number of
  // passes (3) through the boundary handler
  num_comm_round = 6;

  global->field_interval           = field_interval;
  global->restart_interval         = restart_interval;
  global->quota_check_interval = quota_check_interval;
  global->energies_interval  = energies_interval;
  global->velocity_interval        = velocity_interval;
  global->fft_ex_interval     = fft_ex_interval;
  global->vthe                     = uthe;     // c=1
  global->vthi_I2                  = uthi_I2;  // c=1
  global->vthi_I1                   = uthi_I1;   // c=1
  global->wpe1ps                   = wpe1ps;
  global->quota_sec                = quota_sec;
  global->rtoggle                  = 0;
  global->eparticle_interval       = eparticle_interval;
  global->I1particle_interval      = I1particle_interval;
  global->I2particle_interval      = I2particle_interval;
  global->load_particles           = load_particles;
  global->do_collisions            = do_collisions;
  global->self_collisions_only     = self_collisions_only;
  global->tstep_coll               = tstep_coll;
  global->nppc_max                 = (int)nppc_max;
  global->cvar                     = cvar;
  global->mime_I1                   = mime_I1;
  global->mime_I2                  = mime_I2;

  // BJA - add these to global array to eliminate magic numbers in collision operator
  global->Z_e                      = -1;
  global->Z_I1                     = Z_I1;
  global->Z_I2                     = Z_I2;

  // BJA - for collision operator density scaling
  global->inv_nppc                 = 1.0 / (double)(nppc);

  global->topology_x               = topology_x;
  global->topology_y               = topology_y;
  global->topology_z               = topology_z;

  // SETUP THE GRID
  sim_log("Setting up computational grid.");
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = 1;
  grid->eps0 = eps0;

  // Partition a periodic box among the processors sliced uniformly in x:
  define_periodic_grid( 0,         -0.5*Ly,    -0.5*Lz,       // Low corner
                        Lx,         0.5*Ly,     0.5*Lz,       // High corner
                        nx,         ny,         nz,           // Resolution
                        topology_x, topology_y, topology_z ); // Topology

  // Set up Maxwellian reinjection B.C.

  // SETUP THE SPECIES
  sim_log("Setting up species. ");

  //????????????????????????????????????????????????????????????????????????????

  // Ion particle array sizes can be made smaller than those for electrons - this definition
  // below uses the semantics above, that the ions in region II have the highest charge.

  double max_local_np_e            = 2.0*Ne/nproc();
  double max_local_np_i            = max_local_np_e / Z_I2;
  double max_local_nm_e            = max_local_np_e / 10.0;
  double max_local_nm_i            = max_local_nm_e / Z_I2;
  sim_log( "num electron, ion macroparticles: "<<max_local_np_e<<" "<<max_local_np_i );
  species_t *electron=NULL, *ion_I1=NULL, *ion_I2=NULL;
  sim_log("- Creating electron species.");
  electron = define_species("electron", -1, 1, max_local_np_e,
                            max_local_nm_e, ele_sort_freq, 1);
  sim_log("- Creating I1 ion species.");
  ion_I1   = define_species("I1", Z_I1, mime_I1, max_local_np_i,
                            max_local_nm_i, ion_sort_freq, 1);
  sim_log("- Creating I2 ion species.");
  ion_I2   = define_species("I2", Z_I2, mime_I2, max_local_np_i,
                            max_local_nm_i, ion_sort_freq, 1);
  sim_log("Done setting up species.");

  sim_log("Setting up Maxwellian reinjection boundary condition.");

  particle_bc_t * maxwellian_reinjection =
    define_particle_bc( maxwellian_reflux( species_list, entropy ) );
  set_reflux_temp( maxwellian_reinjection, electron, uthe,    uthe    );
  set_reflux_temp( maxwellian_reinjection, ion_I1,   uthi_I1, uthi_I1 );
  set_reflux_temp( maxwellian_reinjection, ion_I2,   uthi_I2, uthi_I2 );
  sim_log( "Done setting up Maxwellian reinjection boundary condition." );


  // Override x boundaries:  Field Boundaries    = absorbing
  //                         Particle Boundaries = "maxwellian reinjection"


  int ix, iy, iz;        // Domain location in mesh
  RANK_TO_INDEX( int(rank()), ix, iy, iz );

  if ( ix == 0 ) {                                 // Leftmost proc.
    set_domain_field_bc( BOUNDARY(-1,0,0), absorb_fields );
    //set_domain_particle_bc( BOUNDARY(-1,0,0), maxwellian_reinjection );
  }
  if ( ix == topology_x-1 ) {                        // Rightmost proc.
    set_domain_field_bc( BOUNDARY( 1,0,0), absorb_fields );
    //set_domain_particle_bc( BOUNDARY( 1,0,0), maxwellian_reinjection );
  }

  // SETUP THE MATERIALS

  sim_log("Setting up materials. ");
  define_material( "vacuum", 1 );
  define_field_array( NULL, damp );

  // Paint the simulation volume with materials and boundary conditions
# define iv_region ( x<hx*iv_thick || x>Lx -hx*iv_thick )  /* all boundaries are i.v. */

  set_region_bc( iv_region, maxwellian_reinjection, maxwellian_reinjection, maxwellian_reinjection );



//????????????????????????????????????????????????###############################################
// tanh load

#if 1

  //
  // Load particles
  //
  //

  // First order Juttner distribution correction applied to rejection method
  // for electron load: theta = T_e[keV] / 511.
  // rejection method: reject particle if mu < u^4 * one_over_eight_theta
  // where mu is sampled from uniform dist (0,1) and u = p/mc

  double one_over_eight_theta = ( 511. / ( 8. * t_e ) );

  // Multiplier on thermal velocity to account for weakly relativistic distribution
  // (necessary because Kevin's load makes < p^2/2m > the same for each particle, not
  // < K.E. >, which is somewhat lower due to relativistic kinematics. This is the
  // multiplier one applies to uthe in the sampler to force < K.E. > to be the same
  //
  // Note:  This is EMPIRICALLY obtained for the Molvig PRL problem and needs to be
  //        recomputed for each problem.
  //
  // FIXME: Need to derive this

//???????????????????????????????????????????????????????????????????????
//double uthe_relativistic_correction = 1.0048;
  double uthe_relativistic_correction = 1.0088;

  if ( load_particles ) {
    sim_log("Loading particles.");
    // Fast load of particles--don't bother fixing artificial domain correlations
    double xmin=grid->x0, xmax=grid->x1;
    double ymin=grid->y0, ymax=grid->y1;
    double zmin=grid->z0, zmax=grid->z1;

    double xi; // skin depths
    double chi = 20; // 1/microns

    // nL1, nR1, nL2, nR2 are all limiting number densities relative to 1.0
    // Note: with the tanh load, we are not initializing in global pressure balance
    double nL1, nR1, nL2, nR2;
//  nL1 = 0;
//  nR1 = n_i_I1 / n_e;
//  nL2 = n_i_I2 / n_e;
//  nR2 = 0;

    // Force the limits to match what we used in the step function profile

    nL1 = 0;
    nR1 =  (Z_I2 + 1.0) / (Z_I1 + 1.0);

    nL2 = 1.0;
    nR2 = 0;

    // Precompute some values for the tanh load
    double an1, bn1, an2, bn2;
    an1 = 0.5 * (nL1 + nR1) - nL1;
    bn1 = 0.5 * (nL1 + nR1);
    an2 = 0.5 * (nL2 + nR2) - nL2;
    bn2 = 0.5 * (nL2 + nR2);

    // Note: chi is in units of (1/microns). simulation spatial values are in units of skin depths.
    //       delta = skin depth in cm


    // Needed for Juttner distribution corrections.

    double vx, vy, vz, u2;


    // Sanity check on semantics
    if ( fabs(qi_I1) > fabs(qi_I2) ) {

      sim_log( "nparticles = "<< Ne/(topology_x*topology_y*topology_z) );

      repeat( Ne/(topology_x*topology_y*topology_z) ) {
        double x = uniform( rng(0), xmin, xmax );
        double y = uniform( rng(0), ymin, ymax );
        double z = uniform( rng(0), zmin, zmax );

        if ( iv_region ) continue;           // Particle fell in iv_region.  Don't load.

        // phase of sinusoidal perturbation
        double zphase = 0;
        if ( grid->nz > 1 ) zphase = 2.0 * M_PI * ( z / Lz );

        // xi represents position from center of box in microns (chi is in units of 1/microns)

        xi =  ( x - (Lx/2.0) - perturb_ampl*cos(zphase) );      // distance from center of box in skin depths
        xi *= delta * 1e4;                                      // distance from center of box in microns

        // Rejection method - use probability of loading an ion of species 1 and Z_I1 electrons at (x,y,z)

        if ( uniform( rng(0), 0, 1 ) < an1*tanh(chi*xi) + bn1 ) {
          inject_particle( ion_I1, x, y, z,
                           normal( rng(0), 0, uthi_I1 ),
                           normal( rng(0), 0, uthi_I1 ),
                           normal( rng(0), 0, uthi_I1 ), fabs(qi_I1), 0, 0 );

          repeat( Z_I1 ) {

            // Leading order correction to Juttner distribution; also fixes uthe
            // to account for weakly relativistic effects

            do {
              vx = normal( rng(0), 0, uthe * uthe_relativistic_correction );
              vy = normal( rng(0), 0, uthe * uthe_relativistic_correction );
              vz = normal( rng(0), 0, uthe * uthe_relativistic_correction );
              u2 = vx*vx + vy*vy + vz*vz;
            } while ( uniform( rng(0), 0, 1 ) < u2*u2*one_over_eight_theta );
            inject_particle( electron, x, y, z, vx, vy, vz, fabs(qe_1), 0, 0 );

          } // repeat
        } // if

        // Rejection method - use probability of loading an ion of species 2 and Z_I2 electrons at (x,y,z)

        if ( uniform( rng(0), 0, 1 ) < an2*tanh(chi*xi) + bn2 ) {
          inject_particle( ion_I2, x, y, z,
                           normal( rng(0), 0, uthi_I2),
                           normal( rng(0), 0, uthi_I2),
                           normal( rng(0), 0, uthi_I2), fabs(qi_I2), 0, 0 );

          repeat( Z_I2 ) {

            // Leading order correction to Juttner distribution; also fixes uthe
            // to account for weakly relativistic effects

            do {
              vx = normal( rng(0), 0, uthe * uthe_relativistic_correction );
              vy = normal( rng(0), 0, uthe * uthe_relativistic_correction );
              vz = normal( rng(0), 0, uthe * uthe_relativistic_correction );
              u2 = vx*vx + vy*vy + vz*vz;
            } while ( uniform( rng(0), 0, 1 ) < u2*u2*one_over_eight_theta );
            inject_particle( electron, x, y, z, vx, vy, vz, fabs(qe_2), 0, 0 );

          } // repeat
        } // if
      } // repeat

    } else { // qi_I1 <= qi_I2 -- not allowed!

      sim_log("ERROR: PARTICLE LOAD SEMANTICS REQUIRE qi_I1 > qi_I2");
      exit(0);

    } // if qi_I1 > qi_I2

  } // if load particles
#endif



// Step function load

#if 0
  // Load particles
  if ( load_particles ) {
    sim_log("Loading particles.");
    // Fast load of particles--don't bother fixing artificial domain correlations
    double xmin=grid->x0, xmax=grid->x1;
    double ymin=grid->y0, ymax=grid->y1;
    double zmin=grid->z0, zmax=grid->z1;

    // Kim's load similarity variable
    double xi;
    double I1_load_probability = (Z_I2 + 1.0) / (Z_I1 + 1.0);

    // Modify loading so that we use uniformly weighted particles throughout.
    // Conditional test on qi_I1 vs. qi_I2 is done so we can avoid putting
    // a conditional in the inner loop.

    if ( fabs(qi_I1) > fabs(qi_I2) ) {

      sim_log( "nparticles = "<< Ne/(topology_x*topology_y*topology_z) );

      repeat( Ne/(topology_x*topology_y*topology_z) ) {
        double x = uniform( rng(0), xmin, xmax );
        double y = uniform( rng(0), ymin, ymax );
        double z = uniform( rng(0), zmin, zmax );

        if ( iv_region ) continue;           // Particle fell in iv_region.  Don't load.

        // Similarity variable
        xi =  ( ((Lx - x) - (Lx/2.0)) / Lx ); // position (note: flipped left-right about box center in x) in units of fractions of half box size (dimensionless)
        xi *= half_box_size_in_xi;            // position (in units of similarity variable xi)

        if ( xi < 0 ) {                               // region 1 - left of layer in Kim's paper
          if ( uniform_rand( 0, 1 ) < I1_load_probability ) {
            inject_particle( ion_I1, x, y, z,
                             normal( rng(0), 0, uthi_I1 ),
                             normal( rng(0), 0, uthi_I1 ),
                             normal( rng(0), 0, uthi_I1 ), qi_I1, 0, 0 );

            repeat( Z_I1 ) {

              // Leading order correction to Juttner distribution; also fixes uthe
              // to account for weakly relativistic effects

              do {
                vx = normal( rng(0), 0, uthe*uthe_relativistic_correction );
                vy = normal( rng(0), 0, uthe*uthe_relativistic_correction );
                vz = normal( rng(0), 0, uthe*uthe_relativistic_correction );
                u2 = vx*vx + vy*vy + vz*vz;
              } while ( uniform( rng(0), 0, 1 ) < u2*u2*one_over_eight_theta );
              inject_particle( electron, x, y, z, vx, vy, vz, -qe_1, 0, 0 );

            }
          }
        } else { // ( xi > 0 )                        // region 2 - right of layer in Kim's paper
          inject_particle( ion_I2, x, y, z,
                           normal( rng(0), 0, uthi_I2 ),
                           normal( rng(0), 0, uthi_I2 ),
                           normal( rng(0), 0, uthi_I2 ), qi_I2, 0, 0 );
          repeat( Z_I2 ) {

            // Leading order correction to Juttner distribution; also fixes uthe
            // to account for weakly relativistic effects

            do {
              vx = normal( rng(0), 0, uthe*uthe_relativistic_correction );
              vy = normal( rng(0), 0, uthe*uthe_relativistic_correction );
              vz = normal( rng(0), 0, uthe*uthe_relativistic_correction );
              u2 = vx*vx + vy*vy + vz*vz;
            } while ( uniform( rng(0), 0, 1 ) < u2*u2*one_over_eight_theta );
            inject_particle( electron, x, y, z, vx, vy, vz, -qe_2, 0, 0 );

          }
        } // if
      } // repeat

    } else { // qi_I1 <= qi_I2 -- not allowed!

      sim_log("ERROR: PARTICLE LOAD SEMANTICS REQUIRE qi_I1 > qi_I2");
      exit(0);

    } // if qi_I1 > qi_I2

  } // if load particles
#endif

// Pade fit - high-Z profile

#define I1_PROFILE_PROBABILITY(xi)                                                   \
  ( 1.03642 *                                                                        \
    ( (0.5184201288791436 - 0.9377618161236364*xi + 0.574908984738253*pow(xi,2) -    \
     0.1546771982667282*pow(xi,3) + 0.048328624369960646*pow(xi,4) -                 \
     0.021607896487013795*pow(xi,5) + 0.003199315832322355*pow(xi,6)) /              \
    (1 - 1.0598281629210033*xi + 0.4034701980854526*pow(xi,2) -                      \
     0.16656733984429437*pow(xi,3) + 0.07810662218354122*pow(xi,4) -                 \
     0.014952297721801492*pow(xi,5) +  0.003891946432949665*pow(xi,6))               \
   + 0.000200856 /* <-- offset to maintain positivity */  ) )

// Pade fit - low-Z profile

#define I2_PROFILE_PROBABILITY(xi)                                                   \
  ( fabs(qi_I2 / qi_I1) * 2.0 *                                                      \
    (0.23118956282383712 - 0.04366718965766797*xi - 0.09636585279725629*pow(xi,2) -  \
     0.003080684826052242*pow(xi,3) + 0.013994021713531602*pow(xi,4) +               \
     0.0037279452459577454*pow(xi,5) + 0.00028706869350233273*pow(xi,6)) /           \
   (1 - 1.059828056957553*xi + 0.4034701206993297*pow(xi,2) -                        \
     0.16656731568727687*pow(xi,3) + 0.07810660441748166*pow(xi,4) -                 \
     0.014952288878112218*pow(xi,5) + 0.0038919434427708804*pow(xi,6)) )

// Pade fit - low-Z mass flux function (see Mathematica worksheet "interfacex.nb"
// and notes from 6 Aug. 2015) This term is equal to:
//   a11(xi) * ( d log x/dxi + (z(xi)-x(xi))(1 + (Ai/AI)*((ZI+1)/(Zi+1)) d (ai/zi)/dxi )
// and to be turned into a speed, must be multiplied by:
//   (mI/mi)(1/ZI^2) Dp (T/1 keV)^5/2 Ai^(-0.5) x0^-1
// to get mass flux and, therefore, divided by rho c to get speed in simulation units
// (vi/c); to get heavy ion speed, set vI = -vi (mi/mI)

#define I2_MASS_FLUX_FUNCTION(xi)                                                    \
  (0.1542759226444911 + 0.03788536348380374*xi + 0.049021909111038944*pow(xi,2) -    \
   0.027055568725443978*pow(xi,3) - 0.02666788955260103*pow(xi,4) -                  \
     0.0059218828537447215*pow(xi,5) - 0.00040972771099082993*pow(xi,6)) /           \
   (1 + 1.5040184829931*xi + 1.1005467498678096*pow(xi,2) +                          \
    0.44355972952928274*pow(xi,3) + 0.10151805763024849*pow(xi,4) +                  \
    0.012694566906927567*pow(xi,5) +  0.0006930803042010016*pow(xi,6))

// Total mass density variation as a function of xi, normalized to rho0

#define NORMALIZED_MASS_DENSITY(xi)                                                  \
  ( 0.7495077162996403 - 0.5186071721292596*xi + 0.21278251708066112*pow(xi,2) -     \
    0.0938868243001879*pow(xi,3) + 0.02887440490448864*pow(xi,4) ) /                 \
  ( 1 - 0.44326556385099014*xi + 0.2013090643184375*pow(xi,2) -                      \
    0.11084677940449869*pow(xi,3) + 0.028403286204103685*pow(xi,4) )

#if 0
  // Load particles
  if ( load_particles ) {
    sim_log("Loading particles.");
    // Fast load of particles--don't bother fixing artificial domain correlations
    double xmin=grid->x0, xmax=grid->x1;
    double ymin=grid->y0, ymax=grid->y1;
    double zmin=grid->z0, zmax=grid->z1;

    // Kim's load similarity variable
    double xi;


    // Modify loading so that we use uniformly weighted particles throughout.
    // Conditional test on qi_I1 vs. qi_I2 is done so we can avoid putting
    // a conditional in the inner loop.

    if ( fabs(qi_I1) > fabs(qi_I2) ) {

      repeat( Ne/(topology_x*topology_y*topology_z) ) {
        double x = uniform_rand( xmin, xmax );
        double y = uniform_rand( ymin, ymax );
        double z = uniform_rand( zmin, zmax );

        // Initial ion diffusion speeds - omit initial electron speeds for now i
        //                                (they carry little momentum and ve_diffusion << vthe)
        double vx1_diffusion = 0;            // Diffusion speed of ion species 1 in simulation units (v/c)
        double vx2_diffusion = 0;            // Diffusion speed of ion species 2 in simulation units (v/c)
        double rho;

        if ( iv_region ) continue;           // Particle fell in iv_region.  Don't load.

        // Similarity variable
        xi =  ( ((Lx - x) - (Lx/2.0)) / Lx ); // position (flipped about box center in x) in units of fractions of half box size (dimensionless)
        xi *= half_box_size_in_xi;            // position (in units of similarity variable xi)


        if ( xi < -5.0 ) {                               // region 1 - left of layer
          if ( uniform_rand( 0, 1 ) < I1_PROFILE_PROBABILITY(-5.0) ) {
            inject_particle( electron, x, y, z,
                             maxwellian_rand(uthe),
                             maxwellian_rand(uthe),
                             maxwellian_rand(uthe), qe_1, 0, 0 );
            inject_particle( ion_I1, x, y, z,
                             maxwellian_rand(uthi_I1) + vx1_diffusion,
                             maxwellian_rand(uthi_I1),
                             maxwellian_rand(uthi_I1), qi_I1, 0, 0 );
          }
        } else if ( xi > 1.65 ) {                        // region 2 - right of layer
          if ( uniform_rand( 0, 1 ) < I2_PROFILE_PROBABILITY(1.65) ) {
            inject_particle( electron, x, y, z,
                             maxwellian_rand(uthe),
                             maxwellian_rand(uthe),
                             maxwellian_rand(uthe), qe_1, 0, 0 );
            inject_particle( ion_I2, x, y, z,
                             maxwellian_rand(uthi_I2) + vx2_diffusion,
                             maxwellian_rand(uthi_I2),
                             maxwellian_rand(uthi_I2), qi_I1, 0, 0 );
          } // if
        } else {                                        // overlap region

          // Compute initial drift speeds of ion species 1 and 2
          rho = rho0 * NORMALIZED_MASS_DENSITY(xi);

          // Aluminum is on the right, D on the left, so diffusion flux of
          // D (species 2) should be positive and diffusion flux of Al (species 1)
          // should be negative in sign
          if ( xi < 1.570 ) {

            vx2_diffusion = + ( A_I1 / (Z_I1*Z_I1*A_I2) )
                              * D_proton
                              * pow( t_e/1000., 2.5 )
                              * pow( A_I2, -0.5 )
                              * I2_MASS_FLUX_FUNCTION(xi)
                              / ( rho * speed_of_light * delta_x0 );

            // If mass fluxes are equal, then
            //   m_i Gamma_i = -m_I Gamma_I
            // however, Gamma_i = n_i v_i and Gamma_I = n_I v_I, so this implies
            //   v_I = - (m_i/m_I) [ Z_i n_i / Z_I n_I ] (Z_I / Z_i) v_i
            // And note that term in [] is the ratio of the PROFILE_PROBABILITY values

            // N.b. the conditional check here is because the Pade fit goes slightly
            // negative at large xi (and we wish to avoid div by zero).

            vx1_diffusion =  -   ( A_I2 / A_I1 )
                               * ( I2_PROFILE_PROBABILITY(xi) / I1_PROFILE_PROBABILITY(xi) )
                               * ( Z_I1 / Z_I2 )
                               * vx2_diffusion;

          } else {
            vx1_diffusion = vx2_diffusion = 0;     // vx2_diffusion fit goes slightly negative at xi > 1.57
          }

          if ( uniform_rand( 0, 1 ) < I1_PROFILE_PROBABILITY(xi) ) {
            inject_particle( electron, x, y, z,
                             maxwellian_rand(uthe),
                             maxwellian_rand(uthe),
                             maxwellian_rand(uthe), qe_1, 0, 0 );
            inject_particle( ion_I1, x, y, z,
                             maxwellian_rand(uthi_I1) + vx1_diffusion,
                             maxwellian_rand(uthi_I1),
                             maxwellian_rand(uthi_I1), qi_I1, 0, 0 );
          } // if

          if ( uniform_rand( 0, 1 ) < I2_PROFILE_PROBABILITY(xi) ) {
            inject_particle( electron, x, y, z,
                             maxwellian_rand(uthe),
                             maxwellian_rand(uthe),
                             maxwellian_rand(uthe), qe_1, 0, 0 );
            inject_particle( ion_I2, x, y, z,
                             maxwellian_rand(uthi_I2) + vx2_diffusion,
                             maxwellian_rand(uthi_I2),
                             maxwellian_rand(uthi_I2), qi_I1, 0, 0 );
          } // if
        } // else

      } // repeat

    } else { // qi_I1 <= qi_I2 -- not allowed!

      sim_log("ERROR: PARTICLE LOAD SEMANTICS REQUIRE qi_I1 > qi_I2");
      exit(0);

    } // else

  } // if
#endif

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

  // global->fdParams.format = band_interleave;
  // sim_log ( "Field output format          : band_interleave" );

  global->hedParams.format = band;
  sim_log ( "Electron hydro output format : band" );

  global->hI1dParams.format = band;
  sim_log ( "Hydrogen hydro output format : band" );

  global->hI2dParams.format = band;
  sim_log ( "Helium hydro output format   : band" );

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
  // //????????????????????????????????????????????????????????????????????????????
  int stride_x = 1, stride_y = 1, stride_z = 1;
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
  // I1 hydro

  // relative path to electron species data from global header
  sprintf(global->hI1dParams.baseDir, "I1hydro");

  // base file name for fields output
  sprintf(global->hI1dParams.baseFileName, "I1_hydro");

  // set I1 hydro strides
  global->hI1dParams.stride_x = stride_x;
  global->hI1dParams.stride_y = stride_y;
  global->hI1dParams.stride_z = stride_z;
  sim_log ( "Ion species x-stride " << global->hI1dParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hI1dParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hI1dParams.stride_z );

  // add I1 hydro parameters to list
  global->outputParams.push_back(&global->hI1dParams);

  //----------------------------------------------------------------------
  // I2 hydro

  // relative path to electron species data from global header
  sprintf(global->hI2dParams.baseDir, "I2hydro");

  // base file name for fields output
  sprintf(global->hI2dParams.baseFileName, "I2_hydro");

  // set I2 hydro strides
  global->hI2dParams.stride_x = stride_x;
  global->hI2dParams.stride_y = stride_y;
  global->hI2dParams.stride_z = stride_z;
  sim_log ( "Ion species x-stride " << global->hI2dParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hI2dParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hI2dParams.stride_z );

  // add I2 hydro parameters to list
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
  global->hI1dParams.output_variables(  current_density  | charge_density |
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

  create_hydro_list(varlist, global->hI1dParams);
  sim_log ( "Ion species variable list: " << varlist );

  create_hydro_list(varlist, global->hI2dParams);
  sim_log ( "Ion species variable list: " << varlist );

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

  //if ( step()%200==0 ) sim_log("Time step: "<<step());
  sim_log("Time step: "<<step());

sim_log( "Beginning diagnostics block." );
#if DEBUG_SYNCHRONIZE_IO
  // DEBUG
  barrier();
  sim_log( "Beginning diagnostics block.  Synchronized." );
#endif

# define should_dump(x) \
  (global->x##_interval>0 && remainder(step(),global->x##_interval)==0)

  // Do a mkdir at time t=0 to ensure we have all the directories we need
  if ( step()==0 && rank()==0 ) {
    sim_log("Making output directories.");
    dump_mkdir("rundata");
    dump_mkdir("fft");
    dump_mkdir("field");
    dump_mkdir("ehydro");
    dump_mkdir("I1hydro");
    dump_mkdir("I2hydro");
    dump_mkdir("restart");
    dump_mkdir("particle");
    dump_mkdir("velocity");
    sim_log("Finished making output directories.");

    // Turn off rundata for now
    // dump_grid("rundata/grid");
    // dump_materials("rundata/materials");
    // dump_species("rundata/species");
    global_header("global", global->outputParams);

  }

  // Energy dumps store all the energies in various directions of E and B
  // and the total kinetic (not including rest mass) energies of each species
  // species in a simple text format. By default, the energies are appended to
  // the file. However, if a "0" is added to the dump_energies call, a new
  // energies dump file will be created. The energies are in the units of the
  // problem and are all time centered appropriately. Note: When restarting a
  // simulation from a restart dump made at a prior time step to the last
  // energies dump, the energies file will have a "hiccup" of intervening
  // time levels. This "hiccup" will not occur if the simulation is aborted
  // immediately following a restart dump. Energies dumps are in a text
  // format and the layout is documented at the top of the file. Only rank 0
  // makes makes an energies dump.
  if( should_dump(energies) ) dump_energies( "rundata/energies", step()==0 ? 0 : 1 );

  // Field and hydro data

  if ( should_dump(field) ) {
    field_dump( global->fdParams );
//  dump_fields( "field/fields", (int)step() );

#if 1
    if ( global->load_particles ) {
      hydro_dump( "electron", global->hedParams );
      hydro_dump( "I1",       global->hI1dParams );
      hydro_dump( "I2",       global->hI2dParams );
    }
#endif

  }

  // Particle dump data
#if 0
  if ( should_dump(particle) && global->load_particles ) {
    dump_particles( "electron", "particle/eparticle" );
    if ( global->mobile_ions ) {
      if ( global->H_present  ) dump_particles( "H",  "particle/Hparticle" );
      if ( global->He_present ) dump_particles( "He", "particle/Heparticle" );
    }
  }
#endif


#define ALLOCATE(A,LEN,TYPE)                                             \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate."));

#if DEBUG_SYNCHRONIZE_IO
  // DEBUG
  barrier();
  sim_log( "setup done; going to standard VPIC output." );
#endif

#if DEBUG_SYNCHRONIZE_IO
  // DEBUG
  barrier();
  sim_log( "standard VPIC output done; going to FFT." );
#endif

//????????????????????????????????????????????????????????????????????????????
#if 0
// kx fft diag.
  // ------------------------------------------------------------------------
  // Custom diagnostic where we write out Ex for each cell in grid.
  // These data are stored for every fft_ex_interval time step.
  //
  // Refactored for 2D LPI problem.
# define FFT_HEADER_SIZE (sizeof(int))

# define WRITE_FFT(SUFF,INTERVAL)                                              \
  BEGIN_PRIMITIVE {                                                            \
    status = fileIO_##SUFF.open( fname_##SUFF, io_read_write);                 \
    if ( status==fail ) ERROR(("Could not open file."));                       \
    fileIO_##SUFF.seek( uint64_t( FFT_HEADER_SIZE                              \
                                  + uint64_t((step()/INTERVAL)*stride*sizeof(float)) ),  \
                        SEEK_SET );                                            \
    fileIO_##SUFF.write( SUFF, stride );                                       \
    fileIO_##SUFF.close();                                                     \
  } END_PRIMITIVE

# define WRITE_FFT_HEADER(SUFF)                                                \
  BEGIN_PRIMITIVE {                                                            \
    status = fileIO_##SUFF.open( fname_##SUFF, io_write);                      \
    if ( status==fail ) ERROR(("Could not open file."));                       \
    fileIO_##SUFF.write( &grid->nx, 1 );                                       \
    fileIO_##SUFF.close();                                                     \
  } END_PRIMITIVE

  static int initted=0;
  static float *ex;
  static char fname_ex[256];
  FileIO       fileIO_ex;
  FileIOStatus status;
  int stride=grid->nx;

  if ( !initted ) {
    // Allocate space for data arrays
    long array_length=grid->nx;
    ALLOCATE(ex,  array_length, float);
    // Define filenames
    sprintf( fname_ex,  "fft/fft_ex.%i",  (int)rank() );
    // On first timestep prepend a header with number of x meshpoints to each file
    if ( step()==0 ) {
      WRITE_FFT_HEADER(ex);
    }
    initted=1;
  }

  // *** Ex ***
  if ( !(step()%global->fft_ex_interval) ) {
    // Store data into array ex
    for ( int i=0; i<grid->nx; ++i ) {
      int k=INDEX_FORTRAN_3(i+1,1,grid->nz/2+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
      ex[i]=field[k].ex;
    }
    // Write array to file
    sim_log("Writing FFT data for Ex fields.");
    DIAG_LOG("Startign to dump FFT Ex data.");
    WRITE_FFT(ex, global->fft_ex_interval);
    DIAG_LOG("Finished dumping FFT Ex data.");
  }
#endif

#if DEBUG_SYNCHRONIZE_IO
  // DEBUG
  barrier();
  sim_log( "FFT done; going to velocity-space diagnostic." );
#endif

//????????????????????????????????????????????????????????????????????????????
#if 0
  // -------------------------------------------------------------------------
  // Diagnostic to write a 2d vx, vz binned velocity distribution
  // Note that we have converted, using the relativistic gammas, from momentum
  // (which the code uses) to velocity prior to binning.
  //
#define NVX 100
#define NVZ 100
/* Set PZMASK to an appropriate value */
// PZMAX is the distance in simulation units on each side from line z = fixed value in sim units
#define PZMAX 10
// distance in sim units
#define PZCENTER (0)  // at 0 micron
#define PMASK ((pz-PZCENTER)*(pz-PZCENTER)<PZMAX*PZMAX)
#define PREPARE_VELOCITY_SPACE_DATA(SUFF, NAME)                           \
    {                                                                     \
      species_t *sp;                                                      \
      for (int i=0; i<NVX; ++i)                                           \
        for (int j=0; j<NVZ; ++j)                                         \
	  f_##SUFF[i*NVZ+j]=0;                                            \
      sp = find_species_name(NAME, species_list);                         \
      for (int ip=0; ip<sp->np; ip++) {                                   \
	particle_t *p=&sp->p[ip];                                         \
        /* Lots of stuff commented because PMASK only has pz */           \
	int nxp2=grid->nx+2;                                              \
	int nyp2=grid->ny+2;                                              \
	/* int nzp2=grid->nz+2;    */                                     \
	/* Turn index i into separate ix, iy, iz indices */               \
	int iz = p->i/(nxp2*nyp2);                                        \
	/* int iy = (p->i-iz*nxp2*nyp2)/nxp2;  */                         \
	/* int ix = p->i-nxp2*(iy+nyp2*iz); */                            \
	/* Compute real particle position from relative coords and grid data */ \
	/* double px = grid->x0+((ix-1)+(p->dx+1)*0.5)*grid->dx; */       \
	/* double py = grid->y0+((iy-1)+(p->dy+1)*0.5)*grid->dy; */       \
	double pz = grid->z0+((iz-1)+(p->dz+1)*0.5)*grid->dz;             \
	float invgamma=1/sqrt(1+p->ux*p->ux+p->uy*p->uy+p->uz*p->uz);     \
	float vx=p->ux*grid->cvac*invgamma;                               \
	float vz=p->uz*grid->cvac*invgamma;                               \
	long ivx=long(vx/dvx_##SUFF+(NVX/2));                             \
	long ivz=long(vz/dvz_##SUFF+(NVZ/2));                             \
	if ( abs(ivx)<NVX && abs(ivz)<NVZ && PMASK ) f_##SUFF[ivx*NVZ+ivz]+=p->q;  \
      }                                                                   \
    }

#define INCLUDE_VELOCITY_HEADER 0
#if INCLUDE_VELOCITY_HEADER
#  define VELOCITY_HEADER_SIZE (2*sizeof(int)+2*sizeof(float))
#  define WRITE_VELOCITY_HEADER(SUFF)                               \
    {                                                               \
      int nvpts[2] = { NVX, NVZ };                                  \
      float dv[2];                                                  \
      dv[0] = dvx_##SUFF; dv[1] = dvz_##SUFF;                       \
      status = fileIO_##SUFF.open( fname_##SUFF, io_write );        \
      if ( status==fail ) ERROR(("Could not open file."));          \
      fileIO_##SUFF.write( &nvpts, 2 );                             \
      fileIO_##SUFF.write( &dv,    2 );                             \
      fileIO_##SUFF.close();                                        \
    }
#else
#  define VELOCITY_HEADER_SIZE 0
#  define WRITE_VELOCITY_HEADER(SUFF)
#endif

// NOTE:  WE DO NOT WRITE ON TIME STEP 0! (see -1 in seek())
#define DUMP_VELOCITY_DATA(SUFF,LEN,HEADER_SIZE)                    \
    {                                                               \
      status = fileIO_##SUFF.open( fname_##SUFF,                    \
                                   ( ( INCLUDE_VELOCITY_HEADER==0 && \
                                       step()==global->velocity_interval ) \
                                     ? io_write                     \
                                     : io_read_write ) );           \
      if ( status==fail ) ERROR(("Could not open file."));          \
      fileIO_##SUFF.seek( uint64_t( HEADER_SIZE +                   \
                                    (step()/global->velocity_interval-1) \
                                    * LEN * sizeof(float)) ,        \
                          SEEK_SET );                               \
      fileIO_##SUFF.write( f_##SUFF, LEN );                         \
      fileIO_##SUFF.close();                                        \
    }

  if ( 1 ) {
    float f_e[NVX*NVZ], f_I2[NVX*NVZ], f_I1[NVX*NVZ];
    float vmax_e  = 10*global->vthe,    dvx_e,  dvz_e;
    float vmax_I1  = 10*global->vthi_I1,  dvx_I1,  dvz_I1;
    float vmax_I2 = 10*global->vthi_I2, dvx_I2, dvz_I2;
    FileIO       fileIO_e, fileIO_I1, fileIO_I2;
    FileIOStatus status;
    static char fname_e[256], fname_I1[256], fname_I2[256];
    dvx_e  = dvz_e  = 2*vmax_e /NVX;
    dvx_I1 = dvz_I1 = 2*vmax_I1/NVX;
    dvx_I2 = dvz_I2 = 2*vmax_I2/NVX;
    sprintf( fname_e,  "velocity/velocity_e.%i",  (int)rank() );
    sprintf( fname_I1, "velocity/velocity_I1.%i",  (int)rank() );
    sprintf( fname_I2, "velocity/velocity_I2.%i", (int)rank() );
    if ( !step() ) {
      WRITE_VELOCITY_HEADER(e);
      WRITE_VELOCITY_HEADER(I1);
      WRITE_VELOCITY_HEADER(I2);
    }

    // NOTE: We don't write on time step 0, as per comment above.
    if ( step()!=0 && (step()%global->velocity_interval)==0 ) {
      PREPARE_VELOCITY_SPACE_DATA(e, "electron");
      PREPARE_VELOCITY_SPACE_DATA(I1, "I1");
      PREPARE_VELOCITY_SPACE_DATA(I2, "I2");
      DUMP_VELOCITY_DATA(e,  NVX*NVZ, VELOCITY_HEADER_SIZE);
      DUMP_VELOCITY_DATA(I1, NVX*NVZ, VELOCITY_HEADER_SIZE);
      DUMP_VELOCITY_DATA(I2, NVX*NVZ, VELOCITY_HEADER_SIZE);
    }

  }

#endif

#if DEBUG_SYNCHRONIZE_IO
  // DEBUG
  barrier();
  sim_log( "Velocity-space diagnostic done; going to restart." );
#endif

  //---------------------------------------------------------------------------

  // Restart dump filenames are by default tagged with the current timestep.
  // If you do not want to do this add "0" to the dump_restart arguments. One
  // reason to not tag the dumps is if you want only one restart dump (the most
  // recent one) to be stored at a time to conserve on file space. Note: A
  // restart dump should be the _very_ _last_ dump made in a diagnostics
  // routine. If not, the diagnostics performed after the restart dump but
  // before the next timestep will be missed on a restart. Restart dumps are
  // in a binary format. Each rank makes a restart dump.

  // Restarting from restart files is accomplished by running the executable
  // with "restart restart" as additional command line arguments.  The executable
  // must be identical to that used to generate the restart dumps or else
  // the behavior may be unpredictable.

  // Note:  restart is not currently working with custom boundary conditions
  // (such as the reflux condition) and has not been tested with emission
  // turned on.

#if 0
  if ( step() && should_dump(restart) ) {
    if ( !global->rtoggle ) {
      DIAG_LOG("Starting restart0 dump");
      global->rtoggle=1;
      dump_restart("restart/restart0",0);
      DIAG_LOG("Restart dump restart0 completed.");
    } else {
      DIAG_LOG("Starting restart1 dump");
      global->rtoggle=0;
      dump_restart("restart/restart1",0);
      DIAG_LOG("Restart dump restart1 completed.");
    }
  }
#endif

  if ( step()>0 && should_dump(restart) ) {
//if ( step()>0 && remainder(step(),2500)==0 ) {
    static const char * restart_fbase[2] = { "restart/restart0", "restart/restart1" };
    double dumpstart = uptime();

    //dump_restart( restart_fbase[global->rtoggle], 0 );
    checkpt( restart_fbase[global->rtoggle], 0 );

    global->rtoggle^=1;
  }


  // Shut down simulation when wall clock time exceeds global->quota_sec.
  // Note that the mp_elapsed() is guaranteed to return the same value for all
  // processors (i.e., elapsed time on proc #0), and therefore the abort will
  // be synchronized across processors.

  if ( step()>0 && global->quota_check_interval && (step()%global->quota_check_interval)==0 ) {
    if ( uptime() > global->quota_sec ) {

      checkpt( "restart/restart", 0 );

      sim_log( "Restart dump restart completed." );
      sim_log( "Allowed runtime exceeded for this job.  Terminating." );
      barrier(); // Just to be safe
      finalize();
      exit(0);
    }
  }

#if DEBUG_SYNCHRONIZE_IO
  // DEBUG
  barrier();
  sim_log( "All diagnostics done in begin_diagnostics" );
#endif

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

/* Datatypes required by collision algorithm. */

/* Key datatype, used for shuffling particles */
typedef struct _pkey_t {
  int i;           /* Particle array index */
} pkey_t;

typedef struct _cdata_t {
  species_t *s;    /* Particle species */
  float     cvar;  /* Base variance of polar angle--see T&A paper */
  float     m;     /* Particle mass / electron mass */
  float     z;     /* Particle core charge */
} cdata_t;

/* Takizuka and Abe model for 2-particle elastic scattering from JCP 25, 205 (1977).
       I, J:    particle indices
       P1, P2:  particle array heads
       M1, M2:  relative masses of particles
       var:     variance of collision operator */

#define SCATTER_PARTICLES(I,J,P1,P2,M1,M2,SQRT_VAR)                                        \
    do {                                                                                   \
      float ux, uy, uz, uperp, u, dux, duy, duz, mratio1, mratio2;                         \
      float sin_theta, one_m_cos_theta, sin_phi, cos_phi, delta, phi;                      \
                                                                                           \
      mratio1=(float)M2/(float)(M1+M2);                                                    \
      mratio2=(float)M1/(float)(M1+M2);                                                    \
      ux=(P1)[I].ux-(P2)[J].ux;                                                            \
      uy=(P1)[I].uy-(P2)[J].uy;                                                            \
      uz=(P1)[I].uz-(P2)[J].uz;                                                            \
      uperp=sqrt(ux*ux+uy*uy);                                                             \
      u=sqrt(uperp*uperp+uz*uz);                                                           \
      delta=normal( rng(0), 0, SQRT_VAR )*pow(u,-1.5);                                     \
      sin_theta=2*delta/(1+delta*delta);                                                   \
      one_m_cos_theta=sin_theta*delta;                                                     \
      phi=2*M_PI*uniform( rng(0), 0, 1 );                                                  \
      sin_phi=sin(phi);                                                                    \
      cos_phi=cos(phi);                                                                    \
      if ( uperp>0 ) {  /* General case */                                                 \
        dux=((ux*uz*sin_theta*cos_phi)-uy*u*sin_theta*sin_phi)/uperp - ux*one_m_cos_theta; \
        duy=((uy*uz*sin_theta*cos_phi)+ux*u*sin_theta*sin_phi)/uperp - uy*one_m_cos_theta; \
        duz=-uperp*sin_theta*cos_phi-uz*one_m_cos_theta;                                   \
      } else {          /* Handle purely z-directed difference vectors separately */       \
        dux=u*sin_theta*cos_phi;                                                           \
        duy=u*sin_theta*sin_phi;                                                           \
        duz=-u*one_m_cos_theta;                                                            \
      }                                                                                    \
      (P1)[I].ux+=mratio1*dux;                                                             \
      (P1)[I].uy+=mratio1*duy;                                                             \
      (P1)[I].uz+=mratio1*duz;                                                             \
      (P2)[J].ux-=mratio2*dux;                                                             \
      (P2)[J].uy-=mratio2*duy;                                                             \
      (P2)[J].uz-=mratio2*duz;                                                             \
    } while (0);

begin_particle_collisions {

  // Put collision model here since particle injection occurs after the particle
  // sort (collision models require sorted particle array).

  // Eventually this needs to be integrated into the VPIC source tree, but for now
  // leave this in the input deck.

  if ( global->do_collisions && global->load_particles && step()>0 && step()%global->tstep_coll==0 ) {
    species_t *s1, *s2;
    int icell, i, j, npar1, npar2, isp1, isp2;
    static int initted=0, nsp;
    static cdata_t cd[3 /* =nsp */];
    static pkey_t *karr1, *karr2;
    float var, sqrt_var, sqrt_half_var;
    particle_t *phead1, *phead2;

    sim_log("Sample collision operator: Begin");

    if ( initted==0 ) {
      karr1=(pkey_t *)malloc( (size_t)global->nppc_max*sizeof(pkey_t) );
      if ( !karr1 ) ERROR(("Could not allocate key array in collision handler."));
      karr2=(pkey_t *)malloc( (size_t)global->nppc_max*sizeof(pkey_t) );
      if ( !karr2 ) ERROR(("Could not allocate key array in collision handler."));

      isp1=0;
      cd[isp1].s    = NULL;
      cd[isp1].s    = find_species( "electron" );
      if ( cd[isp1].s==NULL ) ERROR(("Could not find electron species."));
      cd[isp1].m    = 1;
      cd[isp1].z    = global->Z_e;
      cd[isp1].cvar = global->cvar*pow(cd[isp1].z,4)/pow(cd[isp1].m*0.5,2);
      sim_log("cvar = "<<cd[isp1].cvar);

      ++isp1;
      cd[isp1].s    = NULL;
      cd[isp1].s    = find_species( "I1" );
      if ( cd[isp1].s==NULL ) ERROR(("Could not find I1 species."));
      cd[isp1].m    = global->mime_I1;  /* in m_e */
      cd[isp1].z    = global->Z_I1;               /* in e */
      cd[isp1].cvar = global->cvar*pow(cd[isp1].z,4)/pow(cd[isp1].m*0.5,2);
      sim_log("cvar = "<<cd[isp1].cvar);

      ++isp1;
      cd[isp1].s    = NULL;
      cd[isp1].s    = find_species( "I2" );
      if ( cd[isp1].s==NULL ) ERROR(("Could not find I2 species."));
      cd[isp1].m    = global->mime_I2; /* in m_e */
      cd[isp1].z    = global->Z_I2;               /* in e */
      cd[isp1].cvar = global->cvar*pow(cd[isp1].z,4)/pow(cd[isp1].m*0.5,2);
      sim_log("cvar = %e"<<cd[isp1].cvar);

      nsp=isp1+1;
      initted=1;
    }


    // Self-scattering among species
    for ( isp1=0; isp1<nsp; ++isp1 ) {
      sim_log("Self-Scatter of speciees "<<isp1);
      s1=cd[isp1].s;
      for ( icell=((grid->nx+2)*(grid->ny+2)*(grid->nz+2))-2; icell>=0; --icell ) {

        npar1  = s1->partition[icell+1]-s1->partition[icell];
        phead1 = s1->p;

        var=cd[isp1].cvar * (double)npar1 * global->inv_nppc;  // multiplying by density as per Bill's mods to get density-variation of collision operator
        sqrt_var=sqrt(var);
        sqrt_half_var=sqrt(var/2);


        // Fischer-Yates shuffle on key array to randomize the order of particles in a cell
        for ( i=0; i<npar1; ++i ) {
          int j;
          j = uniform( rng(0), 0, (double)i+1 );                  // j is a random integer satisfying 0 <= j <= i
          if ( j != i ) karr1[i].i = karr1[j].i;
          karr1[j].i = s1->partition[icell] + i;
        }


        if ( npar1>4 ) {
          for ( i=0; i<npar1-4; i+=2 )
            SCATTER_PARTICLES( karr1[i  ].i, karr1[i+1].i, phead1, phead1, 1.0, 1.0, sqrt_var );
          if ( i==npar1-4 ) {
            SCATTER_PARTICLES( karr1[i  ].i, karr1[i+1].i, phead1, phead1, 1.0, 1.0, sqrt_var );
            SCATTER_PARTICLES( karr1[i+2].i, karr1[i+3].i, phead1, phead1, 1.0, 1.0, sqrt_var  );
          } else {
            SCATTER_PARTICLES( karr1[i  ].i, karr1[i+1].i, phead1, phead1, 1.0, 1.0, sqrt_half_var );
            SCATTER_PARTICLES( karr1[i  ].i, karr1[i+2].i, phead1, phead1, 1.0, 1.0, sqrt_half_var );
            SCATTER_PARTICLES( karr1[i+1].i, karr1[i+2].i, phead1, phead1, 1.0, 1.0, sqrt_half_var );
          }
        }
      }
    }

    // Cross-species scattering

    if ( global->self_collisions_only==0 ) {
      for ( isp1=1; isp1<nsp; ++isp1 ) {
        s1=cd[isp1].s;
        for ( isp2=0; isp2<isp1; ++isp2 ) {
          sim_log("Cross-Scatter of species "<<isp1<<" and "<<isp2);
          float m1=cd[isp1].m, m2=cd[isp2].m;
          float rmass=m1*m2/(m1+m2);

          s2=cd[isp2].s;


          for ( icell=((grid->nx+2)*(grid->ny+2)*(grid->nz+2))-2; icell>=0; --icell ) {

            npar1  = s1->partition[icell+1]-s1->partition[icell];
            npar2  = s2->partition[icell+1]-s2->partition[icell];
            phead1 = s1->p;
            phead2 = s2->p;

            // multiply by density as per Bill's mods to get density-variation collisions
            double npar_min = ( npar1 < npar2 ? npar1 : npar2 );
            var=( global->cvar*pow(cd[isp1].z,2)*pow(cd[isp2].z,2)/(rmass*rmass) ) * npar_min * global->inv_nppc;
            sqrt_var=sqrt(var);

            // Fischer-Yates shuffle on key array to randomize the order of particles in a cell
            for ( i=0; i<npar1; ++i ) {
              int j;
              j = uniform( rng(0), 0, (double)i+1 );                  // j is a random integer satisfying 0 <= j <= i
              if ( j != i ) karr1[i].i = karr1[j].i;
              karr1[j].i = s1->partition[icell] + i;
            }
            for ( i=0; i<npar2; ++i ) {
              int j;
              j = uniform( rng(0), 0, (double)i+1 );                  // j is a random integer satisfying 0 <= j <= i
              if ( j != i ) karr2[i].i = karr2[j].i;
              karr2[j].i = s2->partition[icell] + i;
            }

            // Call this with npar1 >= npar2

#           define CROSS_COLLISION_LOOP(NPAR1,NPAR2,P1,P2,K1,K2,M1,M2,SQRT_VAR)  \
              do {                                                               \
                if ( NPAR1>0 && NPAR2>0 ) {                                      \
                  int ii, iimax, index1=-1, index2=-1;                           \
                                                                                 \
                  ii=(NPAR1)/(NPAR2);                                            \
                  iimax=(NPAR1)-(NPAR2)*ii;                                      \
                  for ( i=0; i<iimax; ++i) {                                     \
                    ++index2;                                                    \
                    for ( j=0; j<ii+1; ++j ) {                                   \
                      ++index1;                                                  \
                      SCATTER_PARTICLES( K1[index1].i, K2[index2].i, P1, P2, M1, M2, (SQRT_VAR)/sqrt(ii+1) ); \
                    }                                                            \
                  }                                                              \
                  for ( ; i<(NPAR2); ++i ) {                                     \
                    ++index2;                                                    \
                    for ( j=0; j<ii; ++j ) {                                     \
                      ++index1;                                                  \
                      SCATTER_PARTICLES( K1[index1].i, K2[index2].i, P1, P2, M1, M2, (SQRT_VAR)/sqrt(ii) );  \
                    }                                                            \
                  }                                                              \
                }                                                                \
              } while (0)

            if ( npar1>npar2 ) CROSS_COLLISION_LOOP( npar1, npar2, phead1, phead2, karr1, karr2, m1, m2, sqrt_var );
            else               CROSS_COLLISION_LOOP( npar2, npar1, phead2, phead1, karr2, karr1, m2, m1, sqrt_var );
          }
        }
      }
    }
    sim_log("Sample collision operator: Finish");
  }

}
