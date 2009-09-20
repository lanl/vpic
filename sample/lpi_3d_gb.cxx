//========================================================================
//
// LPI 3D deck - Linearly polarized (in y) plane wave incident from left 
//               boundary 
//
// Adapted from Albright's Lightning 3D LPI deck.  
// B. Albright, X-1-PTA;      28 Jan. 2007
// 
// Only poynting data output (to stdout).  
// Restart, quota checks, divergence cleaning turned off. 
//========================================================================

begin_globals {
  double e0;                   // peak amplitude of oscillating electric field
  double omega;                // angular freq. of the beam
  int field_interval;          // how often to dump fields and hydro
  int particle_interval;       // how often to dump particle data
  int poynting_interval;       // how often to compute poynting flux on boundary
  int restart_interval;        // how often to write restart data
  int quota_check_interval;    // how often to check if quote exceeded
  double quota_sec;            // run quota in sec
  int rtoggle;                 // enable save of last 2 restart files for safety
  int load_particles;          // were particles loaded? 
  double topology_x;           // domain topology needed to normalize Poynting diagnostic 
  double topology_y;
  double topology_z;
  int mobile_ions;
  int H_present; 
  int He_present; 
  int write_poynting_data;     // Whether to write poynting data to file (or just stdout)

  // Parameters for 3d Gaussian wave launch
  double lambda;               
  double waist;                // how wide the focused beam is
  double width;                
  double zcenter;              // center of beam at boundary in z 
  double ycenter;              // center of beam at boundary in y
  double xfocus;               // how far from boundary to focus
  double mask;                 // # gaussian widths from beam center where I is nonzero

};

begin_initialization {
  
  // System of units

  double ec         = 4.8032e-10;          // stat coulomb
  double c_vac      = 2.99792458e10;       // cm/sec
  double m_e        = 9.1094e-28;          // g
  double k_b        = 1.6022e-12;          // erg/eV
  double mec2       = m_e*c_vac*c_vac/k_b; 
  double mpc2       = mec2*1836.0; 

  double cfl_req    = 0.98;                // How close to Courant should we try to run
  double damp       = 0;                   // How much radiation damping
  double iv_thick   = 2;                   // Thickness of impermeable vacuum (in cells)

  // Experimental parameters

  double t_e               = 4000;         // electron temperature, eV
  double t_i               = 2000;         // ion temperature, eV
  double n_e_over_n_crit   = 0.14;         // n_e/n_crit
  double vacuum_wavelength = 351 * 1e-7;   // third micron light (cm)
  double laser_intensity   = 1e17 * 1e7;   // in ergs/cm^2 (note: 1 W = 1e7 ergs)

  // Simulation parameters 

  // Assume f/# = 4
  double Lx                = 35.0* 1e-4;   // In cm (note: 1 micron = 1e-4 cm)   
  double Ly                = 6.0 * 1e-4;                 
  double Lz                = 6.0 * 1e-4;                 

# if 1
  // Smaller problem for 48 core CU
  Lx /= 270;                                // 1/270 of the targeted problem on full RR 
  double nx                = 63;
  double ny                = 84;
  double nz                = 84; 
  double topology_x        = 3;
  double topology_y        = 4;
  double topology_z        = 4;           // 48 proc
#endif 

# if 0
  // Smaller problem for 720 core CU
  Lx /= 18;                                // 1/18 of the targeted problem on full RR 
  double nx                = 110;
  double ny                = 252;
  double nz                = 252; 
  double topology_x        = 5;
  double topology_y        = 12;
  double topology_z        = 12;           // 720 proc
#endif 

# if 0
  double nx                = 2048;
  double ny                = 256;
  double nz                = 256; 
  double topology_x        = 64;
  double topology_y        = 4;
  double topology_z        = 4;            // 1024 proc
# endif 


# if 0
  // Make box size smaller for fewer proc run:  for f/#=4
  Lx /= 4;  nx /= 4;
  Ly /= 2;  ny /= 2;   
  Lz /= 2;  nz /= 2; 

  topology_x = 64; 
  topology_y = 4;  
  topology_z = 4; 
# endif  


# if 0
  // Make box size smaller for test run
  Lx /= 4*32;  nx /= 4*32;
  Ly /= 2*2;   ny /= 2*2;   
  Lz /= 2*2;   nz /= 2*2; 

  topology_x = 2; 
  topology_y = 2;  
  topology_z = 2; 
# endif  


  // For GB:  Max performance = increase nppc until problem doesn't fit in memory anymore! 

  double nppc              = 25;           // Ave. number of particles/cell in ea. species
  int load_particles       = 1;            // Flag to turn on/off particle load 

  int mobile_ions          = 1;            // Whether or not to push ions
  int He_present=1, H_present=1;           // Initialize parameters
  double f_He              = 0.5;          // Ratio of number density of He to total ion density
  double f_H               = 1-f_He;       // Ratio of number density of H  to total ion density
  if ( f_He==1 ) H_present=0; 
  if ( f_He==0 ) He_present=0; 

  // Precompute some useful variables. 
  double A_H               = 1;
  double A_He              = 4;
  double Z_H               = 1;
  double Z_He              = 2; 
  double mic2_H            = mpc2*A_H;
  double mic2_He           = mpc2*A_He;
  double mime_H            = mic2_H /mec2; 
  double mime_He           = mic2_He/mec2; 

  double uth_e             = sqrt(t_e/mec2);      // vthe/c
  double uthi_H            = sqrt(t_i/mic2_H);    // vthi/c for H
  double uthi_He           = sqrt(t_i/mic2_He);   // vthi/c for He

  // Plasma skin deptth in cm
  double delta = (vacuum_wavelength / (2*M_PI) ) / sqrt( n_e_over_n_crit ); 

  double n_e   = c_vac*c_vac*m_e/(4*M_PI*ec*ec*delta*delta); // electron density in cm^-3
  double debye = uth_e*delta;                     // electron Debye length (cm)
  double omega = sqrt( 1/n_e_over_n_crit );       // laser beam freq. in wpe

  // Peak instantaneous E field in "natural units" 
  double e0    = sqrt( 2*laser_intensity / (m_e*c_vac*c_vac*c_vac*n_e) );  

  // Set up local mesh resolution and time step
  Lx /= delta;                                    // Convert box size to skin depths
  Ly /= delta;   
  Lz /= delta;   

  double hx = Lx/nx; 
  double hy = Ly/ny; 
  double hz = Lz/nz; 

  double cell_size_x       = hx*delta/debye;      // Cell size in Debye lengths
  double cell_size_y       = hy*delta/debye;
  double cell_size_z       = hz*delta/debye;

  double f_number          = 4;                   // f/# of beam
  double lambda            = vacuum_wavelength/delta; // vacuum wavelength in c/wpe
  double waist             = f_number*lambda;     // width of beam at focus in c/wpe
  double xfocus            = Lx/2;                // in c/wpe
  double ycenter           = 0;                   // center of spot in y on lhs boundary
  double zcenter           = 0;                   // center of spot in z on lhs boundary
  double mask              = 1.5;                 // set drive I=0 outside r>mask*width at lhs boundary
  double width = waist*sqrt( 1 + (lambda*xfocus/(M_PI*waist*waist))*(lambda*xfocus/(M_PI*waist*waist))); 
  e0                       = e0*(waist/width);    // at entrance (3D Gaussian) 
 
  double dt                = cfl_req*courant_length(Lx,Ly,Lz,nx,ny,nz); // in 1/wpe; n.b. c=1 in nat. units
  double nsteps_cycle      = trunc_granular(2*M_PI/(dt*omega),1)+1; 
  dt                       = 2*M_PI/omega/nsteps_cycle; // nsteps_cycle time steps in one laser cycle

  double t_stop            = 1001;                // Runtime in 1/wpe
  int particle_interval    = 0; 
  int poynting_interval    = 16;                  // Num. steps between dumping poynting flux
  int field_interval       = 0;                   // Num. steps between saving field, hydro data
  int restart_interval     = 0;                   // Num. steps between restart dumps
  int quota_check_interval = 100000;
  double quota_sec         = 12*3600;             // Run quota in sec. 
  int write_poynting_data  = 0;                   // Whether to write poynting data to file (or just stdout)

  double N_e               = nppc*nx*ny*nz;       // Number of macro electrons in box
  double Np_e              = Lx*Ly*Lz;            // "Number" of "physical" electrons in box (nat. units)
  double q_e               = -Np_e/N_e;           // Charge per macro electron
  double N_i               = N_e;                 // Number of macro ions of each species in box 
  double Np_i              = Np_e/(Z_H*f_H+Z_He*f_He); // "Number" of "physical" ions of each sp. in box
  double qi_H              = Z_H *f_H *Np_i/N_i;  // Charge per H  macro ion
  double qi_He             = Z_He*f_He*Np_i/N_i;  // Charge per He macro ion 
  
  // Print simulation parameters

  sim_log("***** Simulation parameters *****");
  sim_log("* Processors:                     "<<nproc());
  sim_log("* Topology:                       "<<topology_x<<" "<<topology_y<<" "<<topology_z); 
  sim_log("* Time step, max time, nsteps:    "<<dt<<" "<<t_stop<<" "<<int(t_stop/(dt))); 
  sim_log("* Debye length, XYZ cell sizes:   "<<debye<<" "<<cell_size_x<<" "<<cell_size_y<<" "<<cell_size_z);
  sim_log("* Real cell sizes (in Debyes):    "<<hx/uth_e<<" "<<hy/uth_e<<" "<<hz/uth_e);
  sim_log("* Lx, Ly, Lz =                    "<<Lx<<" "<<Ly<<" "<<Lz);
  sim_log("* nx, ny, nz =                    "<<nx<<" "<<ny<<" "<<nz);
  sim_log("* Charge/macro electron =         "<<q_e);
  sim_log("* Average particles/processor:    "<<N_e/nproc());
  sim_log("* Average particles/cell:         "<<nppc);
  sim_log("* Omega_0, Omega_pe:              "<<omega<<" "<<1);
  sim_log("* Plasma density, ne/nc:          "<<n_e<<" "<<n_e_over_n_crit);
  sim_log("* Vac wavelength (nm):            "<<vacuum_wavelength*1e7);
  sim_log("* I_laser (W/cm^2):               "<<laser_intensity/1e7);
  sim_log("* T_e, T_i (eV)                   "<<t_e<<" "<<t_i);
  sim_log("* m_e, m_H, m_He                  "<<"1 "<<mime_H<<" "<<mime_He);
  sim_log("* Radiation damping:              "<<damp);
  sim_log("* Fraction of courant limit:      "<<cfl_req);
  sim_log("* vthe/c:                         "<<uth_e);
  sim_log("* vthi_H /c:                      "<<uthi_H);
  sim_log("* vthe_He/c:                      "<<uthi_He);
  sim_log("* emax at entrance:               "<<e0);
  sim_log("* emax at waist:                  "<<e0/(waist/width));
  sim_log("* Poynting interval:              "<<poynting_interval); 
  sim_log("* field interval:                 "<<field_interval); 
  sim_log("* restart interval:               "<<restart_interval); 
  sim_log("* num vacuum edge grids:          "<<iv_thick);
  sim_log("* width, waist, xfocus:           "<<width<<" "<<waist<<" "<<xfocus); 
  sim_log("* ycenter, zcenter, mask:         "<<ycenter<<" "<<zcenter<<" "<<mask); 
  sim_log("*********************************");

  // Set up high level simulation parameters

  sim_log("Setting up high-level simulation parameters."); 
  num_step             = int(t_stop/(dt)); 

  verbose              = 0;
  num_comm_round       = 6;
  // DEBUG:  We need more info to find bugs.
  status_interval      = 20; // 2000;
  sync_shared_interval = status_interval; 
  clean_div_e_interval = status_interval; 
  clean_div_b_interval = status_interval; 

  // On cell, divergence cleaning (of E particularly) is very expensive.  Turn this
  // off for the GB runs. 
# if 0                                            
  sync_shared_interval = status_interval/10; 
  clean_div_e_interval = status_interval;           
  clean_div_b_interval = status_interval/10; 
# endif 
  
  global->e0                   = e0; 
  global->omega                = omega; 
  global->field_interval       = field_interval; 
  global->particle_interval    = particle_interval;
  global->poynting_interval    = poynting_interval;
  global->restart_interval     = restart_interval;
  global->quota_check_interval = quota_check_interval;
  global->quota_sec            = quota_sec;
  global->rtoggle              = 0;
  global->load_particles       = load_particles;
  global->mobile_ions          = mobile_ions; 
  global->H_present            = H_present; 
  global->He_present           = He_present; 
  global->topology_x           = topology_x;  
  global->topology_y           = topology_y;  
  global->topology_z           = topology_z;  
  global->xfocus               = xfocus;  
  global->ycenter              = ycenter; 
  global->zcenter              = zcenter; 
  global->mask                 = mask; 
  global->waist                = waist; 
  global->width                = width; 
  global->lambda               = lambda; 
  global->write_poynting_data  = write_poynting_data; 

  // Set up the species
  // Allow additional local particles in case of non-uniformity.

  // Note:  We have to adjust sort intervals for maximum performance on Cell. 

  sim_log("Setting up species."); 
  double max_local_np = 2.5*N_e/nproc(); 
  double max_local_nm = max_local_np/10;                       // Default is max_local_np/12.5  
  int ele_sort_int    = 200;
  int ion_sort_int    = ele_sort_int*20; 
  species_t * electron = define_species("electron", -1,  1, max_local_np, max_local_nm, ele_sort_int, 1); 

  // Start with two ion species.  We have option to go to Xe and Kr gas fills if 
  // we need a higher ion/electron macroparticle ratio.  

  species_t *ion_H, *ion_He; 
  if ( mobile_ions ) {
    if ( H_present  ) ion_H  = define_species("H",  Z_H,  mime_H,  max_local_np, max_local_nm, ion_sort_int, 1); 
    if ( He_present ) ion_He = define_species("He", Z_He, mime_He, max_local_np, max_local_nm, ion_sort_int, 1); 
  }

  // Set up grid
  sim_log("Setting up computational grid."); 
  define_units( 1, 1 );
  define_timestep( dt );

  sim_log("Setting up periodic mesh."); 
  define_absorbing_grid( 0,         -0.5*Ly,    -0.5*Lz,        // Low corner
                         Lx,         0.5*Ly,     0.5*Lz,        // High corner 
                         nx,         ny,         nz,            // Resolution
                         topology_x, topology_y, topology_z,    // Topology
                         reflect_particles );                   // Default particle boundary condition 

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

  sim_log("Overriding x boundaries to absorb fields."); 
  int ix, iy, iz;        // Domain location in mesh
  RANK_TO_INDEX( int(rank()), ix, iy, iz ); 

  // Set up Maxwellian reinjection B.C. 

  sim_log("Setting up Maxwellian reinjection boundary condition."); 

  particle_bc_t * maxwellian_reinjection =
    define_particle_bc( maxwellian_reflux( species_list, entropy ) );
  set_reflux_temp( maxwellian_reinjection, electron, uth_e, uth_e );
  if( mobile_ions ) {
    if( H_present  ) set_reflux_temp( maxwellian_reinjection, ion_H,  uthi_H,  uthi_H );
    if( He_present ) set_reflux_temp( maxwellian_reinjection, ion_He, uthi_He, uthi_He );
  }

  // Set up materials
  sim_log("Setting up materials."); 
  define_material( "vacuum", 1 );
  define_field_array( NULL, damp ); 
 
  // Paint the simulation volume with materials and boundary conditions
# define iv_region (   x<      hx*iv_thick || x>Lx  -hx*iv_thick  \
                    || y<-Ly/2+hy*iv_thick || y>Ly/2-hy*iv_thick  \
                    || z<-Lz/2+hz*iv_thick || z>Lz/2-hz*iv_thick ) /* all boundaries are i.v. */ 
  set_region_bc( iv_region, maxwellian_reinjection, maxwellian_reinjection, maxwellian_reinjection ); 

  // Load particles 
  if ( load_particles ) {
    sim_log("Loading particles.");
    // Fast load of particles--don't bother fixing artificial domain correlations
    double xmin=grid->x0, xmax=grid->x1;
    double ymin=grid->y0, ymax=grid->y1;
    double zmin=grid->z0, zmax=grid->z1;
    repeat( N_e/(topology_x*topology_y*topology_z) ) {
      double x = uniform( rng(0), xmin, xmax );
      double y = uniform( rng(0), ymin, ymax );
      double z = uniform( rng(0), zmin, zmax );
      if ( iv_region ) continue;           // Particle fell in iv_region.  Don't load.
      inject_particle( electron, x, y, z,
                       normal( rng(0), 0, uth_e ),
                       normal( rng(0), 0, uth_e ),
                       normal( rng(0), 0, uth_e ), q_e, 0, 0 );
      if ( mobile_ions ) {
        if ( H_present )  // Inject an H macroion on top of macroelectron
          inject_particle( ion_H, x, y, z, 
                           normal( rng(0), 0, uthi_H ), 
                           normal( rng(0), 0, uthi_H ), 
                           normal( rng(0), 0, uthi_H ), qi_H, 0, 0 ); 
        if ( He_present ) // Inject an H macroion on top of macroelectron
          inject_particle( ion_He, x, y, z, 
                           normal( rng(0), 0, uthi_He ),
                           normal( rng(0), 0, uthi_He ),
                           normal( rng(0), 0, uthi_He ), qi_He, 0, 0 ); 
      }
    }
  }

  sim_log("***Finished with user-specified initialization ***"); 

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
  if ( step()%200==0 ) sim_log("Time step: "<<step()); 

# define should_dump(x) \
  (global->x##_interval>0 && remainder(step(),global->x##_interval)==0)

  // Turn off these diagnostics for GB runs 
# if 0

 if ( step()==0 ) {
    // A grid dump contains all grid parameters, field boundary conditions,
    // particle boundary conditions and domain connectivity information. This
    // is stored in a binary format. Each rank makes a grid dump
    dump_grid("rundata/grid");

    // A materials dump contains all the materials parameters. This is in a
    // text format. Only rank 0 makes the materials dump
    dump_materials("rundata/materials");

    // A species dump contains the physics parameters of a species. This is in
    // a text format. Only rank 0 makes the species dump
    dump_species("rundata/species");
  }

  // Field and hydro data

  if ( should_dump(field) ) {
    dump_fields("field/fields");
    if ( global->load_particles ) {
      dump_hydro( "electron", "hydro/e_hydro" );
      if ( global->mobile_ions ) {
        if ( global->H_present  ) dump_hydro( "H",  "hydro/H_hydro"  );
        if ( global->He_present ) dump_hydro( "He", "hydro/He_hydro" );
      }

  }

  // Partcile dump data

  if ( should_dump(particle) && global->load_particles ) {
    dump_particles( "electron", "particle/eparticle" );
    if ( global->mobile_ions ) {
      if ( global->H_present  ) dump_particles( "H",  "particle/Hparticle" );
      if ( global->He_present ) dump_particles( "He", "particle/Heparticle" );
    }
  } 

# endif 

  // Ponyting data - send to stdout for GB run
  // Write Poynting flux at left boundary
  // Poynting flux is defined positive if directed in the +x direction. 
  // TODO: It is assumed that Ponyting dumps are infrequent, so we can afford the 
  // mpi_allreduce() here.  This needs to be verified on RR hardware. 
  
#define ALLOCATE(A,LEN,TYPE)                                             \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate.")); 

  static float *pvec=NULL; 
  static double psum, gpsum; 
  static FILE *fp_poynting;
  static char fname_poynting[]="poynting/poynting";
  static int stride, initted=0; 
  
  if ( !initted ) { 
    stride=(grid->ny-1)*(grid->nz-1); 
    ALLOCATE( pvec, stride, float ); 
    initted=1; 
  }

  if ( step()>0 && should_dump(poynting) ) {
    int i, j, k, k1, k2, ix, iy, iz; 
    for ( i=0; i<stride; ++i ) pvec[i]=0;     // Initialize pvec to zero. 
    // FIXME: THIS COULD BE DONE BETTER WITH VPIC.HXX HELPERS
    RANK_TO_INDEX( int(rank()), ix, iy, iz ); 
    if ( ix==0 ) {                            // Compute Poynting for domains on left of box
      for ( j=1; j<grid->ny; ++j ) {
        for ( k=1; k<grid->nz; ++k ) {
          k1 = INDEX_FORTRAN_3(1,j+1,k+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
          k2 = INDEX_FORTRAN_3(2,j+1,k+1,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1);
          pvec[(j-1)*(grid->nz-1)+k-1] = (  field(k2).ey*0.5*(field(k1).cbz+field(k2).cbz)
                                          - field(k2).ez*0.5*(field(k1).cby+field(k2).cby) )
                                          / (grid->cvac*grid->cvac*global->e0*global->e0);
        }
      }
    }                                         // Leave pvec = zero in mp_allsum_d for interior

    // Sum poynting flux on surface
    for ( i=0, psum=0; i<stride; ++i ) psum+=pvec[i]; 
    // Sum over all surfaces
    mp_allsum_d(&psum, &gpsum, 1);
    // Divide by number of mesh points summed over
    gpsum /= stride*global->topology_y*global->topology_z; 
    
    if ( rank()==0 && global->write_poynting_data ) { 
      fp_poynting=fopen( fname_poynting, (step()==global->poynting_interval ? "w" : "rb+") );
      if ( !fp_poynting ) ERROR(("Could not open file."));
      fseek( fp_poynting,
             (step()/global->poynting_interval-1)*sizeof(double), SEEK_SET );
      fwrite( &gpsum, sizeof(double), 1, fp_poynting );
      fclose(fp_poynting);
    } 
    sim_log("** step = "<<step()<<" Poynting = "<<gpsum);  // Dump data to stdout
  }

# if 0
  // Restart dump 

  if ( should_dump(restart) ) {
    char *restart_fbase[] = { "restart/restart0", "restart/restart1" }; 
    dump_restart( restart_fbase[global->rtoggle], 0 ); 
    global->rtoggle^=1; 
  } 

  if ( step()>0 && global->quota_check_interval && (step()%global->quota_check_interval)==0 ) { 
    if ( mp_elapsed(grid->mp) > global->quota_sec ) {
      dump_restart( "restart/restart", 0 ); 
      sim_log( "Restart dump restart completed." ); 
      sim_log( "Allowed runtime exceeded for this job.  Terminating." ); 
      mp_barrier( grid->mp ); // Just to be safe
      mp_finalize( grid->mp ); 
      exit(0); 
    }
  }
# endif 

} 


begin_field_injection { 
  // Inject a light wave from lhs boundary with E aligned along y
  // Use scalar diffraction theory for the Gaussian beam source.  (This is approximate). 

  // For quiet startup (i.e., so that we don't propagate a delta-function noise
  // pulse at time t=0) we multiply by a constant phase term exp(i phi) where: 
  //   phi = k*global->xfocus+atan(h)    (3d) 

  // Inject from the left a field of the form ey = e0 sin( omega t )

# define DY    ( grid->y0 + (iy-0.5)*grid->dy - global->ycenter )
# define DZ    ( grid->z0 + (iz-1  )*grid->dz - global->zcenter )
# define R2    ( DY*DY + DZ*DZ )                                   
# define PHASE ( global->omega*t + h*R2/(global->width*global->width) )
# define MASK  ( R2<=pow(global->mask*global->width,2) ? 1 : 0 )

  if ( grid->x0==0 ) {               // Node is on left boundary
    double alpha      = grid->cvac*grid->dt/grid->dx;
    double emax_coeff = (4/(1+alpha))*global->omega*grid->dt*global->e0;
    double prefactor  = emax_coeff*sqrt(2/M_PI); 
    double t          = grid->dt*step(); 

    // Compute Rayleigh length in c/wpe
    double rl         = M_PI*global->waist*global->waist/global->lambda; 

    double pulse_shape_factor = 1;
    float pulse_length        = 70;  // units of 1/wpe
    float sin_t_tau           = sin(0.5*t*M_PI/pulse_length);
    pulse_shape_factor        = ( t<pulse_length ? sin_t_tau : 1 );
    double h                  = global->xfocus/rl;   // Distance / Rayleigh length

    // Loop over all Ey values on left edge of this node
    for ( int iz=1; iz<=grid->nz+1; ++iz ) 
      for ( int iy=1; iy<=grid->ny; ++iy )  
        field(1,iy,iz).ey += prefactor*cos(PHASE)*exp(-R2/(global->width*global->width))*MASK*pulse_shape_factor; 
  }
}


begin_particle_injection {
  // No particle injection for this simulation
}


begin_current_injection {
  // No current injection for this simulation
}


begin_particle_collisions {
  // No collisions for this simulation
}


