//////////////////////////////////////////////////////
//
//   Harris Sheet Reconnection - Open Boundary Model
//
//////////////////////////////////////////////////////

#include "injection"   //  Subroutine to compute re-injection velocity

//////////////////////////////////////////////////////

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
  double v_A;
  double topology_x;       // domain topology
  double topology_y;
  double topology_z;

  // parameters for the collision model

  int    ee_collisions;         // Flag to signal we want to do e-e collisions.
  int    ei_collisions;         // Flag to signal we want to do e-i collisions.
  int    ii_collisions;         // Flag to signal we want to do i-i collisions.

  double cvar;                  // Base variance (dimensionless) used in particle collision
  double nppc_max;              // Max number of particles/cell (used to define key array size).
  int tstep_coll;            // Collision interval (=multiple of sort interval).
  double Z;
  double mi_me;
  double wpewce;

  //  Variables for Open BC Model
  double nb;        // Background density
  int nsp;          //  Number of Species
  double vth[2];    // Thermal velocity of Harris components
  double vthb[2];    // Thermal velocity of background components
  double q[2];      // Species charge
  double L_de;      //  Initial Harris sheet thickness
  double uf[2];      // Initial Fluid Drift in Harris
  double nfac;      //  Normalization factor to convert particles per cell into density
  double rin[3];     //  Relaxation parameter for inflow boundary moments
  double rout[3];   //  Relaxation parameter for outlfow boundary moments
  double sort[2];   // Intervals where we know particles are sorted
  double edrive;    // Drive field for inflow boundary
  double tdrive;
  int left,right,top,bottom;  // Keep track of boundary domains
  double *nbot, *ubot, *pbot, *bbot, *fbot;         // Moments for bottom injectors
  double *ntop, *utop, *ptop, *btop, *ftop;         // Moments for top injectors
  double *nleft, *uleft, *pleft, *bleft, *fleft;     // Moments for left injectors
  double *nright, *uright, *pright, *bright, *fright; // Moments for right injectors

//  Variables for new output format

  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;
};

begin_initialization {
 // use natural PIC units

 double ec   = 1;         // Charge normalization
 double me   = 1;         // Mass normalization
 double c    = 1;         // Speed of light
 double de   = 1;         // Length normalization (electron inertial length)
 double eps0 = 1;         // Permittivity of space

  double cfl_req   = 0.99;  // How close to Courant should we try to run
  double wpedt_max = 0.36;  // How big a timestep is allowed if Courant is not too restrictive
  double damp      = 0.0;   // Level of radiation damping
  int rng_seed     = 1;     // Random number seed increment

  // Physics parameters

  double mi_me   = 400.0;    // Ion mass / electron mass
  double L_di    = 0.5;    // Sheet thickness / ion inertial length
  double Ti_Te   = 5.0;    // Ion temperature / electron temperature
  double Z   = 1.0;      // Ion charge
  double nb_n0   = 0.228;   // background plasma density
  double Tbe_Te  = 0.7598;  // Ratio of background electron temperature to Harris electron temperature
  double Tbi_Ti  = 0.3039;  // Ratio of background ion temperature to Harris ion temperature
  double wpe_wce = 2.0;      // electron plasma freq / electron cyclotron freq
  double bg = 0.0;        // electron plasma freq / electron cyclotron freq
  double theta   = 0;      // B0 = Bx
  double taui    = 300;    // simulation wci's to run

  double quota   = 11.0;   // run quota in hours
  double quota_sec = quota*3600;  // Run quota in seconds

  double cs   = cos(theta);
  double sn   = sin(theta);

  //derived qunatities
  double mi = me*mi_me;                                   // Ion mass
  double Te = me*c*c/(2*eps0*wpe_wce*wpe_wce*(1+Ti_Te));  // Electron temperature
  double Ti = Te*Ti_Te;                                   // Ion temperature
  double vthe = sqrt(Te/me);                              // Electron thermal velocity
  double vthi = sqrt(Ti/mi);                              // Ion thermal velocity
  double vtheb = sqrt(Tbe_Te*Te/me);             // normalized background electron thermal velocity
  double vthib = sqrt(Tbi_Ti*Ti/mi);             // normalized background ion thermal velocity
  double wci  = 1.0/(mi_me*wpe_wce);                      // Ion cyclotron frequency
  double wce  = wci*mi_me;                                // Electron cyclotron freqeuncy
  double wpe  = wce*wpe_wce;                              // electron plasma frequency
  double wpi  = wpe/sqrt(mi_me);                          // ion plasma frequency
  double di   = c/wpi;                                    // ion inertial length
  double L    = L_di*di;                                  // Harris sheet thickness
  double rhoi_L = sqrt(Ti_Te/(1.0+Ti_Te))/L_di;
  double v_A= (wci/wpi)/sqrt(nb_n0); //  based on nb

  //  Parameters for Open BC model

  double rin[3] =  {0.000, 0.06, 0.000};  //  Relaxation - density, velocity + particle flux, pressure tensor
  double rout[3] = {0.002, 0.002, 0.002};  //  Relaxation - density, velocity + particle flux, pressure tensor
  double edrive = 0.08*v_A/(wpe_wce);    //  Setting edrive = 0 will give undriven limit
  //  double edrive = 0.0099;    //  Setting edrive = 0 will give undriven limit

  double tdrive = 32000.0;
  double sort_interval = 10;        //  Injector moments are also updated at this internal

  // Numerical parameters

  double nppc          = 400; // Average number of macro particle per cell per species

  double Lx            = 20*di; // size of box in x dimension
  double Ly            = 0.00781*di;  // size of box in y dimension
  double Lz            = 20*di; // size of box in z dimension

  int topology_factor = 64;

  double topology_x = 256/topology_factor;  // Number of domains in x, y, and z
  double topology_y = 1;
  double topology_z = 1;  // For load balance, best to keep "1" or "2" for Harris sheet

  double nx = 2560/sqrt(topology_factor);
  double ny = 1;
  double nz = 2560/sqrt(topology_factor);

  double hx = Lx/nx;
  double hy = Ly/ny;
  double hz = Lz/nz;

  double b0            = me*c*wce/ec; // Asymptotic magnetic field strength
  double n0            = me*eps0*wpe*wpe/(ec*ec);  // Peak electron (ion) density
  double vdri          = 2*c*Ti/(ec*b0*L);   // Ion drift velocity
  double vdre          = -vdri/(Ti_Te);           // electron drift velocity

  double Npe_sheet     = 2*n0*Lx*Ly*L*tanh(0.5*Lz/L); // number of physical electrons in sheet
  double Npe_back      = nb_n0*n0*Ly*Lz*Lx;    // Number of physical electrons in background
  double Npe           = Npe_sheet + Npe_back;
  double Ne            = nppc*nx*ny*nz;  // total macro electrons in box
  double Ne_sheet      = Ne*Npe_sheet/Npe;
  double Ne_back       = Ne*Npe_back/Npe;
  Ne_sheet  = trunc_granular(Ne_sheet,nproc()); // Make it divisible by number of processors
  Ne_back  = trunc_granular(Ne_back,nproc()); // Make it divisible by number of processors
  Ne = Ne_sheet + Ne_back;
  double qe_s = -ec*Npe_sheet/Ne_sheet;  // Charge per macro electron
  double qe_b = -ec*Npe_back/Ne_back;  // Charge per macro electron
  double qi_s =  ec*Npe_sheet/Ne_sheet;  // Charge per macro electron
  double qi_b =  ec*Npe_back/Ne_back;  // Charge per macro electron
  double nfac = qi_s/(hx*hy*hz);            // Convert density to particles per cell

  double gdri = 1/sqrt(1-vdri*vdri/(c*c));  // gamma of ion drift frame
  double gdre = 1/sqrt(1-vdre*vdre/(c*c)); // gamma of electron drift frame
  double udri = vdri*gdri;                 // 4-velocity of ion drift frame
  double udre = vdre*gdre;                  // 4-velocity of electron drift frame
  double tanhf = tanh(0.5*Lz/L);
  double Lpert = 1.5*Lx;   // wavelength of perturbation
  double dbz =  0.03*b0; //  Perturbation in Bz relative to Bo (Only change here)
  double dbx = -dbz*Lpert/(2*Lz); // Set Bx perturbation so that div(B) = 0

  // Determine the time step

  double dg = courant_length(Lx,Ly,Lz,nx,ny,nz);        // courant length
  double dt = cfl_req*dg/c;                             // courant limited time step
  if( wpe*dt>wpedt_max) dt=wpedt_max/wpe;               // override timestep if plasma frequency limited

  // Intervals for output

  int restart_interval = 8000;
  int energies_interval = 50;
  int interval = int(0.25/(wci*dt));
  int fields_interval = interval;
  int ehydro_interval = interval;
  int Hhydro_interval = interval;
  int eparticle_interval = 8*interval;
  int Hparticle_interval = 8*interval;
  int quota_check_interval     = 100;

  // Collision parameters
  // In CGS, variance of tan theta = 2 pi e^4 n_e dt_coll loglambda / (m_ab^2 u^3)

  int ii_collisions = 1;                   // flags to do collisions between the corresponding species
  int ee_collisions = 1;                   //
  int ei_collisions = 1;                   //

  int tstep_coll = (int) sort_interval;          // How frequently to apply collision operator
  double dt_coll = dt*(tstep_coll);          // in (1/wpe)
  double nuei_wce = 0.05;
  double cvar = dt_coll*3.0*sqrt(2.0*M_PI)/4.0*pow(vthe,3)*nuei_wce/wpe_wce;
  double nppc_max   = 20*nppc;             // Max possible number of particles/cell of ea. species
                                           // (to size the key array in collision handler)

  //  Determine which domains area along the boundaries - Use macro from grid/partition.c

# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {                   \
    int _ix, _iy, _iz;                                                    \
    _ix  = (rank);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(topology_x);   /* iy = iy+gpy*iz */                    \
    _ix -= _iy*int(topology_x);   /* ix = ix */                           \
    _iz  = _iy/int(topology_y);   /* iz = iz */                           \
    _iy -= _iz*int(topology_y);   /* iy = iy */                           \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \
  } END_PRIMITIVE

  int ix, iy, iz, left=0,right=0,top=0,bottom=0;
  RANK_TO_INDEX( int(rank()), ix, iy, iz );
  if ( ix ==0 ) left=1;
  if ( ix ==topology_x-1 ) right=1;
  if ( iz ==0 ) bottom=1;
  if ( iz ==topology_z-1 ) top=1;

  ///////////////////////////////////////////////
  // Setup high level simulation parameters
  num_step             = int(taui/(wci*dt));
  num_step             = 1000;
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
  global->v_A  = v_A;
  global->edrive  = edrive;
  global->tdrive  = tdrive;

  global->topology_x  = topology_x;
  global->topology_y  = topology_y;
  global->topology_z  = topology_z;

  global->left = left;
  global->right = right;
  global->top = top;
  global->bottom = bottom;

  //  Parameters for the open boundary model

  global->nsp = 2;
  global->nb  = nb_n0;
  global->rin[0]  = rin[0];
  global->rin[1]  = rin[1];
  global->rin[2]  = rin[2];
  global->rout[0]  = rout[0];
  global->rout[1]  = rout[1];
  global->rout[2]  = rout[2];
  global->vth[0]  = sqrt(2)*vthe;
  global->vth[1]  = sqrt(2)*vthi;
  global->vthb[0]  = sqrt(2)*vtheb;
  global->vthb[1]  = sqrt(2)*vthib;
  global->q[0]  = qe_s;
  global->q[1]  = qi_s;
  global->uf[0]  = udre;
  global->uf[1]  = udri;
  global->sort[0]  = sort_interval;
  global->sort[1]  = sort_interval;
  global->nfac  = nfac;
  global->L_de  = L;

  // Collision model parameters

  global->ee_collisions            = ee_collisions;
  global->ii_collisions            = ii_collisions;
  global->ei_collisions            = ei_collisions;

  global->cvar                     = cvar;
  global->nppc_max                 = nppc_max;
  global->tstep_coll               = tstep_coll;
  global->mi_me                    = mi_me;
  global->Z                        = Z;
  global->wpewce                   = wpe_wce;
  global->nfac                     = nfac;

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the grid

  // Setup basic grid parameters
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = c;
  grid->eps0 = eps0;
  grid->damp = damp;

  // Define the grid

  define_periodic_grid(  0, -0.5*Ly, -0.5*Lz,    // Low corner
                          Lx, 0.5*Ly, 0.5*Lz,     // High corner
                          nx, ny, nz,             // Resolution
                          topology_x, topology_y, topology_z); // Topology

 // ***** Set Field Boundary Conditions *****

  //  sim_log("Absorbing fields on X & Z-boundaries");
  //if ( iz==0 )            set_domain_field_bc( BOUNDARY(0,0,-1), absorb_fields );
  //if ( iz==topology_z-1 ) set_domain_field_bc( BOUNDARY( 0,0,1), absorb_fields );
  if ( ix==0 )            set_domain_field_bc( BOUNDARY(-1,0,0), absorb_fields );
  if ( ix==topology_x-1 ) set_domain_field_bc( BOUNDARY( 1,0,0), absorb_fields );

  sim_log("Conducting fields on X & Z-boundaries");
  if ( iz==0 )            set_domain_field_bc( BOUNDARY(0,0,-1), pec_fields );
  if ( iz==topology_z-1 ) set_domain_field_bc( BOUNDARY( 0,0,1), pec_fields );
  //if ( ix==0 )            set_domain_field_bc( BOUNDARY(-1,0,0), pec_fields );
  //if ( ix==topology_x-1 ) set_domain_field_bc( BOUNDARY( 1,0,0), pec_fields );

 // ***** Set Particle Boundary Conditions *****

  sim_log("Absorb particles on X & Z-boundaries");
  if ( iz==0 )            set_domain_particle_bc( BOUNDARY(0,0,-1), absorb_particles );
  if ( iz==topology_z-1 ) set_domain_particle_bc( BOUNDARY(0,0,1), absorb_particles );
  if ( ix==0 )            set_domain_particle_bc( BOUNDARY(-1,0,0), absorb_particles );
  if ( ix==topology_x-1 ) set_domain_particle_bc( BOUNDARY(1,0,0), absorb_particles );

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the species

  sim_log("Setting up species. ");
  species_t *electron = define_species("electron",-ec/me,2.5*Ne/nproc(),-1,sort_interval,0);
  species_t *ion = define_species("ion", ec/mi,2.5*Ne/nproc(),-1,sort_interval,0);

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup materials

  sim_log("Setting up materials. ");

  define_material( "vacuum", 1 );
  define_material( "resistive",1,1,1);

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.

////////////////////////////////////////////////////////////////////////////////////////////
//  Finalize Field Advance

  sim_log("Finalizing Field Advance");

  finalize_field_advance(standard_field_advance);

  //  Define resistive layer surrounding boundary --> set thickness=0 to eliminate this feature

    double thickness = 2;
#define resistive_layer (x < hx*thickness || x > Lx-hx*thickness || z <-Lz/2+hz*thickness  || z > Lz/2-hz*thickness )

    if (thickness > 0) set_region_material(resistive_layer, "resistive", "resistive")    ;

  ///////////////////////////////////////////////////
  // Log diagnostic information about this simulation

  sim_log( "***********************************************" );
  sim_log("* Topology:                       "<<topology_x<<" "<<topology_y<<" "<<topology_z);
  sim_log ( "tanhf = " << tanhf );
  sim_log ( "L_di   = " << L_di );
  sim_log ( "rhoi/L   = " << rhoi_L );
  sim_log ( "Ti/Te = " << Ti_Te ) ;
  sim_log ( "nb/n0 = " << nb_n0 ) ;
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
  sim_log ( "damp = " << damp );
  sim_log ( "courant = " << c*dt/dg );
  sim_log ( "nproc = " << nproc ()  );
  sim_log ( "nppc = " << nppc );
  sim_log ( " b0 = " << b0 );
  sim_log ( " v_A (based on nb) = " << v_A );
  sim_log ( " di = " << di );
  sim_log ( " Ne = " << Ne );
  sim_log ( " Ne_sheet = " << Ne_sheet );
  sim_log ( " Ne_back = " << Ne_back );
  sim_log ( "total # of particles = " << 2*Ne );
  sim_log ( "dt*wpe = " << wpe*dt );
  sim_log ( "dt*wce = " << wce*dt );
  sim_log ( "dt*wci = " << wci*dt );
  sim_log ( " energies_interval: " << energies_interval );
  sim_log ( "dx/de = " << Lx/(de*nx) );
  sim_log ( "dy/de = " << Ly/(de*ny) );
  sim_log ( "dz/de = " << Lz/(de*nz) );
  sim_log ( "dx/rhoi = " << (Lx/nx)/(vthi/wci)  );
  sim_log ( "dx/rhoe = " << (Lx/nx)/(vthe/wce)  );
  sim_log ( "L/debye = " << L/(vthe/wpe)  );
  sim_log ( "dx/debye = " << (Lx/nx)/(vthe/wpe)  );
  sim_log ( "n0 = " << n0 );
  sim_log ( "vthi/c = " << vthi/c );
  sim_log ( "vthe/c = " << vthe/c );
  sim_log ( "vdri/c = " << vdri/c );
  sim_log ( "vdre/c = " << vdre/c );
  sim_log("* nu/wce:                        "<<nuei_wce);
  sim_log("* nu*dt_coll:                    "<<nuei_wce/wpe_wce*dt_coll);

  // Dump simulation information to file "info"
  if (rank() == 0 ) {
    FILE *fp_info;
    if ( ! (fp_info=fopen("info", "w")) ) ERROR(("Cannot open file."));
    fprintf(fp_info, "           ***** Simulation parameters ***** \n");
    fprintf(fp_info, "           L/di   =               %e\n", L_di);
    fprintf(fp_info, "           L/de   =               %e\n", L/de);
    fprintf(fp_info, "           rhoi/L =               %e\n", rhoi_L);
    fprintf(fp_info, "           Ti/Te  =               %e\n", Ti_Te );
    fprintf(fp_info, "           Tbi/Ti =               %e\n", Tbi_Ti );
    fprintf(fp_info, "           Tbe/Te =               %e\n", Tbe_Te );
    fprintf(fp_info, "           nb/n0 =                %e\n", nb_n0 );
    fprintf(fp_info, "           wpe/wce =              %e\n", wpe_wce );
    fprintf(fp_info, "           mi/me =                %e\n", mi_me );
    fprintf(fp_info, "           theta =                %e\n", theta );
    fprintf(fp_info, "           taui =                 %e\n", taui );
    fprintf(fp_info, "           num_step =             %i\n", num_step );
    fprintf(fp_info, "           Lx/de =                %e\n", Lx/de );
    fprintf(fp_info, "           Ly/de =                %e\n", Ly/de );
    fprintf(fp_info, "           Lz/de =                %e\n", Lz/de );
    fprintf(fp_info, "           Lx/di =                %e\n", Lx/di );
    fprintf(fp_info, "           Ly/di =                %e\n", Ly/di );
    fprintf(fp_info, "           Lz/di =                %e\n", Lz/di );
    fprintf(fp_info, "           nx =                   %e\n", nx );
    fprintf(fp_info, "           ny =                   %e\n", ny );
    fprintf(fp_info, "           nz =                   %e\n", nz );
    fprintf(fp_info, "           damp =                 %e\n", damp );
    fprintf(fp_info, "           courant =              %e\n", c*dt/dg );
    fprintf(fp_info, "           nproc =                %e\n", nproc() );
    fprintf(fp_info, "           nppc =                 %e\n", nppc );
    fprintf(fp_info, "           b0 =                   %e\n", b0 );
    fprintf(fp_info, "           v_A (based on nb) =    %e\n", v_A );
    fprintf(fp_info, "           di =                   %e\n", di );
    fprintf(fp_info, "           Ne =                   %e\n", Ne );
    fprintf(fp_info, "           Ne_sheet =             %e\n", Ne_sheet );
    fprintf(fp_info, "           Ne_back =              %e\n", Ne_back );
    fprintf(fp_info, "           total # of particles = %e\n", 2*Ne );
    fprintf(fp_info, "           dt*wpe =               %e\n", wpe*dt );
    fprintf(fp_info, "           dt*wce =               %e\n", wce*dt );
    fprintf(fp_info, "           dt*wci =               %e\n", wci*dt );
    fprintf(fp_info, "           energies_interval:     %i\n", energies_interval);
    fprintf(fp_info, "           dx/de =                %e\n", Lx/(de*nx) );
    fprintf(fp_info, "           dy/de =                %e\n", Ly/(de*ny) );
    fprintf(fp_info, "           dz/de =                %e\n", Lz/(de*nz) );
    fprintf(fp_info, "           L/debye =              %e\n", L/(vthe/wpe) );
    fprintf(fp_info, "           dx/rhoi =              %e\n", (Lx/nx)/(vthi/wci) );
    fprintf(fp_info, "           dx/rhoe =              %e\n", (Lx/nx)/(vthe/wce) );
    fprintf(fp_info, "           dx/debye =             %e\n", (Lx/nx)/(vthe/wpe) );
    fprintf(fp_info, "           n0 =                   %e\n", n0 );
    fprintf(fp_info, "           vthi/c =               %e\n", vthi/c );
    fprintf(fp_info, "           vthe/c =               %e\n", vthe/c );
    fprintf(fp_info, "           vdri/c =               %e\n", vdri/c );
    fprintf(fp_info, "           vdre/c =               %e\n", vdre/c );
    fprintf(fp_info, " tstep_coll:                    %i\n", tstep_coll);
    fprintf(fp_info, " nu/wce:                        %g\n", nuei_wce);
    fprintf(fp_info, " nu*dt_coll:                    %g\n", nuei_wce/wpe_wce*dt_coll);
    fprintf(fp_info, "           ***************************\n");
    fclose(fp_info);
}


  ////////////////////////////
  // Load fields

  sim_log( "Loading fields" );
  set_region_field( everywhere, 0, 0, 0,                    // Electric field
                    cs*b0*tanh(z/L)+dbx*cos(2.0*M_PI*(x-0.5*Lx)/Lpert)*sin(M_PI*z/Lz), //Bx
                    -sn*b0*tanh(z/L) + b0*bg, //By
                    dbz*cos(M_PI*z/Lz)*sin(2.0*M_PI*(x-0.5*Lx)/Lpert) ); // Bz

//  Localized Perturbation to lauch a light wave

//   # define R2 ((x-0.5*Lx)*(x-0.5*Lx) + z*z)/(L*L)
//     # define PERT  0.2*tanh(R2)/cosh(R2)
//     set_region_field( everywhere, 0, 0, 0,                    // Electric field
//                  cs*b0*tanh(z/L) - z*PERT, //Bx
//                  b0*bg, //By
//                  (x-0.5*Lx)*PERT ); // Bz

  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specied as logical equations (i.e. x>0 && x+y<2)

  // LOAD PARTICLES

  sim_log( "Loading particles" );

  // Do a fast load of the particles

  seed_rand( rng_seed*nproc() + rank() );  //Generators desynchronized
  double xmin = grid->x0 , xmax = grid->x0+(grid->dx)*(grid->nx);
  double ymin = grid->y0 , ymax = grid->y0+(grid->dy)*(grid->ny);
  double zmin = grid->z0 , zmax = grid->z0+(grid->dz)*(grid->nz);

  // Load Harris population

  sim_log( "-> Main Harris Sheet" );

  repeat ( Ne_sheet/nproc() ) {
    double x, y, z, ux, uy, uz, d0 ;

    do {
      z = L*atanh(uniform_rand(-1,1)*tanhf);
    } while( z<= zmin || z>=zmax );
    x = uniform_rand(xmin,xmax);
    y = uniform_rand(ymin,ymax);

    // inject_particles() will return an error for particles no on this
    // node and will not inject particle locally

    ux = maxwellian_rand(vthe);
    uy = maxwellian_rand(vthe);
    uz = maxwellian_rand(vthe);
    d0 = gdre*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udre;
    uy = d0*cs - ux*sn;
    ux = d0*sn + ux*cs;

    inject_particle(electron, x, y, z, ux, uy, uz, qe_s, 0, 0 );

    ux = maxwellian_rand(vthi);
    uy = maxwellian_rand(vthi);
    uz = maxwellian_rand(vthi);
    d0 = gdri*uy + sqrt(ux*ux + uy*uy + uz*uz + 1)*udri;
    uy = d0*cs - ux*sn;
    ux = d0*sn + ux*cs;

    inject_particle(ion, x, y, z, ux, uy, uz, qi_s, 0, 0 );

  }

  sim_log( "-> Background Population" );

  repeat ( Ne_back/nproc() ) {

  double z = uniform_rand(zmin,zmax);
  double x = uniform_rand(xmin,xmax);
  double y = uniform_rand(ymin,ymax);

  inject_particle( electron, x, y, z,
                   maxwellian_rand(vtheb),
                   maxwellian_rand(vtheb),
                   maxwellian_rand(vtheb),qe_b, 0, 0);

  inject_particle( ion, x, y, z,
                   maxwellian_rand(vthib),
                   maxwellian_rand(vthib),
                   maxwellian_rand(vthib),qi_b, 0 ,0 );
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

        //      global->fdParams.output_variables( electric | magnetic );
        //global->hedParams.output_variables( current_density | charge_density );
        //global->hHdParams.output_variables( current_density | charge_density );

        global->fdParams.output_variables( all );
        global->hedParams.output_variables( all );
        global->hHdParams.output_variables( all );

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
        (global->x##_interval>0 && remainder(step, global->x##_interval) == 0)

#include <FileIO.hxx>

begin_diagnostics {

  const int nsp=global->nsp;
  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;

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
        if(step==0) {
                dump_mkdir("fields");
                dump_mkdir("hydro");
                dump_mkdir("rundata");
                dump_mkdir("injectors");
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
                dump_energies("rundata/energies", step == 0 ? 0 : 1);
        } // if

        /*--------------------------------------------------------------------------
         * Field data output
         *------------------------------------------------------------------------*/

        if(step == 1 || should_dump(fields)) field_dump(global->fdParams);

        /*--------------------------------------------------------------------------
         * Electron species output
         *------------------------------------------------------------------------*/

        if(should_dump(ehydro)) hydro_dump("electron", global->hedParams);

        /*--------------------------------------------------------------------------
         * Ion species output
         *------------------------------------------------------------------------*/

        if(should_dump(Hhydro)) hydro_dump("ion", global->hHdParams);

        /*--------------------------------------------------------------------------
         * Restart dump
         *------------------------------------------------------------------------*/

        if(step && !(step%global->restart_interval)) {
                if(!global->rtoggle) {
                        global->rtoggle = 1;
                        dump_restart("restart1/restart", 0);
                        DUMP_INJECTORS(1);
                }
                else {
                        global->rtoggle = 0;
                        dump_restart("restart2/restart", 0);
                        DUMP_INJECTORS(2);
                } // if
        } // if

  // Dump particle data

        char subdir[36];
        if ( should_dump(eparticle) && step !=0 && step > 56*(global->fields_interval)  ) {
                  //    if ( should_dump(eparticle) && step !=0 ) {
          sprintf(subdir,"particle/T.%d",step);
          dump_mkdir(subdir);
          sprintf(subdir,"particle/T.%d/eparticle",step);
          dump_particles("electron",subdir);
        }


//   if ( should_dump(Hparticle) ) {
//     dump_particles("ion",  "Hparticle");
//   }

  // Shut down simulation when wall clock time exceeds global->quota_sec.
  // Note that the mp_elapsed() is guaranteed to return the same value for all
  // processors (i.e., elapsed time on proc #0), and therefore the abort will
  // be synchronized across processors. Note that this is only checked every
  // few timesteps to eliminate the expensive mp_elapsed call from every
  // timestep. mp_elapsed has an ALL_REDUCE in it!

    if( step>0 && global->quota_check_interval>0 && (step&global->quota_check_interval)==0 ) {
    if( mp_elapsed( grid->mp ) > global->quota_sec ) {
      sim_log( "Allowed runtime exceeded for this job.  Terminating....\n");
      dump_restart("restart0/restart",0);
      sim_log( "Restart dump restart completed." );
      DUMP_INJECTORS(0);
      mp_barrier( grid->mp ); // Just to be safe
      mp_finalize( grid->mp );
      exit(0); // Exit or abort?
    }
  }

} // end diagnostics

// *******************  PARTICLE INJECTION  - OPEN BOUNDARY ***************************

begin_particle_injection {
  int inject;
  double x, y, z, age, flux, vtherm, vd;
  double uv[3];
  double zcell;
  const int nsp=global->nsp;
  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;
  const double rin[3]={global->rin[0],global->rin[1],global->rin[2]};
  const double rout[3]={global->rout[0],global->rout[1],global->rout[2]};
  const double sqpi =1.772453850905516;
  const double dt=grid->dt;
  const double hx=grid->dx;
  const double hy=grid->dy;
  const double hz=grid->dz;
  const double nb=global->nb;
  const double nfac=global->nfac;

  // Initialize the injectors on the first call

    static int initted=0;
    if ( !initted ) {

      initted=1;

      if (rank() == 0) MESSAGE(("----------------Initializing the Particle Injectors-----------------"));

      // Intialize injectors for Harris Sheet with a uniform background

      if (global->right) {

        DEFINE_INJECTOR(right,ny,nz);

        if (step ==0) {
          for ( int n=1; n<=nsp; n++ ) {
            double cn = (uf(2)/vth(2))/(vth(n)/vth(2));
            for ( int k=1;k<=nz; k++ ) {
              for ( int j=1;j<=ny; j++ ) {
                bright(n,k,j) = 0;
                zcell = (grid->z0 + k*hz-hz/2)/(global->L_de);
                nright(n,k,j) = (nb + 1/(cosh(zcell)*cosh(zcell)))/nfac;
                fright(n,k,j) = (nb*vthb(n) + vth(n)/(cosh(zcell)*cosh(zcell)))/(2*hx*sqpi*nfac);
                uright(1,n,k,j) = 0;
                uright(2,n,k,j) = uf(n)/(1+nb*cosh(zcell)*cosh(zcell));
                uright(3,n,k,j) = 0;
                pright(1,2,n,k,j)=pright(2,1,n,k,j)=pright(1,3,n,k,j)=pright(3,1,n,k,j)=pright(2,3,n,k,j)=pright(3,2,n,k,j)=0;
                pright(1,1,n,k,j) = (nb*vthb(n)*vthb(n) + vth(n)*vth(n)/(cosh(zcell)*cosh(zcell)))/(2*nfac);
                pright(2,2,n,k,j) = (nb*vthb(n)*vthb(n) + vth(n)*vth(n)*(1/(cosh(zcell)*cosh(zcell))+2*nb*cn*cn/(1+nb*cosh(zcell)*cosh(zcell))))/(2*nfac);
                pright(3,3,n,k,j) = pright(1,1,n,k,j);
              }
            }
          }  // end for
        } // endif

        else
          READ_INJECTOR(right,ny,nz,0);

      } //end right boundary

      if (global->left) {

        DEFINE_INJECTOR(left,ny,nz);

        if (step==0) {
          for ( int n=1; n<=nsp; n++ ) {
            double cn = (uf(2)/vth(2))/(vth(n)/vth(2));
            for ( int k=1;k<nz+1; k++ ) {
              for ( int j=1;j<=ny; j++ ) {
                bleft(n,k,j) = 0;
                zcell = (grid->z0 + k*hz-hz/2)/(global->L_de);
                nleft(n,k,j) = (nb + 1/(cosh(zcell)*cosh(zcell)))/nfac;
                fleft(n,k,j) = (nb*vthb(n) + vth(n)/(cosh(zcell)*cosh(zcell)))/(2*hx*sqpi*nfac);
                uleft(1,n,k,j) = 0;
                uleft(2,n,k,j) = uf(n)/(1+nb*cosh(zcell)*cosh(zcell));
                uleft(3,n,k,j) = 0;
                pleft(1,2,n,k,j)=pleft(2,1,n,k,j)=pleft(1,3,n,k,j)=pleft(3,1,n,k,j)=pleft(2,3,n,k,j)=pleft(3,2,n,k,j)=0;
                pleft(1,1,n,k,j) = (nb*vthb(n)*vthb(n) + vth(n)*vth(n)/(cosh(zcell)*cosh(zcell)))/(2*nfac);
                pleft(2,2,n,k,j) = (nb*vthb(n)*vthb(n) + vth(n)*vth(n)*(1/(cosh(zcell)*cosh(zcell))+2*nb*cn*cn/(1+nb*cosh(zcell)*cosh(zcell))))/(2*nfac);
                pleft(3,3,n,k,j) = pleft(1,1,n,k,j);
              }
            }
          } // end for
        } //endif

        else
          READ_INJECTOR(left,ny,nz,0);

      } // end left boundary

      if (global->top) {

        DEFINE_INJECTOR(top,ny,nx);

        if (step==0) {
          for ( int n=1; n<=nsp; n++ ) {
            for ( int i=1;i<=nx; i++ ) {
              for ( int j=1;j<=ny; j++ ) {
                btop(n,i,j) = 0;
                ntop(n,i,j) = nb/nfac;
                ftop(n,i,j) = ntop(n,i,j)*vthb(n)/(2*hz*sqpi);
                utop(1,n,i,j) = 0;
                utop(2,n,i,j) = 0;
                utop(3,n,i,j) = 0;
                ptop(1,2,n,i,j)=ptop(2,1,n,i,j)=ptop(1,3,n,i,j)=ptop(3,1,n,i,j)=ptop(2,3,n,i,j)=ptop(3,2,n,i,j)=0;
                ptop(1,1,n,i,j) = ntop(n,i,j)*vthb(n)*vthb(n)/2;;
                ptop(2,2,n,i,j) = ntop(n,i,j)*vthb(n)*vthb(n)/2;;
                ptop(3,3,n,i,j) = ntop(n,i,j)*vthb(n)*vthb(n)/2;;
              }
            }
          } // end for
        } //endif

        else
          READ_INJECTOR(top,ny,nx,0);

      } // end top boundary

      if (global->bottom) {

        DEFINE_INJECTOR(bot,ny,nx);

        if (step ==0) {
          for ( int n=1; n<=nsp; n++ ) {
            for ( int i=1;i<=nx; i++ ) {
              for ( int j=1;j<=ny; j++ ) {
                bbot(n,i,j) = 0;
                nbot(n,i,j) = nb/nfac;
                fbot(n,i,j) = nbot(n,i,j)*vthb(n)/(2*hz*sqpi);
                ubot(1,n,i,j) = 0.0;
                ubot(2,n,i,j) = 0.0;
                ubot(3,n,i,j) = 0.0;
                pbot(1,2,n,i,j)=pbot(2,1,n,i,j)=pbot(1,3,n,i,j)=pbot(3,1,n,i,j)=pbot(2,3,n,i,j)=pbot(3,2,n,i,j)=0;
                pbot(1,1,n,i,j) = nbot(n,i,j)*vthb(n)*vthb(n)/2;
                pbot(2,2,n,i,j) = nbot(n,i,j)*vthb(n)*vthb(n)/2;
                pbot(3,3,n,i,j) = nbot(n,i,j)*vthb(n)*vthb(n)/2;
              }
            }
          } // end for
        } //endif

        else
          READ_INJECTOR(bot,ny,nx,0);

      }  // end bottom boundary

      if (rank() == 0) MESSAGE(("-------------------------------------------------------------------"));

    }// End of Intialization

    //  Inject particles on Left Boundary

    if (global->left) {
      for ( int n=1; n<=nsp; n++ ) {
        species_t * species = find_species_id(n-1,species_list );
        for ( int k=1;k<=nz; k++ ) {
          for ( int j=1;j<=ny; j++ ) {
            bleft(n,k,j) = bleft(n,k,j) + dt*fleft(n,k,j);
            inject = (int) bleft(n,k,j);
            bleft(n,k,j) = bleft(n,k,j) - (double) inject;
            double uflow[3] = {uleft(1,n,k,j),uleft(2,n,k,j),uleft(3,n,k,j)};
            double press[9] = {pleft(1,1,n,k,j),pleft(1,2,n,k,j),pleft(1,3,n,k,j),pleft(2,1,n,k,j),pleft(2,2,n,k,j),pleft(2,3,n,k,j),pleft(3,1,n,k,j),pleft(3,2,n,k,j),pleft(3,3,n,k,j)};
            repeat(inject) {
              compute_injection(uv,nleft(n,k,j),uflow,press,1,2,3,rng);
              x = grid->x0;
              y = grid->y0 + hy*(j-1) + hy*uniform_rand(0,1);
              z = grid->z0 + hz*(k-1) + hz*uniform_rand(0,1);
              age = 0;
              inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n), age, 1 );
            }
          }
        }
      }
    } // end left injector

    //  Inject particles on Right Boundary

    if (global->right) {
      for ( int n=1; n<=nsp; n++ ) {
        species_t * species = find_species_id(n-1,species_list );
        for ( int k=1;k<=nz; k++ ) {
          for ( int j=1;j<=ny; j++ ) {
            bright(n,k,j) = bright(n,k,j) + dt*fright(n,k,j);
            inject = (int) bright(n,k,j);
            bright(n,k,j) = bright(n,k,j) - (double) inject;
            double uflow[3] = {uright(1,n,k,j),uright(2,n,k,j),uright(3,n,k,j)};
            double press[9] = {pright(1,1,n,k,j),pright(1,2,n,k,j),pright(1,3,n,k,j),pright(2,1,n,k,j),pright(2,2,n,k,j),pright(2,3,n,k,j),pright(3,1,n,k,j),pright(3,2,n,k,j),pright(3,3,n,k,j)};
            repeat(inject) {
              compute_injection(uv,nright(n,k,j),uflow,press,-1,2,3,rng);
              x = grid->x1;
              y = grid->y0 + hy*(j-1) + hy*uniform_rand(0,1);
              z = grid->z0 + hz*(k-1) + hz*uniform_rand(0,1);
              age = 0;
              inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n), age, 1 );
            }
          }
        }
      }
    } // end right injector

    //  Inject particles on Top Boundary

    if (global->top) {
      for ( int n=1; n<=nsp; n++ ) {
        species_t * species = find_species_id(n-1,species_list );
        for ( int i=1;i<=nx; i++ ) {
          for ( int j=1;j<=ny; j++ ) {

            vtherm = sqrt(2.0*ptop(3,3,n,i,j)/ntop(n,i,j));
            double t=grid->dt*step;
            double tau = global->tdrive;
            double vexb = (global->edrive)*(1-exp(-t/tau))/field(i,j,nz).cbx;
            vd = vexb/vtherm;
            btop(n,i,j) = btop(n,i,j) + dt*ntop(n,i,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hz);

            //      btop(n,i,j) = btop(n,i,j) + dt*ftop(n,i,j);
            inject = (int) btop(n,i,j);
            btop(n,i,j) = btop(n,i,j)- (double) inject;
            double uflow[3] = {0 ,0 , -vexb};

            //      double uflow[3] = {utop(1,n,i,j),utop(2,n,i,j),utop(3,n,i,j)};
            double press[9] = {ptop(1,1,n,i,j),ptop(1,2,n,i,j),ptop(1,3,n,i,j),ptop(2,1,n,i,j),ptop(2,2,n,i,j),ptop(2,3,n,i,j),ptop(3,1,n,i,j),ptop(3,2,n,i,j),ptop(3,3,n,i,j)};
            repeat(inject) {
              compute_injection(uv,ntop(n,i,j),uflow,press,-3,2,1,rng);
              x = grid->x0 + hx*(i-1) + hx*uniform_rand(0,1) ;
              y = grid->y0 + hy*(j-1) + hy*uniform_rand(0,1);
              z = grid->z1;
              age=0;
              inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n), age, 1 );
            }
          }
        }
      }
    }  // end top injector

    //  Inject particles on Bottom Boundary

    if (global->bottom) {
      for ( int n=1; n<=nsp; n++ ) {
        species_t * species = find_species_id(n-1,species_list );
        for ( int i=1;i<=nx; i++ ) {
          for ( int j=1;j<=ny; j++ ) {

            vtherm = sqrt(2.0*pbot(3,3,n,i,j)/nbot(n,i,j));
            double t=grid->dt*step;
            double tau = global->tdrive;
            double vexb = (global->edrive)*(1-exp(-t/tau))/field(i,j,1).cbx;
            vd =   -vexb/vtherm;
            bbot(n,i,j) = bbot(n,i,j) + dt*nbot(n,i,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hz);

            //      bbot(n,i,j) = bbot(n,i,j) + dt*fbot(n,i,j);


            inject = (int) bbot(n,i,j);
            bbot(n,i,j) = bbot(n,i,j)- (double) inject;
            //      double uflow[3] = {ubot(1,n,i,j),ubot(2,n,i,j),ubot(3,n,i,j)};
            double uflow[3] = {0 ,0 , vexb };
            double press[9] = {pbot(1,1,n,i,j),pbot(1,2,n,i,j),pbot(1,3,n,i,j),pbot(2,1,n,i,j),pbot(2,2,n,i,j),pbot(2,3,n,i,j),pbot(3,1,n,i,j),pbot(3,2,n,i,j),pbot(3,3,n,i,j)};
            repeat(inject) {
              compute_injection(uv,nbot(n,i,j),uflow,press,3,2,1,rng);
              x = grid->x0 + hx*(i-1) + hx*uniform_rand(0,1);
              y = grid->y0 + hy*(j-1) + hy*uniform_rand(0,1);
              z = grid->z0;
              age = 0;
              inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n), age, 1 );
            }
          }
        }
      }
    } // end bottom injector


//  *******  Update the injector moments at every sort interval *********

    double v[3];
    double u[3];
    double p[9];

#define icell(i,j,k) INDEX_FORTRAN_3(i,j,k,0,nx+1,0,ny+1,0,nz+1)
#define v(i) v[INDEX_FORTRAN_1(i,1,3)]
#define u(i) u[INDEX_FORTRAN_1(i,1,3)]
#define p(i,j) p[INDEX_FORTRAN_2(i,j,1,3,1,3)]

//  Parameters for measuring moments  - BE CAREFUL - don't make too big or we will go off the node

 int noff = 1;   // Offset from edge - to measure moments
                 // noff = 0  --> start with cell directly on boundary

 int nav = 2;    // How many cells to include in the "inward" direction in the averging
 int navin = 2;    // How many cells to include in the "inward" direction in the averging

//  Right boundary Moments


 if (global->right) {
   for ( int n=1; n<=nsp; n++ ) {
     species_t * species = find_species_id(n-1,species_list );
     particle_t * part;
     if (remainder(step, global->sort[n-1]) == 0) {
       double npart;
       for ( int k=1;k<=nz; k++ ) {
         for ( int j=1;j<=ny; j++ ) {
           npart = 0;
           flux = 0;
           u[0] = u[1] = u[2] = 0;
           p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
           for ( int i=nx-noff; i>nx-nav-noff; i-- ){
             int nstart = species->partition[icell(i,j,k)];
             int nstop  = species->partition[icell(i,j,k)+1];
             int ncell  = nstop - nstart;
             npart = npart + ncell;
             for (int np=nstart; np<nstop ; np++) {
               part=&species->p[np];
               double gamma = sqrt(1.0+part->ux*part->ux + part->uy*part->uy + part->uz*part->uz);
               v(1) = part->ux;
               v(2) = part->uy;
               v(3) = part->uz;
               if (v(1) < 0) flux = flux - v(1)/gamma;
               for ( int a=1;a<=3; a++ ) {
                 u(a) = u(a) + v(a);
                 for ( int b=1;b<=3; b++ ) {
                   p(a,b) = p(a,b) + v(a)*v(b);
                 }
               }
             } // end particle loop for single cell
           } // end cells included for these moments

           if ( npart > 0 ) {
             fright(n,k,j) = (1.0-rout[1])*fright(n,k,j) + rout[1]*flux/(nav*hx);
             nright(n,k,j) = (1.0-rout[0])*nright(n,k,j) + rout[0]*npart/nav;
             for ( int a=1; a<=3; a++ ) {
               uright(a,n,k,j) = (1-rout[1])*uright(a,n,k,j) + rout[1]*u(a)/npart;
               for ( int b=1;b<=3; b++ ) {
                 p(a,b) = (p(a,b) - u(a)*u(b)/npart)/nav;
                 pright(a,b,n,k,j) = (1-rout[2])*pright(a,b,n,k,j) + rout[2]*p(a,b);
               }
             }
           }
         }
       }
     }
   }
 }  // end right moment update

  //    Left boundary Moments

 if (global->left) {
   for ( int n=1; n<=nsp; n++ ) {
     species_t * species = find_species_id(n-1,species_list );
     particle_t * part;
     if (remainder(step, global->sort[n-1]) == 0) {
       double npart;
       for ( int k=1;k<=nz; k++ ) {
         for ( int j=1;j<=ny; j++ ) {
           npart = 0;
           flux = 0;
           u[0] = u[1] = u[2] = 0;
           p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
           for ( int i=1+noff;i<=nav+noff; i++ ) {
             int nstart = species->partition[icell(i,j,k)];
             int nstop  = species->partition[icell(i,j,k)+1];
             int ncell  = nstop - nstart;
             npart = npart + ncell;
             for (int np=nstart; np < nstop ; np++) {
               part=&species->p[np];
               double gamma = sqrt(1.0+part->ux*part->ux + part->uy*part->uy + part->uz*part->uz);
               v(1) = part->ux;
               v(2) = part->uy;
               v(3) = part->uz;
               if (v(1) > 0) flux = flux + v(1)/gamma;
               for ( int a=1; a<=3; a++ ) {
                 u(a) = u(a) + v(a);
                 for ( int b=1;b<=3; b++ ) {
                   p(a,b) = p(a,b) + v(a)*v(b);
                 }
               }
             } // end particle loop for single cell
           } // end cells included for these moments

           if ( npart > 0 ) {
             fleft(n,k,j) = (1.0-rout[1])*fleft(n,k,j) + rout[1]*flux/(nav*hx);
             nleft(n,k,j) = (1.0-rout[0])*nleft(n,k,j) + rout[0]*npart/nav;
             for ( int a=1; a<=3; a++ ) {
               uleft(a,n,k,j) = (1-rout[1])*uleft(a,n,k,j) + rout[1]*u(a)/npart;
               for ( int b=1;b<=3; b++ ) {
                 p(a,b) = (p(a,b) - u(a)*u(b)/npart)/nav;
                 pleft(a,b,n,k,j) = (1-rout[2])*pleft(a,b,n,k,j) + rout[2]*p(a,b);
               }
             }
           }
         }
       }
     }
   }
}  // end left moment update

  //    Top boundary Moments

 if (global->top) {
   for ( int n=1; n<=nsp; n++ ) {
     species_t * species = find_species_id(n-1,species_list );
     particle_t * part;
     if (remainder(step, global->sort[n-1]) == 0) {
       double npart;
       for ( int i=1;i<=nx; i++ )
         for ( int j=1;j<=ny; j++ ) {
           npart = 0;
           flux = 0;
           u[0] = u[1] = u[2] = 0;
           p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
           for ( int k=nz-noff; k>nz-navin-noff; k-- ) {
             int nstart = species->partition[icell(i,j,k)];
             int nstop  = species->partition[icell(i,j,k)+1];
             int ncell  = nstop - nstart;
             npart = npart + ncell;
             for (int np=nstart; np < nstop ; np++) {
               part=&species->p[np];
               double gamma = sqrt(1.0+part->ux*part->ux + part->uy*part->uy + part->uz*part->uz);
               v(1) = part->ux;
               v(2) = part->uy;
               v(3) = part->uz;
               if (v(3) < 0) flux = flux - v(3)/gamma;
               for ( int a=1; a<=3; a++ ) {
                 u(a) = u(a) + v(a);
                 for ( int b=1;b<=3; b++ ) {
                   p(a,b) = p(a,b) + v(a)*v(b);
                 }
               }
             } // end particle loop for single cell
           } // end cells included for these moments

           if ( npart > 0 ) {
             ftop(n,i,j) = (1.0-rin[1])*ftop(n,i,j) + rin[1]*flux/(navin*hz);

//           if ( npart/navin < ntop(n,i,j) ) ftop(n,i,j) = (ntop(n,i,j)- npart/navin)/dt; else ftop(n,i,j)=0;

             ntop(n,i,j) = (1-rin[0])*ntop(n,i,j) + rin[0]*npart/navin;
             for ( int a=1;a<=3; a++ ) {
               utop(a,n,i,j) = (1-rin[1])*utop(a,n,i,j) + rin[1]*u(a)/npart;
               for ( int b=1;b<=3; b++ ) {
                 p(a,b) = (p(a,b) - v(a)*v(b)/npart)/navin;
                 ptop(a,b,n,i,j) = (1-rin[2])*ptop(a,b,n,i,j) + rin[2]*p(a,b);
               }
             }
           }
         }
     }
   }
 }  // end bottom moment update

 if (global->bottom) {
   for ( int n=1; n<=nsp; n++ ) {
     species_t * species = find_species_id(n-1,species_list );
     particle_t * part;
     if (remainder(step, global->sort[n-1]) == 0) {
       double npart;
       for ( int i=1;i<=nx; i++ )
         for ( int j=1;j<=ny; j++ ) {
           npart = 0;
           flux = 0;
           u[0] = u[1] = u[2] = 0;
           p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
           for ( int k=1+noff;k<=navin+noff; k++ ) {
             int nstart = species->partition[icell(i,j,k)];
             int nstop  = species->partition[icell(i,j,k)+1];
             int ncell  = nstop - nstart;
             npart = npart + ncell;
             for (int np=nstart; np < nstop ; np++) {
               part=&species->p[np];
               double gamma = sqrt(1.0+part->ux*part->ux + part->uy*part->uy + part->uz*part->uz);
               v(1) = part->ux;
               v(2) = part->uy;
               v(3) = part->uz;
               if (v(3) > 0) flux = flux + v(3)/gamma;
               for ( int a=1; a<=3; a++ ) {
                 u(a) = u(a) + v(a);
                 for ( int b=1;b<=3; b++ ) {
                   p(a,b) = p(a,b) + v(a)*v(b);
                 }
               }
             } // end particle loop for single cell
           } // end cells included for these moments

           if ( npart > 0) {
             fbot(n,i,j) = (1.0-rin[1])*fbot(n,i,j) + rin[1]*flux/(navin*hz);
             //      if ( npart/navin < nbot(n,i,j) ) fbot(n,i,j) = (nbot(n,i,j)- npart/navin)/dt; else fbot(n,i,j)=0;
             nbot(n,i,j) = (1-rin[0])*nbot(n,i,j) + rin[0]*npart/navin;
             for ( int a=1;a<=3; a++ ) {
               ubot(a,n,i,j) = (1-rin[1])*ubot(a,n,i,j) + rin[1]*u(a)/npart;
               for ( int b=1;b<=3; b++ ) {
                 p(a,b) = (p(a,b) - v(a)*v(b)/npart)/navin;
                 pbot(a,b,n,i,j) = (1-rin[2])*pbot(a,b,n,i,j) + rin[2]*p(a,b);
               }
             }
           }
         }
     }
   }
 }  // end bottom moment update


//  Periodically save injector moments on outflow boundaries
//  Only do this on the outflow boundaries, and only if we have
//  a single domain on each boundary - otherwise would have to combine data files.

 if ( global->topology_y == 1 && global->topology_z ==1 ) {

// How often to write moments to file -

   int nskip = 10;

   if (global->left) {
     int j = 1;
     for ( int n=1; n<=nsp; n++ ) {
       if (remainder(step, nskip*global->sort[n-1]) == 0) {
         char buffer[20];
         sprintf(buffer, "injectors/left%i.dat", n);
         FileIO fileIO;
         FileIOStatus status;
         status= fileIO.open(buffer, io_append);
         for ( int k=1;k<=nz; k++ ) {
           zcell = (grid->z0 + k*hz-hz/2)/(global->L_de);
           fileIO.print("%6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g \n", zcell,nleft(n,k,j),uleft(1,n,k,j),uleft(2,n,k,j),uleft(3,n,k,j),pleft(1,1,n,k,j),pleft(2,2,n,k,j),pleft(3,3,n,k,j),pleft(1,2,n,k,j),pleft(1,3,n,k,j),pleft(2,3,n,k,j));
         } //end for
         fileIO.print("  \n \n");
         fileIO.close();
       }
     }
   }  //end left output

   if (global->right) {
     int j = 1;
     for ( int n=1; n<=nsp; n++ ) {
       if (remainder(step, nskip*global->sort[n-1]) == 0) {
         char buffer[20];
         sprintf(buffer, "injectors/right%i.dat", n);
         FileIO fileIO;
         FileIOStatus status;
         status= fileIO.open(buffer, io_append);
         for ( int k=1;k<=nz; k++ ) {
           zcell = (grid->z0 + k*hz-hz/2)/(global->L_de);
           fileIO.print("%6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g %6.4g \n", zcell,nright(n,k,j),uright(1,n,k,j),uright(2,n,k,j),uright(3,n,k,j),pright(1,1,n,k,j),pright(2,2,n,k,j),pright(3,3,n,k,j),pright(1,2,n,k,j),pright(1,3,n,k,j),pright(2,3,n,k,j));
         } //end for
         fileIO.print("  \n \n");
         fileIO.close();
       }
     }
   }  //end right output

 } //end output for injector moments

} // end particle injection

//   *******************  CURRENT INJECTION ***************************

begin_current_injection {

  // No current injection for this simulation

}

//   *******************  FIELD INJECTION ***************************

begin_field_injection {
  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;
  double t=grid->dt*step;
  double tau = global->tdrive;
  int x,y,z;

  //   There macros are from local.c to apply boundary conditions

#define XYZ_LOOP(xl,xh,yl,yh,zl,zh)             \
  for( z=zl; z<=zh; z++ )                       \
    for( y=yl; y<=yh; y++ )                     \
      for( x=xl; x<=xh; x++ )

#define xy_EDGE_LOOP(z) XYZ_LOOP(1,nx,1,ny+1,z,z)
#define yx_EDGE_LOOP(z) XYZ_LOOP(1,nx+1,1,ny,z,z)

  // Top Boundary

  if (global->top) {

    yx_EDGE_LOOP(nz+1) field(x,y,z).ey = (global->edrive)*(1-exp(-t/tau));
    xy_EDGE_LOOP(nz+1) field(x,y,z).ex = -(global->edrive)*(1-exp(-t/tau))*global->bg;

  }

  // Bottom Boundary

  if (global->bottom) {

    yx_EDGE_LOOP(1) field(x,y,z).ey = (global->edrive)*(1-exp(-t/tau));
    xy_EDGE_LOOP(1) field(x,y,z).ex = (global->edrive)*(1-exp(-t/tau))*global->bg;

  }

}  // end field injection

//  Collisional algorithm
#include "collisions"
