////////////////////////////////////////////////////////////////
//
//   Moon simulation with simple parameters
//
//   Start with initially magnetised plasma ExB drifting from the
//   left to the right
//
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////

#include "injection.cp"   //  Routines to compute re-injection velocity

//////////////////////////////////////////////////////

#define NUM_TURNSTILES 4096

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
  int write_restart;
  int write_end_restart;
  double quota_sec;        // Run quota in seconds
  double b0;               // B0 (IMF)
  double sn;               // sin theta (angle between IMF and Vsw toward +x)
  double cs;               // cos theta

  int dis_nav;         //  : number of steps to average over
  int dis_interval;    //  : number of steps between outputs. 
  int dis_iter;        // : iteration count. 0 means we are not avraging at the moment 
  int dis_begin_int;   // : the first time step of the interval

  double topology_x;       // domain topology 
  double topology_y;
  double topology_z;

  //  Variables for Open BC Model
  int nsp;          //  Number of Species
  double rout[3]; // Relaxation parameter for outlfow boundary moments
  double sort[2]; // Intervals where we know particles are sorted
  double npleft[2];  // Left Densities
  double npright[2];  // Right Densities
  double vth[2];    // Thermal velocity 
  double q[2];      // Species charge
  double ur;       //  Fluid velociy on right
  double ul;       //  Fluid velociy on left
  double nfac;      //  Normalization factor to convert particles per cell into density

  int left,right,top_y,bot_y;  // Keep track of boundary domains
  double *nleft, *uleft, *pleft, *bleft, *fleft;     // Moments for left injectors
  double *nright, *uright, *pright, *bright, *fright; // Moments for right injectors
  double *ntop_y, *utop_y, *ptop_y, *btop_y, *ftop_y;     // Moments for left injectors
  double *nbot_y, *ubot_y, *pbot_y, *bbot_y, *fbot_y;     // Moments for left injectors


  // Some params for LMM 
  double Dp;

  double x_dp; // x-location of dipole
  double y_dp; // y-location of dipole

  double c_d;  // distance between wires
  double c_r;  // radius of each wire
  double c_j0; // current density in each wire
  double t_dipole_start; // time to start ramping dipole
  double t_dipole_rise; // time scale for dipole ramp up

//  Variables for new output format

  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;
};

begin_initialization {

 // Use natural PIC units 

  double ec   = 1;         // Charge normalization
  double me   = 1;         // Mass normalization
  double c    = 1;         // Speed of light
  double de   = 1;         // Length normalization (electron inertial length) 
  double wpe  = 1;          // Time normalized by electron plasma frequency
  double eps0 = 1;         // Permittivity of space
  double n0 = 1.0;        // Scale to normalize density

 //  Note - these choices imply that time is normalized to wpe

  //  Some basic numerical parameters

  double cfl_req   = 0.99;  // How close to Courant should we try to run
  double wpedt_max = 0.36;  // How big a timestep is allowed if Courant is not too restrictive
  double damp      = 0.0;   // Level of radiation damping
  int rng_seed     = 1;     // Random number seed increment

  // Set Physics parameters

  double mime   = 1836.0;//100.0;  // Ion mass / electron mass
  double wpe_wce = 2; // elecron plasma freq/ cyclotron freq
  double M_A = 6.5; // Alfven Mach Number
  double beta_e = 0.6;//0.5;   // beta of solar wind electrons
  double beta_i = 0.6;//0.383;   // beta of upstream ions (solar wind)
  double theta = -30.0*M_PI/180.; // Angle between B field and +x unit vector

  double DP_di = 0.55; // Stand-off distance in d_i.

  double cs = cos(theta);
  double sn = sin(theta); 

  // Simulation time and quota

  double tau    = 20;    // Maximum simulation time in units of t*wci
  double quota   = 15.5;   // run quota in hours
  double quota_sec = quota*3600;  // Run quota in seconds
  
  // Derived qunatities

  double mi = me*mime;        // Ion mass
  double wci = 1.0/(mime*wpe_wce);
  double wce = wpe/wpe_wce; 
  double mu0 = 1/(eps0*c*c);
  double b0 = me*c*wce/ec;         //  Initial magnetic field
  // Beta is 2*mu0*n0*T/B0^2
  // c^2 = 1/(mu0eps0) => mu0 = 1/(eps0*c*c)
  double Te = b0*b0/(2.0*n0*(1/(eps0*c*c)))*beta_e; // Electron temperature in units of me*c**2 in plasma rest frame
  double Ti = b0*b0/(2.0*n0*(1/(eps0*c*c)))*beta_i;// Ion temperature in units of me*c**2 in plasma rest frame
  double vthe = sqrt(Te/me);  // Electron thermal velocity 
  double vthi = sqrt(Ti/mi);  // Ion thermal velocity 
  double v_A = b0/sqrt(n0*mi*mu0); 
  double Vflow = M_A*v_A;  // Flow Velocity
  double gam = 1/sqrt(1-Vflow*Vflow/c/c);     // Gamma of plasma flow
  double wpi  = wpe/sqrt(mime);   // ion plasma frequency
  double di   = c/wpi;            // ion inertial length
  double rhoi = vthi/wci;       // ion gyroradius
  double rhoe = vthe/wce;      // electron gyroradius

  double rout[3] = {0.000, 0.000, 0.000};

  double debye = vthe;
  double R_moon = 27.0*di;
  double x_moon = 25.0*di;

  double x_dp = x_moon-18.93*di;//28.0*di;
  double y_dp = -18.93*di;
  double c_d  = 0.15*di;
  double c_r  = 0.05*di;
  double Bstar = b0*M_A; // CHECK!
  double Dp = DP_di*di;
  double c_j0 = 2.0*Dp*Dp/c_d*Bstar/(c_r*c_r); // CHECK!

  double t_dipole_start = 5.0/wci; 
  double t_dipole_rise = 5.0/wci;

  // Simulation  parameters

  double nppc = 50;        // Average number of macro particle per cell per species

  double Lx  = 12.5*di;   // size of box in x dimension
  double ybot = -25*di;
  double ytop = -12.5*di;
  double Ly  = ytop-ybot;     // size of box in y dimension  
  double Lz  = 100*de;     // size of box in z dimension  
  
  double topology_x = 1024;  // Number of domains in x, y, and z
  double topology_y = 8; 
  double topology_z = 1; 

  double nx = 8192;   //  Number of cells in x-direction
  double ny = 8192;     //  Number of cells in y-direction
  double nz = 1;   //  Number of cells in z-direction

  double hx = Lx/nx;   // cell size in x
  double hy = Ly/ny;   // cell size in y
  double hz = Lz/nz;   // cell size in z
  
  // For 2D - set the transverse cell sizes to be same

  //  hy = hx;
  //  hy = hx;
  hz = hx;
  //  Ly = hx;
  Lz = hx;

  double Npart   = nppc*nx*ny*nz;  // total macro electrons in box
  Npart = trunc_granular(Npart,nproc()); // Make it divisible by number of processors
  double qe = -ec*Lx*Ly*Lz/Npart;  // Charge per macro electron
  double qi =  ec*Lx*Ly*Lz/Npart;  // Charge per macro electron       
  double nfac = qi/(hx*hy*hz);    // Convert density to particles per cell

  double dg = courant_length(Lx,Ly,Lz,nx,ny,nz);        // Courant length
  double dt = cfl_req*dg/c;                             // Courant limited time step
  dt = hx/2.0;
  if( wpe*dt>wpedt_max) dt=wpedt_max/wpe;               // override timestep if plasma frequency limited

  //  Intervals for various output and checks

  int restart_interval = 4000;
  int energies_interval = 500;
  int interval = int(0.05/(wci*dt)); //FIXME
  int fields_interval = interval;
  int ehydro_interval = interval;
  int Hhydro_interval = interval;
  int eparticle_interval = int(interval*50);
  int Hparticle_interval = int(interval*50);
  int quota_check_interval = 100;

  ///////////////////////////////////////////////
  // Setup high level simulation parameters
  num_step             = int(tau/(wci*dt)); //Fixme
  status_interval      = 200;
  sync_shared_interval = status_interval/2;
  clean_div_e_interval = 2*status_interval;
  clean_div_b_interval = 2*status_interval;

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
  global->sn  = sn;
  global->cs  = cs;

  global->topology_x  = topology_x;  
  global->topology_y  = topology_y;  
  global->topology_z  = topology_z;  

 // From grid/partition.c: used to determine which domains are on edge
 // This will return the logical "brick" layout of the MPI domains (ix,iy,iz)

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

  int ix, iy, iz, left=0,right=0,top_y=0,bot_y=0;
  RANK_TO_INDEX( int(rank()), ix, iy, iz ); 
  if ( ix ==0 ) left=1;
  if ( ix ==topology_x-1 ) right=1;

  if ( iy == 0 ) bot_y=1;
  if ( iy == topology_y-1 ) top_y = 1;

  //  Parameters for the open boundary model 

  global->nsp = 2;

  global->rout[0]  = rout[0];
  global->rout[1]  = rout[1];
  global->rout[2]  = rout[2];

  global->ur  = Vflow;  
  global->ul  = Vflow;  

  global->q[0]  = -qe;  
  global->q[1]  = qi;  

  global->npleft[0]  = n0;
  global->npleft[1]  = n0;

  global->npright[0]  = n0;
  global->npright[1]  = n0;

  global->vth[0]  = sqrt(2.)*vthe;
  global->vth[1]  = sqrt(2.)*vthi;

  global->left = left;
  global->right = right;
  global->bot_y = bot_y;
  global->top_y = top_y;
  global->nfac = nfac; 


  // Parameters for LMM
  global->Dp = Dp;
  global->x_dp = x_dp;
  global->y_dp = y_dp;
  global->c_d = c_d;
  global->c_r = c_r;
  global->c_j0 = c_j0;
  global->t_dipole_start = t_dipole_start;
  global->t_dipole_rise = t_dipole_rise;


  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the grid

  // Setup basic grid parameters
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = c;
  grid->eps0 = eps0;

  // Define periodic the grid
  
  define_periodic_grid(  0, ybot, -0.5*Lz,    // Low corner
			  Lx, ytop, 0.5*Lz,     // High corner
			  nx, ny, nz,             // Resolution
			  topology_x, topology_y, topology_z); // Topology

  //  Over-ride the peridoic boundary conditions as desired:

  // Set the boundary at x=0 to be absorbing?

   if ( ix==0 ) set_domain_field_bc( BOUNDARY( -1,0,0), pec_fields ); 
   if ( ix==0 ) set_domain_particle_bc( BOUNDARY( -1,0,0), absorb_particles ); 

  // Set the boundary at x=Lx to be absorbing for particles and fields

  if ( ix==topology_x-1 ) set_domain_field_bc( BOUNDARY(1,0,0), absorb_fields ); 
  if ( ix==topology_x-1 ) set_domain_particle_bc( BOUNDARY(1,0,0), absorb_particles ); 


  // Set the boundary at x=Lx to be absorbing for particles and fields

  if ( iy==0 ) set_domain_field_bc( BOUNDARY(0,-1,0), pec_fields ); 
  if ( iy==0 ) set_domain_particle_bc( BOUNDARY(0,-1,0), absorb_particles ); 

  if ( iy==topology_y-1 ) set_domain_field_bc( BOUNDARY(0,1,0), pec_fields ); 
  if ( iy==topology_y-1 ) set_domain_particle_bc( BOUNDARY(0,1,0), absorb_particles );

  // we should have open BC on y and z boundaries. Though periodic seems ok

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the species

  sim_log("Setting up species. ");

  double electron_sort_interval = 25;
  double ion_sort_interval = 25;

  double nmax = 20.0*Npart/nproc();
  double nmovers = 0.1*nmax;
  double sort_method = 0; 
  species_t *electron = define_species("electron",-ec, me, nmax, nmovers, electron_sort_interval, sort_method); 
  species_t *ion = define_species("ion",ec, mi, nmax, nmovers, ion_sort_interval, sort_method); 

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup materials

  sim_log("Setting up materials. "); 

  define_material( "vacuum", 1 );
  //  define_units(3e8 , 8.85e-12);
  // what is the moon?
  material_t * moonstone = define_material( "moon", 1, 1, 0);//#10000);
  material_t * resistive = define_material( "resistive",1,1,0);//10000);

#define MOON ( (x-x_moon)*(x-x_moon) + y*y + z*z < R_moon*R_moon)

# define top_wire (x-global->x_dp-global->c_d/2)*(x-global->x_dp-global->c_d/2)+(y-global->y_dp)*(y-global->y_dp) < global->c_r*global->c_r
# define bot_wire (x-global->x_dp+global->c_d/2)*(x-global->x_dp+global->c_d/2)+(y-global->y_dp)*(y-global->y_dp) < global->c_r*global->c_r

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.

////////////////////////////////////////////////////////////////////////////////////////////
//  Finalize Field Advance

  sim_log("Finalizing Field Advance"); 

  define_field_array(NULL); // second argument is damp, default to zero


  set_region_bc( MOON, absorb_particles, absorb_particles, absorb_particles ); 
  set_region_material( MOON, moonstone, moonstone ); 

  double thickness = 0;
#define resistive_layer ( x > Lx-hx*thickness)//			\
//                         || z <-Lz/2+hz*thickness  || z > Lz/2-hz*thickness )

  if (thickness > 0)
    set_region_material(resistive_layer, resistive, resistive);

  sim_log( "Finished setting moon material" );
  ///////////////////////////////////////////////////
  // Log diagnostic information about this simulation

  sim_log( "***********************************************" );
  sim_log("* Topology:                       "<<topology_x<<" "<<topology_y<<" "<<topology_z); 
  sim_log ( "mi/me = " << mime );
  sim_log ( "beta_e = " << beta_e );
  sim_log ( "beta_i = " << beta_i );
  sim_log ( "Te/(me*c**2) = " << Te );
  sim_log ( "Ti/(me*c**2) = " << Ti );
  sim_log ( "wpe/wce = " << wpe_wce );
  sim_log ( "wpi/wci = " << wpi/wci );
  sim_log ( "Gamma = " << gam );
  sim_log ( "Vflow/c = " << Vflow );
  sim_log ( "V_A/c = " << v_A );
  sim_log ( "tau*wpe = " << tau );
  sim_log ( "wci = " << wci );
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
  sim_log ( " Npart = " << Npart );
  sim_log ( "total # of particles = " << 2*Npart );
  sim_log ( "dt*wpe = " << wpe*dt ); 
  sim_log ( " energies_interval: " << energies_interval );
  sim_log ( "dx/de = " << Lx/(de*nx) );
  sim_log ( "dy/de = " << Ly/(de*ny) );
  sim_log ( "dz/de = " << Lz/(de*nz) );
  sim_log ( "dx/debye = " << (Lx/nx)/(vthe/wpe)  );
  sim_log ( "vthi/c = " << vthi/c );
  sim_log ( "vthe/c = " << vthe/c );
  sim_log ( "dx/rhoe = " << Lx/(rhoe*nx) );
  sim_log ( "dx/rhoi = " << Lx/(rhoi*nx) );
  
  // Dump simulation information to file "info"
  if (rank() == 0 ) {
    FILE *fp_info;
    if ( ! (fp_info=fopen("info", "w")) ) ERROR(("Cannot open file."));
    fprintf(fp_info, "           ***** Simulation parameters ***** \n");
    fprintf(fp_info, "		 mi/me =		%e\n", mime );
    fprintf(fp_info, "		 beta_e =		%e\n", beta_e );
    fprintf(fp_info, "		 beta_i =		%e\n", beta_i );
    fprintf(fp_info, "		 Te/(me*c**2) =		%e\n", Te );
    fprintf(fp_info, "		 Ti/(me*c**2) =		%e\n", Ti );
    fprintf(fp_info, "		 wpe/wce =		%e\n", wpe_wce );
    fprintf(fp_info, "		 Gamma = 		%e\n", gam );
    fprintf(fp_info, "		 Vflow/c = 		%e\n", Vflow );
    fprintf(fp_info, "		 V_A/c = 		%e\n", v_A );
    fprintf(fp_info, "		 tau*wpe =		%e\n", tau );
    fprintf(fp_info, "		 wci =		%e\n", wci );
    fprintf(fp_info, "		 num_step = 		%i\n", num_step );
    fprintf(fp_info, "		 Lx/de = 		%e\n", Lx/de );
    fprintf(fp_info, "		 Ly/de = 		%e\n", Ly/de );
    fprintf(fp_info, "		 Lz/de =		%e\n", Lz/de );
    fprintf(fp_info, "		 Lx/di = 		%e\n", Lx/di );
    fprintf(fp_info, "		 Ly/di = 		%e\n", Ly/di );
    fprintf(fp_info, "		 Lz/di =		%e\n", Lz/di );
    fprintf(fp_info, "		 nx = 			%e\n", nx );
    fprintf(fp_info, "		 ny = 			%e\n", ny );
    fprintf(fp_info, "		 nz =			%e\n", nz );
    fprintf(fp_info, "		 damp =			%e\n", damp );
    fprintf(fp_info, "		 courant = 		%e\n", c*dt/dg );
    fprintf(fp_info, "		 nproc = 		%e\n", nproc() );
    fprintf(fp_info, "		 nppc = 		%e\n", nppc );
    fprintf(fp_info, "		 b0 =			%e\n", b0 );
    fprintf(fp_info, "		 di = 			%e\n", di );
    fprintf(fp_info, "		 Npart =        	%e\n", Npart );
    fprintf(fp_info, "		 total # of particles = %e\n", 2*Npart );
    fprintf(fp_info, "		 dt*wpe = 		%e\n", wpe*dt );
    fprintf(fp_info, "		 energies_interval: 	%i\n", energies_interval);
    fprintf(fp_info, "		 dx/de =		%e\n", Lx/(de*nx) );
    fprintf(fp_info, "		 dy/de =		%e\n", Ly/(de*ny) );
    fprintf(fp_info, "		 dz/de =		%e\n", Lz/(de*nz) );
    fprintf(fp_info, "		 dx/debye = 		%e\n", (Lx/nx)/(vthe/wpe) );
    fprintf(fp_info, "		 vthi/c =		%e\n", vthi/c );
    fprintf(fp_info, "		 vthe/c =		%e\n", vthe/c );
    fprintf(fp_info, "		 dx/rhoe =		%e\n", Lx/(rhoe*nx) );
    fprintf(fp_info, "		 dx/rhoi =		%e\n", Lx/(rhoi*nx) );
    fprintf(fp_info, "		 ***************************\n");
    fclose(fp_info);
  }

    // for the fortran translation routine
    // write binary info file
                
    FileIO fp_info;
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
    fp_info.write(&mime, 1 );
    fp_info.write(&vthe, 1 );
    fp_info.write(&vthi, 1 );    
    fp_info.close();

  ////////////////////////////
  // Load fields

  sim_log( "Loading fields" );
  set_region_field( everywhere, 0, 0, Vflow*-b0*sn,       // Electric field
		               b0*cs, b0*sn ,0  );    // Magnetic field

  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specied as logical equations (i.e. x>0 && x+y<2)

  // Wondering what E in the moon is.
  set_region_field( MOON, 0, 0, Vflow*-b0*sn,       // Electric field
		               b0*cs, b0*sn ,0  );    // Magnetic field

  // LOAD PARTICLES

  sim_log( "Loading particles" );

  // Do a fast load of the particles

  seed_entropy( rank() ); 


  double xmin = grid->x0 , xmax = grid->x0+(grid->dx)*(grid->nx);
  double ymin = grid->y0 , ymax = grid->y0+(grid->dy)*(grid->ny);
  double zmin = grid->z0 , zmax = grid->z0+(grid->dz)*(grid->nz);

  sim_log("abs(qe)=" << abs(qe) );

  double x,y,z,ux,uy,uz,ep;
  repeat ( Npart/nproc() ) {

    // Now choose random position 

    z = uniform(rng(0), zmin, zmax);
    x = uniform(rng(0), xmin, xmax);
    y = uniform(rng(0), ymin, ymax);
    //  Choose random momentum in rest frame of flow

    ux = normal(rng(0), 0, vthe);
    uy = normal(rng(0), 0, vthe); 
    uz = normal(rng(0), 0, vthe); 

    // Now boost particle momentum into flow frame

    ep = sqrt(1.0 + ux*ux + uy*uy + uz*uz);   // particle energy in rest frame of the flow
    ux = gam*(ux + Vflow*ep);  // Boost x-component of momentum

    // Inject electron into the simulation
    inject_particle( electron, x, y, z, ux, uy, uz, -qe, 0, 0); 

    //  Pick a new momentum for the ion - but load at the same position as electron

    //    x = uniform_rand(xmin,xmax);
    //y = uniform_rand(ymin,ymax);
    //z = uniform_rand(zmin,zmax);

    ux = normal(rng(0), 0, vthi)+Vflow;
    uy = normal(rng(0), 0, vthi);
    uz = normal(rng(0), 0, vthi);


    //ep = sqrt(1.0 + ux*ux + uy*uy + uz*uz);   // particle energy in rest frame of the flow
    //ux = gam*(ux + Vflow*ep);  // Boost x-component of momentum

    inject_particle( ion, x, y, z, ux, uy, uz, qi, 0, 0);

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

	//  These have just the most useful things turned on
	//   global->fdParams.output_variables( electric | magnetic | current | div_e_err );
	//    global->hedParams.output_variables( current_density | charge_density | ke_density | stress_tensor);
	//    global->hHdParams.output_variables( current_density | charge_density | ke_density | stress_tensor);

	
	// Here we are dumping everthing

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
        (global->x##_interval>0 && remainder(step(), global->x##_interval) == 0)


begin_diagnostics {
  global->eparticle_interval=global->fields_interval;
  
	  const int nsp = global->nsp;
	  const int nx = grid->nx;
	  const int ny = grid->ny;
	  const int nz = grid->nz;

	  if ( (step()%100)==0 ) sim_log( "Time step: " << step()); 

  /*--------------------------------------------------------------------------
   * Normal rundata dump
   *------------------------------------------------------------------------*/
	if(step() == 0) {

		dump_mkdir("fields");
		dump_mkdir("hydro");
		dump_mkdir("rundata");
		dump_mkdir("data");
		dump_mkdir("restore0");
		dump_mkdir("restore1");  // 1st backup
		dump_mkdir("restore2");  // 2nd backup
		dump_mkdir("particle");

		dump_grid("rundata/grid");
		dump_materials("rundata/materials");
		dump_species("rundata/species");
		global_header("global", global->outputParams);
	} // if


#include "time_average_master.cxx"

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
	 * Restart dump
	 *------------------------------------------------------------------------*/

	if(step() && !(step()%global->restart_interval)) {
	  global->write_restart = 1; // set restart flag. the actual restart files are written during the next step
	} else {
	  if (global->write_restart) {
	    global->write_restart = 0; // reset restart flag

	    double dumpstart = uptime();
	    BEGIN_TURNSTILE(NUM_TURNSTILES) { 
	      if(!global->rtoggle) {
		global->rtoggle = 1;
		checkpt("restore1/restore", 0);
		DUMP_INJECTORS(1);
	      }
	      else {
		global->rtoggle = 0;
		checkpt("restore2/restore", 0);
		DUMP_INJECTORS(2);
	      } // if
	    } END_TURNSTILE;

	    double dumpelapsed = uptime() - dumpstart;
	    sim_log("restart duration " << dumpelapsed); 
	  } // if

	}

  // Dump particle data
	// #ifdef FOURZEROSEVEN
	//   if ( should_dump(eparticle) ) {
	//dump_particles("electron",  "particles/electron");
	// }
	//if ( should_dump(Hparticle) ) {
	// dump_particles("ion",  "particles/ion");    
	// }
	// #else
	char subdir[36]; 
	if (should_dump(eparticle) && step() != 0) { 
	  BEGIN_TURNSTILE(NUM_TURNSTILES) { 
	    sprintf(subdir,"particle/T.%d", step()); 
	    dump_mkdir(subdir); 

	    sprintf(subdir,"particle/T.%d/electron", step()); 
	    dump_particles("electron", subdir);

	    sprintf(subdir,"particle/T.%d/ion", step()); 
	    dump_particles("ion", subdir);

	  } END_TURNSTILE; 
	}
	// #endif 

  // Shut down simulation when wall clock time exceeds global->quota_sec. 
  // Note that the mp_elapsed() is guaranteed to return the same value for all
  // processors (i.e., elapsed time on proc #0), and therefore the abort will 
  // be synchronized across processors. Note that this is only checked every
  // few timesteps to eliminate the expensive mp_elapsed call from every
  // timestep. mp_elapsed has an ALL_REDUCE in it!

	/*  if( step()>0 && global->quota_check_interval>0 && (step()&global->quota_check_interval)==0 ) {
     if( uptime() > global->quota_sec ) {
	  sim_log( "Allowed runtime exceeded for this job.  Terminating....\n");
	  double dumpstart = uptime(); 
	  BEGIN_TURNSTILE(NUM_TURNSTILES) {  
	    checkpt("restore0/restore",0);
	    DUMP_INJECTORS(0);
	  } END_TURNSTILE; 
	  barrier( ); // Just to be safe
	  sim_log( "Restart dump restart completed." );
	  double dumpelapsed = uptime() - dumpstart; 
	  sim_log("Restart duration " <<dumpelapsed); 
	  finalize(  ); 
	  exit(0); // Exit or abort?
     }
     }*/
	if ( ( step() % 400 == 0 ) || (global->write_end_restart) ) {

	  if ( (global->write_end_restart) ) {
	    global->write_end_restart = 0; // reset restart flag
                                                                                                                                                                                                                
	    sim_log( "Allowed runtime exceeded for this job.  Terminating....\n");
	    double dumpstart = uptime();

	    BEGIN_TURNSTILE(NUM_TURNSTILES) {
	      checkpt( "restart0/checkpt", step() );
	    } END_TURNSTILE;

	    FileIO last_restart;
	    last_restart.open("last_restart", io_write);
	    last_restart.print("%d\n",step());
	    last_restart.close();

	    barrier();

	    sim_log( "Restart dump restart completed." );
	    double dumpelapsed = uptime() - dumpstart;
	    sim_log("Restart duration "<< dumpelapsed);

	    finalize();
	    exit(0);
	  }

	  if( uptime() > global->quota_sec ) global->write_end_restart = 1;
	}



} // end diagnostics


 
// *******************  PARTICLE INJECTION  - OPEN BOUNDARY ***************************
 
begin_particle_injection {
  int inject;
  double x, y, z, age, vtherm, vd;
  double uv[3];
  double flux;
  double nfac = global->nfac;
  const int nsp=global->nsp;
  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;
  const double sqpi =1.772453850905516;
  const double dt=grid->dt;
  const double hx=grid->dx;
  const double hy=grid->dy;
  const double hz=grid->dz;
  const double rout[3]={global->rout[0],global->rout[1],global->rout[2]};

  // Initialize the injectors on the first call

    static int initted=0;
    if ( !initted ) {

      initted=1;

      if (rank() == 0) MESSAGE(("----------------Initializing the Particle Injectors-----------------")); 
      
      // MESSAGE(("------rank=%g    right=%i     left=%i    nsp=%i",rank(),global->right,global->left,nsp)); 
      // Intialize injectors 

      if (global->right) {
	      if (rank() == 0) MESSAGE(("----------------Initializing the Right Particle Injectors-----------------")); 
	DEFINE_INJECTOR(right,ny,nz);
	if (step() == 0) { 
	  for ( int n=1; n<=nsp; n++ ) { 
	    for ( int k=1;k<=nz; k++ ) {
	      for ( int j=1;j<=ny; j++ ) { 
		bright(n,k,j) = 0;
		nright(n,k,j) = npright(n)/nfac;
		uright(1,n,k,j) = global->ur;
		uright(2,n,k,j) = 0;
		uright(3,n,k,j) = 0;
		pright(1,2,n,k,j)=pright(2,1,n,k,j)=pright(1,3,n,k,j)=pright(3,1,n,k,j)=pright(2,3,n,k,j)=pright(3,2,n,k,j)=0;
		pright(1,1,n,k,j) = npright(n)*vth(n)*vth(n)/(2.0*nfac);
		pright(2,2,n,k,j) = pright(1,1,n,k,j);
		pright(3,3,n,k,j) = pright(1,1,n,k,j);
	      }      
	    }
	  }  // end for	
	} // endif
	else {

      if (rank() == 0) MESSAGE(("----------------Reading the Particle Injectors-----------------")); 
      READ_INJECTOR(right, ny, nz, 0);
	}
      } //end right boundary
	       
      if (global->left) {

	DEFINE_INJECTOR(left,ny,nz);
	if (step() == 0) {
	  for ( int n=1; n<=nsp; n++ ) {
	    for ( int k=1;k<nz+1; k++ ) {
	      for ( int j=1;j<=ny; j++ ) {   
		bleft(n,k,j) = 0;
		nleft(n,k,j) =  npleft(n)/nfac;
		uleft(1,n,k,j) = global->ul;
		uleft(2,n,k,j) = 0;
		uleft(3,n,k,j) = 0;
		pleft(1,2,n,k,j)=pleft(2,1,n,k,j)=pleft(1,3,n,k,j)=pleft(3,1,n,k,j)=pleft(2,3,n,k,j)=pleft(3,2,n,k,j)=0;
		pleft(1,1,n,k,j) = npleft(n)*vth(n)*vth(n)/(2.0*nfac);
		pleft(2,2,n,k,j) = pleft(1,1,n,k,j);
		pleft(3,3,n,k,j) = pleft(1,1,n,k,j);
	      }     
	    } 
	  } // end for
	} //endif
	else { 
	  READ_INJECTOR(left, ny, nz, 0); 
	}
      } // end left boundary
	    

      if (global->top_y) {
	MESSAGE(("----------------Initializing the Top y Particle Injectors-----------------")); 
	DEFINE_INJECTOR(top_y,nx,nz);
	if (step() == 0) { 
	  for ( int n=1; n<=nsp; n++ ) { 
	    for ( int k=1;k<=nz; k++ ) {
	      for ( int j=1;j<=nx; j++ ) { 
		btop_y(n,k,j) = 0;
		ntop_y(n,k,j) = npright(n)/nfac;
		utop_y(1,n,k,j) = global->ur;
		utop_y(2,n,k,j) = 0;
		utop_y(3,n,k,j) = 0;
		ptop_y(1,2,n,k,j)=ptop_y(2,1,n,k,j)=ptop_y(1,3,n,k,j)=ptop_y(3,1,n,k,j)=ptop_y(2,3,n,k,j)=ptop_y(3,2,n,k,j)=0;
		ptop_y(1,1,n,k,j) = npright(n)*vth(n)*vth(n)/(2.0*nfac);
		ptop_y(2,2,n,k,j) = ptop_y(1,1,n,k,j);
		ptop_y(3,3,n,k,j) = ptop_y(1,1,n,k,j);
	      }      
	    }
	  }  // end for	
	  MESSAGE(("----------------Finished initializing the Top y Particle Injectors-----------------")); 
	} // endif
	else {

      if (rank() == 0) MESSAGE(("----------------Reading the Particle Injectors-----------------")); 
      READ_INJECTOR(top_y, nx, nz, 0);
	}
      } //end top_y boundary

      if (global->bot_y) {
	      if (rank() == 0) MESSAGE(("----------------Initializing the Bot_y Particle Injectors-----------------")); 
	DEFINE_INJECTOR(bot_y,nx,nz);
	if (step() == 0) { 
	  for ( int n=1; n<=nsp; n++ ) { 
	    for ( int k=1;k<=nz; k++ ) {
	      for ( int j=1;j<=nx; j++ ) { 
		bbot_y(n,k,j) = 0;
		nbot_y(n,k,j) = npright(n)/nfac;
		ubot_y(1,n,k,j) = global->ur;
		ubot_y(2,n,k,j) = 0;
		ubot_y(3,n,k,j) = 0;
		pbot_y(1,2,n,k,j)=pbot_y(2,1,n,k,j)=pbot_y(1,3,n,k,j)=pbot_y(3,1,n,k,j)=pbot_y(2,3,n,k,j)=pbot_y(3,2,n,k,j)=0;
		pbot_y(1,1,n,k,j) = npright(n)*vth(n)*vth(n)/(2.0*nfac);
		pbot_y(2,2,n,k,j) = pbot_y(1,1,n,k,j);
		pbot_y(3,3,n,k,j) = pbot_y(1,1,n,k,j);
	      }      
	    }
	  }  // end for	
	} // endif
	else {

      if (rank() == 0) MESSAGE(("----------------Reading the Particle Injectors-----------------")); 
      READ_INJECTOR(bot_y, nx, nz, 0);
	}
      } //end bot_y boundary

      
      if (rank() == 0) MESSAGE(("-------------------------------------------------------------------")); 

    }// End of Intialization

    //            MESSAGE(("n=%i   npleft=%g    nfac=%g     vth=%g ",n,npleft(n),nfac(n),vth(n))); 

    //  Inject particles on Left Boundary

    if (global->left) {
      for ( int n=1; n<=nsp; n++ ) { 
	//	MESSAGE((" Injecting Left  --> n= %i    q=%e    vth=%e   nfac=%e",n,q(n),vth(n),nfac)); 
	species_t * species = find_species_id(n-1,species_list );  
	for ( int k=1;k<=nz; k++ ) {
    	  for ( int j=1;j<=ny; j++ ) { 
    	    if ( nleft(n,k,j) > 0.0 ) {
    	      vtherm = sqrt(2.0*pleft(1,1,n,k,j)/nleft(n,k,j));
    	      vd =  uleft(1,n,k,j)/vtherm;
    	      bleft(n,k,j) = bleft(n,k,j) + dt*nleft(n,k,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hx);
    	      inject = (int) bleft(n,k,j);
    	      bleft(n,k,j) = bleft(n,k,j) - (double) inject;
    	      double uflow[3] = {uleft(1,n,k,j),uleft(2,n,k,j),uleft(3,n,k,j)};
    	      double press[9] = {pleft(1,1,n,k,j),pleft(1,2,n,k,j),pleft(1,3,n,k,j),pleft(2,1,n,k,j),pleft(2,2,n,k,j),pleft(2,3,n,k,j),pleft(3,1,n,k,j),pleft(3,2,n,k,j),pleft(3,3,n,k,j)};	     
    	      //MESSAGE((" Injecting left  --> n= %i    inject=%i   nleft=%e    vth=%e  vd=%e",n,inject,nleft(n,k,j),vtherm,vd)); 
    	      repeat(inject) {	      
    		compute_injection(uv,nleft(n,k,j),uflow,press,1,2,3,rng(0));
    		x = grid->x0; 
		y = grid->y0 + hy*(j-1) + hy*uniform(rng(0), 0, 1); 
		z = grid->z0 + hz*(k-1) + hz*uniform(rng(0), 0, 1); 	    
    		age = 0;
		inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n) , age, 0 );   
    	      }
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
	    vtherm = sqrt(2.0*pright(1,1,n,k,j)/nright(n,k,j));
	    vd = -(global->ur)/vtherm;
	    bright(n,k,j) = bright(n,k,j)+ dt*nright(n,k,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hx);
	    inject = (int) bright(n,k,j);
	    bright(n,k,j) = bright(n,k,j) - (double) inject;
	    double uflow[3] = {uright(1,n,k,j),uright(2,n,k,j),uright(3,n,k,j)};
	    double press[9] = {pright(1,1,n,k,j),pright(1,2,n,k,j),pright(1,3,n,k,j),pright(2,1,n,k,j),pright(2,2,n,k,j),pright(2,3,n,k,j),pright(3,1,n,k,j),pright(3,2,n,k,j),pright(3,3,n,k,j)};	     

	    // MESSAGE((" Injecting right  --> n= %i    inject=%i   nright=%e    vth=%e  vd=%e",n,inject,nright(n,k,j),vtherm,vd)); 
	      // MESSAGE((" Injecting right  --> n= %i    inject=%i",n,inject)); 
	    repeat(inject) {
	      compute_injection(uv,nright(n,k,j),uflow,press,-1,2,3,rng(0));
	      //   MESSAGE((" Injecting right  --> n= %i    uvx=%e",n,uv[0])); 
	      x = grid->x1; 
	      y = grid->y0 + hy*(j-1) + hy*uniform(rng(0), 0, 1); 
	      z = grid->z0 + hz*(k-1) + hz*uniform(rng(0), 0, 1); 	    
	      age = 0;
	      //	      inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n) , age, 0 );
	    }
	  }
	}
      }
    } // end right injector
    
    // Inject particles on the top boundary
    if (global->top_y) {
      for ( int n=1; n<=nsp; n++ ) { 
	species_t * species = find_species_id(n-1,species_list );  
	for ( int k=1;k<=nz; k++ ) {
	  for ( int j=1;j<=nx; j++ ) {
	    vtherm = sqrt(2.0*ptop_y(1,1,n,k,j)/ntop_y(n,k,j));
	    vd = 0.0; // should be zero inflow
	    btop_y(n,k,j) = btop_y(n,k,j)+ dt*ntop_y(n,k,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hy);
	    inject = (int) btop_y(n,k,j);
	    btop_y(n,k,j) = btop_y(n,k,j) - (double) inject;
	    double uflow[3] = {utop_y(1,n,k,j),utop_y(2,n,k,j),utop_y(3,n,k,j)};
	    double press[9] = {ptop_y(1,1,n,k,j),ptop_y(1,2,n,k,j),ptop_y(1,3,n,k,j),ptop_y(2,1,n,k,j),ptop_y(2,2,n,k,j),ptop_y(2,3,n,k,j),ptop_y(3,1,n,k,j),ptop_y(3,2,n,k,j),ptop_y(3,3,n,k,j)};	     

	    // MESSAGE((" Injecting top_y  --> n= %i    inject=%i   ntop_y=%e    vth=%e  vd=%e",n,inject,ntop_y(n,k,j),vtherm,vd)); 
	      // MESSAGE((" Injecting top_y  --> n= %i    inject=%i",n,inject)); 
	    repeat(inject) {
	      compute_injection(uv,ntop_y(n,k,j),uflow,press,-2,1,3,rng(0));
	      //   MESSAGE((" Injecting top_y  --> n= %i    uvx=%e",n,uv[0])); 
	      x = grid->x0 + hx*(j-1) + hx*uniform(rng(0), 0, 1); 
	      y = grid->y1; 
	      z = grid->z0 + hz*(k-1) + hz*uniform(rng(0), 0, 1); 	    
	      age = 0;
	      inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n) , age, 0 );
	    }
	  }
	}
      }
    } // end top_y injector

    // Inject particles on the bottom_y boundary
    if (global->bot_y) {
      for ( int n=1; n<=nsp; n++ ) { 
	species_t * species = find_species_id(n-1,species_list );  
	for ( int k=1;k<=nz; k++ ) {
	  for ( int j=1;j<=nx; j++ ) {
	    vtherm = sqrt(2.0*pbot_y(1,1,n,k,j)/nbot_y(n,k,j));
	    vd = 0.0; // should be zero inflow
	    bbot_y(n,k,j) = bbot_y(n,k,j)+ dt*nbot_y(n,k,j)*vtherm*(exp(-vd*vd)/sqpi+vd*(erf(vd)+1))/(2*hy);
	    inject = (int) bbot_y(n,k,j);
	    bbot_y(n,k,j) = bbot_y(n,k,j) - (double) inject;
	    double uflow[3] = {ubot_y(1,n,k,j),ubot_y(2,n,k,j),ubot_y(3,n,k,j)};
	    double press[9] = {pbot_y(1,1,n,k,j),pbot_y(1,2,n,k,j),pbot_y(1,3,n,k,j),pbot_y(2,1,n,k,j),pbot_y(2,2,n,k,j),pbot_y(2,3,n,k,j),pbot_y(3,1,n,k,j),pbot_y(3,2,n,k,j),pbot_y(3,3,n,k,j)};	     

	    // MESSAGE((" Injecting bot_y  --> n= %i    inject=%i   nbot_y=%e    vth=%e  vd=%e",n,inject,nbot_y(n,k,j),vtherm,vd)); 
	      // MESSAGE((" Injecting bot_y  --> n= %i    inject=%i",n,inject)); 
	    repeat(inject) {
	      compute_injection(uv,nbot_y(n,k,j),uflow,press,2,3,1,rng(0));
	      //   MESSAGE((" Injecting bot_y  --> n= %i    uvx=%e",n,uv[0])); 
	      x = grid->x0 + hx*(j-1) + hx*uniform(rng(0), 0, 1); 
	      y = grid->y0; 
	      z = grid->z0 + hz*(k-1) + hz*uniform(rng(0), 0, 1); 	    
	      age = 0;
	      inject_particle(species, x, y, z, uv[0], uv[1], uv[2], q(n) , age, 0 );
	    }
	  }
	}
      }
    } // end bot_y injector
    
  // *******  Update the injector moments at every sort interval *********
  double v[3];
  double u[3];
  double p[9];

#define icell(i,j,k) INDEX_FORTRAN_3(i,j,k,0,nx+1,0,ny+1,0,nz+1)
#define v(i) v[INDEX_FORTRAN_1(i,1,3)]
#define u(i) u[INDEX_FORTRAN_1(i,1,3)]
#define p(i,j) p[INDEX_FORTRAN_2(i,j,1,3,1,3)]

  // Parameters for measuring moments  - BE CAREFUL -
  // don't make too big or we will go off the node

  int noff = 1;   // Offset from edge - to measure moments
                  // noff = 0  --> start with cell directly on boundary

  int nav = 2;  // How many cells to include in the "inward" direction in avg
  int navin = 2;// How many cells to include in the "inward" direction in avg


  // Right boundary Moments
  if (global->right) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      particle_t * part;
      if (remainder(step(), global->sort[n-1]) == 0) {
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
                double gamma = sqrt(1.0+part->ux*part->ux + part->uy*part->uy
                                    + part->uz*part->uz);
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
              fright(n,k,j) = (1.0-rout[1])*fright(n,k,j)
                + rout[1]*flux/(nav*hx);
              nright(n,k,j) = (1.0-rout[0])*nright(n,k,j) + rout[0]*npart/nav;
              for ( int a=1; a<=3; a++ ) {
                uright(a,n,k,j) = (1-rout[1])*uright(a,n,k,j)
                  + rout[1]*u(a)/npart;
                for ( int b=1;b<=3; b++ ) {
                  p(a,b) = (p(a,b) - u(a)*u(b)/npart)/nav;
                  pright(a,b,n,k,j) = (1-rout[2])*pright(a,b,n,k,j)
                    + rout[2]*p(a,b);
                }
              }
            }
          }
        }
      }
    }
  }  // end right moment update


  
  // top_y boundary Moments
  if (global->top_y) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      particle_t * part;
      if (remainder(step(), global->sort[n-1]) == 0) {
        double npart;
        for ( int k=1;k<=nz; k++ ) {
          for ( int i=1;i<=nx; i++ ) {
            npart = 0;
            flux = 0;
            u[0] = u[1] = u[2] = 0;
            p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
            for ( int j=ny-noff; j>ny-nav-noff; j-- ){
              int nstart = species->partition[icell(i,j,k)];
              int nstop  = species->partition[icell(i,j,k)+1];
              int ncell  = nstop - nstart;
              npart = npart + ncell;
              for (int np=nstart; np<nstop ; np++) {
                part=&species->p[np];
                double gamma = sqrt(1.0+part->ux*part->ux + part->uy*part->uy
                                    + part->uz*part->uz);
                v(1) = part->ux;
                v(2) = part->uy;
                v(3) = part->uz;
                if (v(2) < 0) flux = flux - v(2)/gamma;
                for ( int a=1;a<=3; a++ ) {
                  u(a) = u(a) + v(a);
                  for ( int b=1;b<=3; b++ ) {
                    p(a,b) = p(a,b) + v(a)*v(b);
                  }
                }
              } // end particle loop for single cell
            } // end cells included for these moments

            if ( npart > 0 ) {
              ftop_y(n,k,i) = (1.0-rout[1])*ftop_y(n,k,i)
                + rout[1]*flux/(nav*hy);
              ntop_y(n,k,i) = (1.0-rout[0])*ntop_y(n,k,i) + rout[0]*npart/nav;
              for ( int a=1; a<=3; a++ ) {
                utop_y(a,n,k,i) = (1-rout[1])*utop_y(a,n,k,i)
                  + rout[1]*u(a)/npart;
                for ( int b=1;b<=3; b++ ) {
                  p(a,b) = (p(a,b) - u(a)*u(b)/npart)/nav;
                  ptop_y(a,b,n,k,i) = (1-rout[2])*ptop_y(a,b,n,k,i)
                    + rout[2]*p(a,b);
                }
              }
            }
          }
        }
      }
    }
  }  // end top_y moment update

  // bot_y boundary Moments
  if (global->bot_y) {
    for ( int n=1; n<=nsp; n++ ) {
      species_t * species = find_species_id(n-1,species_list );
      particle_t * part;
      if (remainder(step(), global->sort[n-1]) == 0) {
        double npart;
        for ( int k=1;k<=nz; k++ ) {
          for ( int i=1;i<=nx; i++ ) {
            npart = 0;
            flux = 0;
            u[0] = u[1] = u[2] = 0;
            p[0] = p[1] = p[2] = p[3]= p[4] = p[5] = p[6]= p[7] = p[8] = 0;
            for ( int j=1+noff; j<nav+noff; j++ ){
              int nstart = species->partition[icell(i,j,k)];
              int nstop  = species->partition[icell(i,j,k)+1];
              int ncell  = nstop - nstart;
              npart = npart + ncell;
              for (int np=nstart; np<nstop ; np++) {
                part=&species->p[np];
                double gamma = sqrt(1.0+part->ux*part->ux + part->uy*part->uy
                                    + part->uz*part->uz);
                v(1) = part->ux;
                v(2) = part->uy;
                v(3) = part->uz;
                if (v(2) > 0) flux = flux + v(2)/gamma;
                for ( int a=1;a<=3; a++ ) {
                  u(a) = u(a) + v(a);
                  for ( int b=1;b<=3; b++ ) {
                    p(a,b) = p(a,b) + v(a)*v(b);
                  }
                }
              } // end particle loop for single cell
            } // end cells included for these moments

            if ( npart > 0 ) {
              fbot_y(n,k,i) = (1.0-rout[1])*fbot_y(n,k,i)
                + rout[1]*flux/(nav*hy);
              nbot_y(n,k,i) = (1.0-rout[0])*nbot_y(n,k,i) + rout[0]*npart/nav;
              for ( int a=1; a<=3; a++ ) {
                ubot_y(a,n,k,i) = (1-rout[1])*ubot_y(a,n,k,i)
                  + rout[1]*u(a)/npart;
                for ( int b=1;b<=3; b++ ) {
                  p(a,b) = (p(a,b) - u(a)*u(b)/npart)/nav;
                  pbot_y(a,b,n,k,i) = (1-rout[2])*pbot_y(a,b,n,k,i)
                    + rout[2]*p(a,b);
                }
              }
            }
          }
        }
      }
    }
  }  // end bot_y moment update


} // end particle injection

//   *******************  CURRENT INJECTION ***************************                                                                                                       

//  This injects current into a given region                                                                                                                                  


//  This injects current into a given region                                                                                                                                  

/*#define set_region_current(rgn,					\
			   eqn_ex,eqn_ey,eqn_ez) BEGIN_PRIMITIVE {        \
  const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;          \
  const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;          \
  const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;          \
  field_t *_f0 = field;                                                 \
  for( int _k=0; _k<=_nz+1; _k++ ) {                                    \
    const double _zl = _z0 + _dz*(_k-1.5), _ze = _z0 + _dz*_k, _zc = _z0 + _dz*(_k-0.5); \
    for( int _j=0; _j<=_ny+1; _j++ ) {                                  \
      const double _yl = _y0 + _dy*(_j-1.5), _ye = _y0 + _dy*_j, _yc = _y0 + _dy*(_j-0.5); \
      field_t *_f = _f0 + _LOCAL_CELL_ID(0,_j,_k);                      \
      for( int _i=0; _i<=_nx+1; _i++ ) {                                \
        const double _xl = _x0 + _dx*(_i-1.5), _xe = _x0 + _dx*_i, _xc = _x0 + _dx*(_i-0.5); \
        int _rccc, _rlcc, _rclc, _rllc, _rccl, _rlcl, _rcll, _rlll;     \
        double x, y, z;                                                 \
        x = _xc; y = _yc; z = _zc; _rccc = (rgn);                       \
        x = _xl;                   _rlcc = (rgn);                       \
        x = _xc; y = _yl;          _rclc = (rgn);                       \
        x = _xl;                   _rllc = (rgn);                       \
        x = _xc; y = _yc; z = _zl; _rccl = (rgn);                       \
        x = _xl;                   _rlcl = (rgn);                       \
        x = _xc; y = _yl;          _rcll = (rgn);                       \
        x = _xl;                   _rlll = (rgn);                       \
        x = _xc; y = _ye; z = _ze; if( _rccc || _rclc || _rccl || _rcll ) _f->jfx  = (eqn_ex); \
        x = _xe; y = _yc; z = _ze; if( _rccc || _rccl || _rlcc || _rlcl ) _f->jfy  = (eqn_ey); \
        x = _xe; y = _ye; z = _zc; if( _rccc || _rlcc || _rclc || _rllc ) _f->jfz  = (eqn_ez); \
        x = x;   y = y;   z = z;                                        \
        _f++;                                                           \
      }                                                                 \
    }                                                                   \
  }                                                                     \
  } END_PRIMITIVE */


// NOTE: Current injection adds current and does not set current. The reason                                                                                                        
//       for this is that user_current_injection is called _after_ accumulators                                                                                                     
//       are unloaded. If we overwrite particle currents by accident (e.g., on                                                                                                      
//       boundaries) then we will violate charge conservation.                                                                                                                      
#define set_region_current( rgn, eqn_jx, eqn_jy, eqn_jz ) do {          \
  const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;        \
  const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;        \
  const double _c  = grid->cvac;                                      \
  const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;        \
  for( int _k=0; _k<_nz+2; _k++ ) {                                   \
  const double _zl = _z0 + _dz*(_k-1.5);                            \
  const double _ze = _z0 + _dz*(_k-1.0);                            \
  const double _zc = _z0 + _dz*(_k-0.5);                            \
  for( int _j=0; _j<_ny+2; _j++ ) {                                 \
    const double _yl = _y0 + _dy*(_j-1.5);                          \
    const double _ye = _y0 + _dy*(_j-1.0);                          \
    const double _yc = _y0 + _dy*(_j-0.5);                          \
    field_t *_f = &field(0,_j,_k);                                  \
    for( int _i=0; _i<_nx+2; _i++ ) {                               \
      const double _xl = _x0 + _dx*(_i-1.5);                        \
      const double _xe = _x0 + _dx*(_i-1.0);                        \
      const double _xc = _x0 + _dx*(_i-0.5);                        \
      double x, y, z;                                               \
      int _rccc, _rlcc, _rclc, _rllc, _rccl, _rlcl, _rcll;          \
      x = _xc; y = _yc; z = _zc; _rccc = (rgn);                     \
      x = _xl;                   _rlcc = (rgn);                     \
      x = _xc; y = _yl;          _rclc = (rgn);                     \
      x = _xl;                   _rllc = (rgn);                     \
      x = _xc; y = _yc; z = _zl; _rccl = (rgn);                     \
      x = _xl;                   _rlcl = (rgn);                     \
      x = _xc; y = _yl;          _rcll = (rgn);                     \
      x = _xc; y = _ye; z = _ze; if( _rccc || _rclc || _rccl || _rcll ) _f->jfx += (eqn_jx); \
      x = _xe; y = _yc; z = _ze; if( _rccc || _rccl || _rlcc || _rlcl ) _f->jfy += (eqn_jy); \
      x = _xe; y = _ye; z = _zc; if( _rccc || _rlcc || _rclc || _rllc ) _f->jfz += (eqn_jz); \
      _f++;                                                         \
    }}}                                                             \
       } while(0)

begin_current_injection {

  // inject the coil current                                                                                                                                                  

  double t = grid->dt*step();

  if (t < global->t_dipole_start ) {

    double jt = 0;
    set_region_current(top_wire, 0, 0, jt);
    set_region_current(bot_wire, 0, 0, -jt);

  } else {

    double jt = global->c_j0*(1.0-exp(-(t-global->t_dipole_start)/(global->t_dipole_rise)));                                                                                                                  

    set_region_current(top_wire, 0, 0, jt);
    set_region_current(bot_wire, 0, 0, -jt);

  }
}


begin_field_injection {

  // No field injection for this simulation


  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;
  int x,y,z;
  double b0 = global->b0;
  double sn = global->sn;
  double cs = global->cs;
  double Vflow = global->ur;
  // There macros are from local.c to apply boundary conditions
#define XYZ_LOOP(xl,xh,yl,yh,zl,zh)             \
  for( z=zl; z<=zh; z++ )                       \
    for( y=yl; y<=yh; y++ )                     \
      for( x=xl; x<=xh; x++ )

#define yz_EDGE_LOOP(x) XYZ_LOOP(x,x,1,ny+1,1,1+nz)
#define xz_EDGE_LOOP(y) XYZ_LOOP(1,nx+1,y,y,1,1+nz)

  // Right Boundary
  //  if (global->right) {
  //    yz_EDGE_LOOP(nx+1) field(x,y,z).ez = Vflow*-b0*sn;
  //    yz_EDGE_LOOP(nx+1) field(x,y,z).cby = b0*sn; 
  //  }
  // no need to inject field on the right
  if (global->left) {
    yz_EDGE_LOOP(1) field(x,y,z).ez = Vflow*-b0*sn;
    //    yz_EDGE_LOOP(1) field(x,y,z).cby = b0*sn; 
    yz_EDGE_LOOP(1) field(x,y,z).cbx = b0*cs;
    // should inject bx as well?
  }

  // we have the weird case of a flowing boundary so we need 
  // to set ez and B on the top and bottom boundaries.
  // do we need to set Bx as well?
  if (global->bot_y) {
    xz_EDGE_LOOP(1) field(x,y,z).ez = Vflow*-b0*sn;
    xz_EDGE_LOOP(1) field(x,y,z).cby = b0*sn;
    //    xz_EDGE_LOOP(1) field(x,y,z).cbx = b0*cs;
  }

  if (global->top_y) {
    xz_EDGE_LOOP(ny+1) field(x,y,z).ez = Vflow*-b0*sn;
    xz_EDGE_LOOP(ny+1) field(x,y,z).cby = b0*sn;
    //    xz_EDGE_LOOP(ny) field(x,y,z).cbx = b0*cs;
    //    xz_EDGE_LOOP(ny+1) field(x,y,z).cbx = b0*cs;
    }


}  // end field injection


begin_particle_collisions {

  // No particle collisions in this simulation

}
