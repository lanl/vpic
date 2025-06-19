/*
 * This piece of code computes time-averaged quantites
 *
 * Global variables:
 *
 * global->dis_nav  : number of steps to average over                             : initialized here
 * global->dis_interval  : number of steps between outputs.                       : initialized here
 * global->dis_iter : iteration count. 0 means we are not avraging at the moment  : initialized here
 * global->dis_begin_int : the first time step of the interval                   : initialized here
 *
 * The adjustable parameters of these, dis_interval and dis_nav, are now set by the variables
 * AVG_SPACING and AVG_TOTAL_STEPS in the beginning so they can be more easily found and edited
 *
 * The data is saved in the interval 
 *
 * j*dis_interval <= step < j*dis_interval + dis_nav
 *
 * where j is an integer. the output is assigned time index corresponding to the start of the interval
 *
 * The code is restart-aware. The associated global variables are
 *
 * global->write_restart ;                 // global flag for all to write restart files
 *
 */


// fortran-like indices for arrays (convenience)
#define ALLOCATE(A,LEN,TYPE)                                             \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate."));


#define LOCAL_CELL_ID(x,y,z) \
  INDEX_FORTRAN_3(x,y,z,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1)


#define hydro_array(x,y,z) hydro_array->h[LOCAL_CELL_ID(x,y,z)]
#define Ex(i,j,k)   Ex[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Ey(i,j,k)   Ey[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Ez(i,j,k)   Ez[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Bx(i,j,k)   Bx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define By(i,j,k)   By[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Bz(i,j,k)   Bz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Jex(i,j,k)  Jex[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Jey(i,j,k)  Jey[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Jez(i,j,k)  Jez[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
//#define Jfx(i,j,k)  Jfx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
//#define Jfy(i,j,k)  Jfy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
//#define Jfz(i,j,k)  Jfz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define Jix(i,j,k)  Jix[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Jiy(i,j,k)  Jiy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Jiz(i,j,k)  Jiz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define Pe_xx(i,j,k)  Pe_xx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pe_yy(i,j,k)  Pe_yy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pe_zz(i,j,k)  Pe_zz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pe_xy(i,j,k)  Pe_xy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pe_xz(i,j,k)  Pe_xz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pe_yz(i,j,k)  Pe_yz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define Pi_xx(i,j,k)  Pi_xx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pi_yy(i,j,k)  Pi_yy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pi_zz(i,j,k)  Pi_zz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pi_xy(i,j,k)  Pi_xy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pi_xz(i,j,k)  Pi_xz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Pi_yz(i,j,k)  Pi_yz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define ne(i,j,k)  ne[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define ni(i,j,k)  ni[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

// Ohm's law terms
#define neEx(i,j,k)  neEx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]//E-field
#define neEy(i,j,k)  neEy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define neEz(i,j,k)  neEz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define edotj(i,j,k)  edotj[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define jexbx(i,j,k)  jexbx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)] //je x b
#define jexby(i,j,k)  jexby[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define jexbz(i,j,k)  jexbz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define neUUe_xx(i,j,k)  neUUe_xx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)] //Reynolds
#define neUUe_yy(i,j,k)  neUUe_yy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define neUUe_zz(i,j,k)  neUUe_zz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define neUUe_xy(i,j,k)  neUUe_xy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define neUUe_xz(i,j,k)  neUUe_xz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define neUUe_yz(i,j,k)  neUUe_yz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define djedtx(i,j,k)  djedtx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define djedty(i,j,k)  djedty[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define djedtz(i,j,k)  djedtz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define dBdtx(i,j,k)  dBdtx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define dBdty(i,j,k)  dBdty[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define dBdtz(i,j,k)  dBdtz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

//#define mix_e(i,j,k)  mix_e[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
//#define mix_i(i,j,k)  mix_i[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

//Ion force balance terms
#define niEx(i,j,k)  niEx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]//E-field
#define niEy(i,j,k)  niEy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define niEz(i,j,k)  niEz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define jixbx(i,j,k)  jixbx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)] //je x b
#define jixby(i,j,k)  jixby[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define jixbz(i,j,k)  jixbz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define niUUi_xx(i,j,k)  niUUi_xx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)] //Reynolds
#define niUUi_yy(i,j,k)  niUUi_yy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define niUUi_zz(i,j,k)  niUUi_zz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define niUUi_xy(i,j,k)  niUUi_xy[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define niUUi_xz(i,j,k)  niUUi_xz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define niUUi_yz(i,j,k)  niUUi_yz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define djidtx(i,j,k)  djidtx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define djidty(i,j,k)  djidty[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define djidtz(i,j,k)  djidtz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

/*#define Te(i,j,k) Te[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]

#define Fcollx(i,j,k) Fcollx[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Fcolly(i,j,k) Fcolly[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]
#define Fcollz(i,j,k) Fcollz[INDEX_FORTRAN_3(i,j,k,1,nx,1,ny,1,nz)]*/


#define HydroToCCijk(term) (hydro_array(i,j,k).term + hydro_array(i,j,k+1).term + hydro_array(i,j+1,k).term + hydro_array(i,j+1,k+1).term + hydro_array(i+1,j,k).term + hydro_array(i+1,j,k+1).term + hydro_array(i+1,j+1,k).term + hydro_array(i+1,j+1,k+1).term)*0.125

// vpic definition:
//  float txx, tyy, tzz;   // Stress diagonal            => <p_i v_j f>, i==j
//  float tyz, tzx, txy;   // Stress off-diagonal        => <p_i v_j f>, i!=j


// loop over all cells inside domain (exclude ghost cells)
#define for_all_cells							\
  for(k=1; k<=nz; k++)							\
    for(j=1; j<=ny; j++)						\
      for(i=1; i<=nx; i++)						

{ // start the diagnostics

  static int dis_initialized = 0;   // the flag that signals intialization
  static float *Ex,*Ey,*Ez,*Bx,*By,*Bz;  // current and e*j arrays
  static float *Pe_xx,*Pe_yy,*Pe_zz,*Pe_xy,*Pe_xz,*Pe_yz;  // Pe
  static float *Pi_xx,*Pi_yy,*Pi_zz,*Pi_xy,*Pi_xz,*Pi_yz;  // Pi
  static float *Jex,*Jey,*Jez;
  //  static float *Jfx,*Jfy,*Jfz;
  static float *Jix,*Jiy,*Jiz;
  //  static float *mix_e,*mix_i;
  static float *ne, *ni;  // density

  //Ohm's law terms
  //(pressure is already being recorded)
  static float *neEx, *neEy, *neEz, *edotj ; //n_e E
  static float *jexbx,*jexby,*jexbz; //j_e x b
  static float *neUUe_xx, *neUUe_yy, *neUUe_zz, *neUUe_xy, *neUUe_xz, *neUUe_yz; //Reynolds
  static float *djedtx, *djedty, *djedtz; //inertia
  static float *dBdtx, *dBdty, *dBdtz; //B field change

  static float *niEx, *niEy, *niEz; //n_i E
  static float *jixbx,*jixby,*jixbz; //j_i x b
  static float *niUUi_xx, *niUUi_yy, *niUUi_zz, *niUUi_xy, *niUUi_xz, *niUUi_yz; //Ion Reynolds
  static float *djidtx, *djidty, *djidtz; //ion inertia

  /*  static float *Te;
      static float *Fcollx, *Fcolly, *Fcollz;*/

  int i,j,k,isp;
  FileIO fp;
  char fname [256];
    
  const int nx = grid->nx;
  const int ny = grid->ny;
  const int nz = grid->nz;

  static int dis_sample;
  
  //PARAMETERS FOR AVERAGING INTERVAL
  const int AVG_SPACING = 1;//30;      //18; //Spacing between averagings. 1 should mean every time the 
  //program writes non-averaged data, 2 every other, etc.
  const int AVG_TOTAL_STEPS = 100;//60;  //Total number of steps to average over (NOT dt's)
  const int AVG_SAMPLE = 5;       // Sample frequency. Total average interval = AVG_SAMPLE*AVG_TOTAL_STEPS. 

  if (!dis_initialized) { // initialization upon start or re-start
    
    ALLOCATE(Ex,nx*ny*nz,float);
    ALLOCATE(Ey,nx*ny*nz,float);
    ALLOCATE(Ez,nx*ny*nz,float);

    ALLOCATE(Bx,nx*ny*nz,float);
    ALLOCATE(By,nx*ny*nz,float);
    ALLOCATE(Bz,nx*ny*nz,float);

    ALLOCATE(Jex,nx*ny*nz,float);
    ALLOCATE(Jey,nx*ny*nz,float);
    ALLOCATE(Jez,nx*ny*nz,float);

    //    ALLOCATE(Jfx,nx*ny*nz,float);
    //    ALLOCATE(Jfy,nx*ny*nz,float);
    //    ALLOCATE(Jfz,nx*ny*nz,float);

    ALLOCATE(Jix,nx*ny*nz,float);
    ALLOCATE(Jiy,nx*ny*nz,float);
    ALLOCATE(Jiz,nx*ny*nz,float);

    ALLOCATE(Pe_xx,nx*ny*nz,float);
    ALLOCATE(Pe_yy,nx*ny*nz,float);
    ALLOCATE(Pe_zz,nx*ny*nz,float);
    ALLOCATE(Pe_xy,nx*ny*nz,float);
    ALLOCATE(Pe_xz,nx*ny*nz,float);
    ALLOCATE(Pe_yz,nx*ny*nz,float);

    ALLOCATE(Pi_xx,nx*ny*nz,float);
    ALLOCATE(Pi_yy,nx*ny*nz,float);
    ALLOCATE(Pi_zz,nx*ny*nz,float);
    ALLOCATE(Pi_xy,nx*ny*nz,float);
    ALLOCATE(Pi_xz,nx*ny*nz,float);
    ALLOCATE(Pi_yz,nx*ny*nz,float);

    ALLOCATE(ne,nx*ny*nz,float);
    ALLOCATE(ni,nx*ny*nz,float);

    //    ALLOCATE(mix_e,nx*ny*nz,float);
    // ALLOCATE(mix_i,nx*ny*nz,float);

    //Ohm's law terms    
    ALLOCATE(neEx,nx*ny*nz,float); //E-field
    ALLOCATE(neEy,nx*ny*nz,float);
    ALLOCATE(neEz,nx*ny*nz,float);

    ALLOCATE(edotj,nx*ny*nz,float);
    
    ALLOCATE(jexbx,nx*ny*nz,float); //j_e x b
    ALLOCATE(jexby,nx*ny*nz,float);
    ALLOCATE(jexbz,nx*ny*nz,float);

    ALLOCATE(neUUe_xx,nx*ny*nz,float); //Reynolds
    ALLOCATE(neUUe_yy,nx*ny*nz,float);
    ALLOCATE(neUUe_zz,nx*ny*nz,float);
    ALLOCATE(neUUe_xy,nx*ny*nz,float);
    ALLOCATE(neUUe_xz,nx*ny*nz,float);
    ALLOCATE(neUUe_yz,nx*ny*nz,float);

    ALLOCATE(djedtx,nx*ny*nz,float); //inertia
    ALLOCATE(djedty,nx*ny*nz,float);
    ALLOCATE(djedtz,nx*ny*nz,float);

    ALLOCATE(dBdtx,nx*ny*nz,float); //B field change
    ALLOCATE(dBdty,nx*ny*nz,float);
    ALLOCATE(dBdtz,nx*ny*nz,float);

    //Ion force terms    
    ALLOCATE(niEx,nx*ny*nz,float); //E-field
    ALLOCATE(niEy,nx*ny*nz,float);
    ALLOCATE(niEz,nx*ny*nz,float);
    
    ALLOCATE(jixbx,nx*ny*nz,float); //j_i x b
    ALLOCATE(jixby,nx*ny*nz,float);
    ALLOCATE(jixbz,nx*ny*nz,float);

    ALLOCATE(niUUi_xx,nx*ny*nz,float); //Ion Reynolds
    ALLOCATE(niUUi_yy,nx*ny*nz,float);
    ALLOCATE(niUUi_zz,nx*ny*nz,float);
    ALLOCATE(niUUi_xy,nx*ny*nz,float);
    ALLOCATE(niUUi_xz,nx*ny*nz,float);
    ALLOCATE(niUUi_yz,nx*ny*nz,float);

    ALLOCATE(djidtx,nx*ny*nz,float); //ion inertia
    ALLOCATE(djidty,nx*ny*nz,float);
    ALLOCATE(djidtz,nx*ny*nz,float);

    /*    ALLOCATE(Te,nx*ny*nz,float);

    ALLOCATE(Fcollx,nx*ny*nz,float);
    ALLOCATE(Fcolly,nx*ny*nz,float);
    ALLOCATE(Fcollz,nx*ny*nz,float);*/

    dump_mkdir("hydro-ave");
    sprintf(fname,"hydro-ave/%d",NUMFOLD);
    dump_mkdir(fname);


    if (step()==0) {

      //      dump_mkdir("hydro-ave"); // MOVE ABOVE: NOT SURE FOLDERS WILL BE THERE AFTER STAGE-OUT STAGE-IN.
      // sprintf(fname,"hydro-ave/%d",NUMFOLD);
      // dump_mkdir(fname);

      global->dis_iter = 0;      // initialize iteration count
      global->dis_begin_int = 0; // initialize the start of the averaging interval
    
    }

      //These are labelled at the top
      global->dis_nav = AVG_TOTAL_STEPS;   
      global->dis_interval = AVG_SPACING*global->fields_interval;
      dis_sample = AVG_SAMPLE;

      if ((global->dis_begin_int) < (global->dis_interval)) {
	global->dis_begin_int = global->dis_interval;  // Skip t=0 averaging. 
      }

    // do we need to restart
    if (global->dis_iter>0) {

      //      MESSAGE(("-----dis_iter = %d  %d\n",global->dis_iter,global->dis_begin_int)); 

      sprintf(fname,"hydro-ave/%d/T.%d/%s.%d.%d",NUMFOLD,global->dis_begin_int,"data-ave",global->dis_begin_int,(int) rank());
      if ( ! (fp.open(fname, io_read)==ok) ) ERROR(("Cannot open file with time-averaged data."));
      sim_log( "reading the in-place averaging diagnostics Bx" );
      fp.read(Bx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics By" );
      fp.read(By, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Bz" );
      fp.read(Bz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Ex" );
      fp.read(Ex, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Ey" );
      fp.read(Ey, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Ez" );
      fp.read(Ez, nx*ny*nz );
      //      sim_log( "reading the in-place averaging diagnostics Jfx" );
      //      fp.read(Jfx, nx*ny*nz );
      //      sim_log( "reading the in-place averaging diagnostics Jfy" );
      //      fp.read(Jfy, nx*ny*nz );
      //      sim_log( "reading the in-place averaging diagnostics Jfz" );
      //      fp.read(Jfz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Jex" );
      fp.read(Jex, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Jey" );
      fp.read(Jey, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Jez" );
      fp.read(Jez, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Jix" );
      fp.read(Jix, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Jiy" );
      fp.read(Jiy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Jiz" );
      fp.read(Jiz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pe-xx" );
      fp.read(Pe_xx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pe-yy" );
      fp.read(Pe_yy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pe-zz" );
      fp.read(Pe_zz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pe-xy" );
      fp.read(Pe_xy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pe-xz" );
      fp.read(Pe_xz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pe-yz" );
      fp.read(Pe_yz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pi-xx" );
      fp.read(Pi_xx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pi-yy" );
      fp.read(Pi_yy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pi-zz" );
      fp.read(Pi_zz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pi-xy" );
      fp.read(Pi_xy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pi-xz" );
      fp.read(Pi_xz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :Pi-yz" );
      fp.read(Pi_yz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :ne" );
      fp.read(ne, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics :ni" );
      fp.read(ni, nx*ny*nz );
      //     sim_log( "reading the in-place averaging diagnostics :mix_e" );
      //      fp.read(mix_e, nx*ny*nz );
      //sim_log( "reading the in-place averaging diagnostics :mix_i" );
      //fp.read(mix_i, nx*ny*nz );

      //Ohm's law terms      
      sim_log( "reading the in-place averaging diagnostics ne*Ex" ); // n_e E
      fp.read(neEx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics ne*Ey" );
      fp.read(neEy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics ne*Ez" );
      fp.read(neEz, nx*ny*nz );

      sim_log( "reading the in-place averaging diagnostics E*J" );
      fp.read(edotj, nx*ny*nz );

      sim_log( "reading the in-place averaging diagnostics (je x B)_x" ); // j_e x B
      fp.read(jexbx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (je x B)_y" );
      fp.read(jexby, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (je x B)_z" );
      fp.read(jexbz, nx*ny*nz );

      sim_log( "reading the in-place averaging diagnostics (neUeUe)_xx" ); //Reynolds
      fp.read(neUUe_xx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (neUeUe)_yy" );
      fp.read(neUUe_yy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (neUeUe)_zz" );
      fp.read(neUUe_zz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (neUeUe)_xy" );
      fp.read(neUUe_xy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (neUeUe)_xz" );
      fp.read(neUUe_xz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (neUeUe)_yz" );
      fp.read(neUUe_yz, nx*ny*nz );

      sim_log( "reading the in-place averaging diagnostics d(n_e U_ex)/dt" ); //Inertia
      fp.read(djedtx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics d(n_e U_ey)/dt" );
      fp.read(djedty, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics d(n_e U_ez)/dt" );
      fp.read(djedtz, nx*ny*nz );

      sim_log( "reading the in-place averaging diagnostics d(B_x)/dt" ); //Inertia
      fp.read(dBdtx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics d(B_y)/dt" );
      fp.read(dBdty, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics d(B_z)/dt" );
      fp.read(dBdtz, nx*ny*nz );

      //Ion force terms
      sim_log( "reading the in-place averaging diagnostics ni*Ex" ); // n_i E
      fp.read(niEx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics ni*Ey" );
      fp.read(niEy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics ni*Ez" );
      fp.read(niEz, nx*ny*nz );

      sim_log( "reading the in-place averaging diagnostics (ji x B)_x" ); // j_i x B
      fp.read(jixbx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (ji x B)_y" );
      fp.read(jixby, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (ji x B)_z" );
      fp.read(jixbz, nx*ny*nz );

      sim_log( "reading the in-place averaging diagnostics (niUiUi)_xx" ); //Ion Reynolds
      fp.read(niUUi_xx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (niUiUi)_yy" );
      fp.read(niUUi_yy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (niUiUi)_zz" );
      fp.read(niUUi_zz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (niUiUi)_xy" );
      fp.read(niUUi_xy, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (niUiUi)_xz" );
      fp.read(niUUi_xz, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics (niUiUi)_yz" );
      fp.read(niUUi_yz, nx*ny*nz );

      sim_log( "reading the in-place averaging diagnostics d(n_i U_ix)/dt" ); //Ion Inertia
      fp.read(djidtx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics d(n_i U_iy)/dt" );
      fp.read(djidty, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics d(n_i U_iz)/dt" );
      fp.read(djidtz, nx*ny*nz );

      /*      sim_log( "reading the in-place averaging diagnostics Te" );
      fp.read(Te, nx*ny*nz );

      sim_log( "reading the in-place averaging diagnostics Fcollx" );
      fp.read(Fcollx, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Fcolly" );
      fp.read(Fcolly, nx*ny*nz );
      sim_log( "reading the in-place averaging diagnostics Fcollz" );
      fp.read(Fcollz, nx*ny*nz );*/


      fp.close();	    
    } // end of file read

    dis_initialized = 1;
	  
  } // end of initialization

  //  This way of computing remainder has issues for some reason - Bill
  //  int r=remainder(step, global->dis_interval);

  int r= (step()%global->dis_interval);
  int sample=r%dis_sample;

  if (( r >=0) && (r < global->dis_nav*dis_sample && step()>= global->dis_begin_int )) { // check that we are inside the averaging interval

    if ( sample == 0 ) {

    // if (rank() == 0) MESSAGE(("Adding Slice to Average ---> r = %i   step=%d\n",r,step));     

      if (global->dis_iter == 0) { //  this is the first teration
	for_all_cells Ex(i,j,k) = 0.0; // resets E 
	for_all_cells Ey(i,j,k) = 0.0; // resets  
	for_all_cells Ez(i,j,k) = 0.0; // resets  

	for_all_cells Bx(i,j,k) = 0.0; // resets B 
	for_all_cells By(i,j,k) = 0.0; // resets  
	for_all_cells Bz(i,j,k) = 0.0; // resets  

	//      for_all_cells Jfx(i,j,k) = 0.0; // resets Jfree
	//      for_all_cells Jfy(i,j,k) = 0.0; // resets  
	//      for_all_cells Jfz(i,j,k) = 0.0; // resets  
	
	for_all_cells Jex(i,j,k) = 0.0; // resets Je
	for_all_cells Jey(i,j,k) = 0.0; // resets  
	for_all_cells Jez(i,j,k) = 0.0; // resets  
	
	for_all_cells Jix(i,j,k) = 0.0; // resets Ji
	for_all_cells Jiy(i,j,k) = 0.0; // resets  
	for_all_cells Jiz(i,j,k) = 0.0; // resets  
	
	for_all_cells Pe_xx(i,j,k) = 0.0; // resets P
	for_all_cells Pe_yy(i,j,k) = 0.0; // resets P
	for_all_cells Pe_zz(i,j,k) = 0.0; // resets P
	for_all_cells Pe_xy(i,j,k) = 0.0; // resets P
	for_all_cells Pe_xz(i,j,k) = 0.0; // resets P
	for_all_cells Pe_yz(i,j,k) = 0.0; // resets P
	
	for_all_cells Pi_xx(i,j,k) = 0.0; // resets P
	for_all_cells Pi_yy(i,j,k) = 0.0; // resets P
	for_all_cells Pi_zz(i,j,k) = 0.0; // resets P
	for_all_cells Pi_xy(i,j,k) = 0.0; // resets P
	for_all_cells Pi_xz(i,j,k) = 0.0; // resets P
	for_all_cells Pi_yz(i,j,k) = 0.0; // resets P
	
	for_all_cells ne(i,j,k) = 0.0; // resets ne
	for_all_cells ni(i,j,k) = 0.0; // resets ni
	
	//      for_all_cells mix_e(i,j,k) = 0.0; // resets mix
	// for_all_cells mix_i(i,j,k) = 0.0; // 
	
	//Ohm's law terms
	
	for_all_cells neEx(i,j,k) = 0.0; // resets n_e*E 
	for_all_cells neEy(i,j,k) = 0.0; // |  
	for_all_cells neEz(i,j,k) = 0.0; // v
	
	for_all_cells edotj(i,j,k) = 0.0; // reset E*J
	
	for_all_cells jexbx(i,j,k) = 0.0; // resets j_e x B 
	for_all_cells jexby(i,j,k) = 0.0; // |  
	for_all_cells jexbz(i,j,k) = 0.0; // v
	
	for_all_cells neUUe_xx(i,j,k) = 0.0; // resets n_e U_e U_e
	for_all_cells neUUe_yy(i,j,k) = 0.0; // |
	for_all_cells neUUe_zz(i,j,k) = 0.0; // |
	for_all_cells neUUe_xy(i,j,k) = 0.0; // |
	for_all_cells neUUe_xz(i,j,k) = 0.0; // |
	for_all_cells neUUe_yz(i,j,k) = 0.0; // v
	
	//Ion Force Terms
	
	for_all_cells niEx(i,j,k) = 0.0; // resets n_i*E 
	for_all_cells niEy(i,j,k) = 0.0; // |  
	for_all_cells niEz(i,j,k) = 0.0; // v
	
	for_all_cells jixbx(i,j,k) = 0.0; // resets j_i x B 
	for_all_cells jixby(i,j,k) = 0.0; // |  
	for_all_cells jixbz(i,j,k) = 0.0; // v
	
	for_all_cells niUUi_xx(i,j,k) = 0.0; // resets n_i U_i U_i
	for_all_cells niUUi_yy(i,j,k) = 0.0; // |
	for_all_cells niUUi_zz(i,j,k) = 0.0; // |
	for_all_cells niUUi_xy(i,j,k) = 0.0; // |
	for_all_cells niUUi_xz(i,j,k) = 0.0; // |
	for_all_cells niUUi_yz(i,j,k) = 0.0; // v
	
	for_all_cells djedtx(i,j,k) = 0.0; // reset inertial term
	for_all_cells djedty(i,j,k) = 0.0; // reset inertial term
	for_all_cells djedtz(i,j,k) = 0.0; // reset inertial term

	/*        for_all_cells Te(i,j,k) = 0.0;     // reset electron temp. term

        for_all_cells Fcollx(i,j,k) = 0.0; // reset friction term
        for_all_cells Fcolly(i,j,k) = 0.0; // reset friction term
        for_all_cells Fcollz(i,j,k) = 0.0; // reset friction term*/

	
	//The d(n_e u_e)/dt term needs to be dealt with differently so is not included here
	
	global->dis_begin_int = step();    // mark the beginning of the interval
	sim_log( "Time averaging:start accumulation----dis_begin_int="<<global->dis_begin_int );
	
      }
    
      species_t * sp;
      
      sp = find_species_name("electron", species_list);
      clear_hydro_array(hydro_array);
      accumulate_hydro_p(hydro_array, sp, interpolator_array);
      synchronize_hydro_array(hydro_array);
      for_all_cells {
	// need to check for zeros here
	Jex(i,j,k) = Jex(i,j,k) + HydroToCCijk(jx);
      }
      for_all_cells {
	// need to check for zeros here
	Jey(i,j,k) = Jey(i,j,k) + HydroToCCijk(jy);
      }
      for_all_cells {
	// need to check for zeros here
	Jez(i,j,k) = Jez(i,j,k) + HydroToCCijk(jz);
      }
      
      for_all_cells {
	Pe_xx(i,j,k) = Pe_xx(i,j,k) + HydroToCCijk(txx);
	Pe_yy(i,j,k) = Pe_yy(i,j,k) + HydroToCCijk(tyy);
	Pe_zz(i,j,k) = Pe_zz(i,j,k) + HydroToCCijk(tzz);
	Pe_xy(i,j,k) = Pe_xy(i,j,k) + HydroToCCijk(txy);
	Pe_xz(i,j,k) = Pe_xz(i,j,k) + HydroToCCijk(tzx);
	Pe_yz(i,j,k) = Pe_yz(i,j,k) + HydroToCCijk(tyz);      
      }
      
      for_all_cells {
	ne(i,j,k) = ne(i,j,k) - HydroToCCijk(rho);
	//mix_e(i,j,k) = mix_e(i,j,k) - hydro(i,j,k).rho;
      }
      
      //Ohm's law terms
      
      //n_e*E
      for_all_cells {
	neEx(i,j,k) += HydroToCCijk(rho)*(field(i,j,k).ex + field(i,j+1,k).ex + field(i,j,k+1).ex + field(i,j+1,k+1).ex)*0.25;
	neEy(i,j,k) += HydroToCCijk(rho)*(field(i,j,k).ey + field(i+1,j,k).ey + field(i,j,k+1).ey + field(i+1,j,k+1).ey)*0.25;
	neEz(i,j,k) += HydroToCCijk(rho)*(field(i,j,k).ez + field(i+1,j,k).ez + field(i,j+1,k).ez + field(i+1,j+1,k).ez)*0.25;
	
	edotj(i,j,k) += (
			 field(i,j,k).ex*field(i,j,k).jfx + //x component centered
			 field(i,j,k+1).ex*field(i,j,k+1).jfx +
			 field(i,j+1,k).ex*field(i,j+1,k).jfx +
			 field(i,j+1,k+1).ex*field(i,j+1,k+1).jfx +
			 field(i,j,k).ey*field(i,j,k).jfy + //y component centered
			 field(i,j,k+1).ey*field(i,j,k+1).jfy +
			 field(i+1,j,k).ey*field(i+1,j,k).jfy +
			 field(i+1,j,k+1).ey*field(i+1,j,k+1).jfy + 
			 field(i,j,k).ez*field(i,j,k).jfz + //z component centered
			 field(i,j+1,k).ez*field(i,j+1,k).jfz +
			 field(i+1,j,k).ez*field(i+1,j,k).jfz +
			 field(i+1,j+1,k).ez*field(i+1,j+1,k).jfz
			 )*0.25; //Only compute when looping over the top electrons
      }
    


      // j_e x B
      for_all_cells {
	
	// Assume current in hydro is cell centered
	
	//jexbx(i,j,k) += hydro(i,j,k).jy*(field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5 - hydro(i,j,k).jz*(field(i,j,k).cby + field(i,j+1,k).cby)*0.5;
	//jexby(i,j,k) += hydro(i,j,k).jz*(field(i,j,k).cbx + field(i+1,j,k).cbx)*0.5 - hydro(i,j,k).jx*(field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5;
	//jexbz(i,j,k) += hydro(i,j,k).jx*(field(i,j,k).cby + field(i,j+1,k).cby)*0.5 - hydro(i,j,k).jy*(field(i,j,k).cbx + field(i+1,j,k).cbx)*0.5;
	
	//  Assume both B and J are on Yee Mesh
	
	jexbx(i,j,k) += (HydroToCCijk(jy))*(field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5;
	jexbx(i,j,k) -= (HydroToCCijk(jz))*(field(i,j,k).cby + field(i,j+1,k).cby)*0.5;
	
	jexby(i,j,k) += (HydroToCCijk(jz))*(field(i,j,k).cbx + field(i+1,j,k).cbx)*0.5;
	jexby(i,j,k) -= (HydroToCCijk(jx))*(field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5;
	
	jexbz(i,j,k) += (HydroToCCijk(jx))*(field(i,j,k).cby + field(i,j+1,k).cby)*0.5;
	jexbz(i,j,k) -= (HydroToCCijk(jy))*(field(i,j,k).cbx + field(i+1,j,k).cbx)*0.5;
	
      }   
      
      // ne U_e U_e
      for_all_cells{
	if (fabsf(HydroToCCijk(rho)) > 0) {
	  neUUe_xx(i,j,k) += HydroToCCijk(jx) * HydroToCCijk(jx) / fabsf(HydroToCCijk(rho));
	  neUUe_yy(i,j,k) += HydroToCCijk(jy) * HydroToCCijk(jy) / fabsf(HydroToCCijk(rho));
	  neUUe_zz(i,j,k) += HydroToCCijk(jz) * HydroToCCijk(jz) / fabsf(HydroToCCijk(rho));
	  neUUe_xy(i,j,k) += HydroToCCijk(jx) * HydroToCCijk(jy) / fabsf(HydroToCCijk(rho));
	  neUUe_xz(i,j,k) += HydroToCCijk(jx) * HydroToCCijk(jz) / fabsf(HydroToCCijk(rho));
	  neUUe_yz(i,j,k) += HydroToCCijk(jy) * HydroToCCijk(jz) / fabsf(HydroToCCijk(rho));
	}
      }

      /*      for_all_cells {
        if (fabsf(HydroToCCijk(rho)) > 0) {
	  Te(i,j,k) += (Pe_xx(i,j,k)-neUUe_xx(i,j,k)+Pe_yy(i,j,k)-neUUe_yy(i,j,k)+Pe_zz(i,j,k)-neUUe_zz(i,j,k))/(3.0*fabsf(HydroToCCijk(rho)));
	}
      }

      for_all_cells {
	if (fabsf(Te(i,j,k)) > 0) {
	  Fcollx(i,j,k) = 0.125*(field(i,j,k).jfx+field(i,j,k+1).jfx+field(i,j+1,k).jfx+field(i,j+1,k+1).jfx)*(field(i,j,k).cbx + field(i+1,j,k).cbx)
(field(i,j,k).cby + field(i,j+1,k).cby)*0.5
(field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5

}*/


      // d( n_e u_e)/dt
      if (global->dis_iter == 0) { //  this is the first iteration
	for_all_cells {
	  djedtx(i,j,k) += -HydroToCCijk(px);
	  djedty(i,j,k) += -HydroToCCijk(py);
	  djedtz(i,j,k) += -HydroToCCijk(pz);
	  
	  dBdtx(i,j,k) = -(field(i,j,k).cbx + field(i+1,j,k).cbx)*0.5;
	  dBdty(i,j,k) = -(field(i,j,k).cby + field(i,j+1,k).cby)*0.5;
	  dBdtz(i,j,k) = -(field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5;
	  
	}
	sim_log("first part of du_e/dt added");
	sim_log("first part of dB/dt added");
      }
      else if(global->dis_iter == (global->dis_nav - 1) ){ //this will be the last iteration
	for_all_cells {
	  djedtx(i,j,k) += HydroToCCijk(px);
	  djedty(i,j,k) += HydroToCCijk(py);
	  djedtz(i,j,k) += HydroToCCijk(pz);
	  
	  dBdtx(i,j,k) += (field(i,j,k).cbx + field(i+1,j,k).cbx)*0.5;
	  dBdty(i,j,k) += (field(i,j,k).cby + field(i,j+1,k).cby)*0.5;
	  dBdtz(i,j,k) += (field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5;
	  
	}
	sim_log("second part of du_e/dt added");
	sim_log("second part of dB/dt added");
      }
      //else... it's a middle cell - won't matter until window functions are put in
      //Minutiae - this isn't quite right, none of these are centered differences
      //When I throw in the window function, I can do the center ones cell-centered
      //and the start and end times with forward and backward differencing
      //will need to change the normalization also
      
      // now ions
      
      sp = find_species_name("ion", species_list);
      
      //clear_hydro(hydro, grid);
      //accumulate_hydro_p(hydro, sp->p, sp->np, sp->q_m, interpolator, grid);
      //synchronize_hydro(hydro, grid);
      
      clear_hydro_array(hydro_array);
      accumulate_hydro_p(hydro_array, sp, interpolator_array);
      synchronize_hydro_array(hydro_array);
      
      for_all_cells {
	Jix(i,j,k) = Jix(i,j,k) + HydroToCCijk(jx);
      }
      
      for_all_cells {
	Jiy(i,j,k) = Jiy(i,j,k) + HydroToCCijk(jy);
      }
      
      for_all_cells {
	Jiz(i,j,k) = Jiz(i,j,k) + HydroToCCijk(jz);
      }
      
      for_all_cells {
	Pi_xx(i,j,k) = Pi_xx(i,j,k) + HydroToCCijk(txx);
	Pi_yy(i,j,k) = Pi_yy(i,j,k) + HydroToCCijk(tyy);
	Pi_zz(i,j,k) = Pi_zz(i,j,k) + HydroToCCijk(tzz);
	Pi_xy(i,j,k) = Pi_xy(i,j,k) + HydroToCCijk(txy);
	Pi_xz(i,j,k) = Pi_xz(i,j,k) + HydroToCCijk(tzx);
	Pi_yz(i,j,k) = Pi_yz(i,j,k) + HydroToCCijk(tyz);
      }
      
      for_all_cells {
	
	ni(i,j,k) = ni(i,j,k) + HydroToCCijk(rho);
	//mix_i(i,j,k) = mix_i(i,j,k) + hydro(i,j,k).rho;
      }
      
      // E & B
      for_all_cells {
	Bx(i,j,k) = Bx(i,j,k) + (field(i,j,k).cbx + field(i+1,j,k).cbx)*0.5;
	By(i,j,k) = By(i,j,k) + (field(i,j,k).cby + field(i,j+1,k).cby)*0.5;
	Bz(i,j,k) = Bz(i,j,k) + (field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5;
	
	Ex(i,j,k) = Ex(i,j,k) + (field(i,j,k).ex+field(i,j+1,k).ex+field(i,j,k+1).ex+field(i,j+1,k+1).ex)*0.25;
	Ey(i,j,k) = Ey(i,j,k) + (field(i,j,k).ey+field(i+1,j,k).ey+field(i,j,k+1).ey+field(i+1,j,k+1).ey)*0.25;
	Ez(i,j,k) = Ez(i,j,k) + (field(i,j,k).ez+field(i+1,j,k).ez+field(i,j+1,k).ez+field(i+1,j+1,k).ez)*0.25;
	
	//      Jfx(i,j,k) = Jfx(i,j,k) + field(i,j,k).jfx;
	//      Jfy(i,j,k) = Jfy(i,j,k) + field(i,j,k).jfy;
	//      Jfz(i,j,k) = Jfz(i,j,k) + field(i,j,k).jfz;
      }
      

      // n_i*E
      for_all_cells {
	niEx(i,j,k) += HydroToCCijk(rho)*(field(i,j,k).ex + field(i,j+1,k).ex + field(i,j,k+1).ex + field(i,j+1,k+1).ex)*0.25;
	niEy(i,j,k) += HydroToCCijk(rho)*(field(i,j,k).ey + field(i+1,j,k).ey + field(i,j,k+1).ey + field(i+1,j,k+1).ey)*0.25;
	niEz(i,j,k) += HydroToCCijk(rho)*(field(i,j,k).ez + field(i+1,j,k).ez + field(i,j+1,k).ez + field(i+1,j+1,k).ez)*0.25;
      }
      
      // j_i x B
      for_all_cells {
	
	//  Assume both B and J are on Yee Mesh
	
	jixbx(i,j,k) += (HydroToCCijk(jy))*(field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5;
	jixbx(i,j,k) -= (HydroToCCijk(jz))*(field(i,j,k).cby + field(i,j+1,k).cby)*0.5;
	
	jixby(i,j,k) += (HydroToCCijk(jz))*(field(i,j,k).cbx + field(i+1,j,k).cbx)*0.5;
	jixby(i,j,k) -= (HydroToCCijk(jx))*(field(i,j,k).cbz + field(i,j,k+1).cbz)*0.5;

	jixbz(i,j,k) += (HydroToCCijk(jx))*(field(i,j,k).cby + field(i,j+1,k).cby)*0.5;
	jixbz(i,j,k) -= (HydroToCCijk(jy))*(field(i,j,k).cbx + field(i+1,j,k).cbx)*0.5;

      }
      
    // ni U_i U_i
      for_all_cells{
	if (fabsf(HydroToCCijk(rho)) > 0) {
	  niUUi_xx(i,j,k) += HydroToCCijk(jx) * HydroToCCijk(jx) / fabsf(HydroToCCijk(rho));
	  niUUi_yy(i,j,k) += HydroToCCijk(jy) * HydroToCCijk(jy) / fabsf(HydroToCCijk(rho));
	  niUUi_zz(i,j,k) += HydroToCCijk(jz) * HydroToCCijk(jz) / fabsf(HydroToCCijk(rho));
	  niUUi_xy(i,j,k) += HydroToCCijk(jx) * HydroToCCijk(jy) / fabsf(HydroToCCijk(rho));
	  niUUi_xz(i,j,k) += HydroToCCijk(jx) * HydroToCCijk(jz) / fabsf(HydroToCCijk(rho));
	  niUUi_yz(i,j,k) += HydroToCCijk(jy) * HydroToCCijk(jz) / fabsf(HydroToCCijk(rho));
	}
      }
      
      // d( n_i u_i)/dt
      if (global->dis_iter == 0) { //  this is the first iteration
	for_all_cells {
	  djidtx(i,j,k) = -HydroToCCijk(px);
	  djidty(i,j,k) = -HydroToCCijk(py);
	  djidtz(i,j,k) = -HydroToCCijk(pz);
	}
	sim_log("first part of du_i/dt added");
      }
      else if(global->dis_iter == (global->dis_nav - 1) ){ //this will be the last iteration
	for_all_cells {
	  djidtx(i,j,k) += HydroToCCijk(px);
	  djidty(i,j,k) += HydroToCCijk(py);
	  djidtz(i,j,k) += HydroToCCijk(pz);
	}
	sim_log("second part of du_i/dt added");
      }   
      
      
      global->dis_iter++;  // iteration count
      
    } // CLOSES sample==0 loop -> Allows restart dumping when sample!=0.

    // write files if it's the last iteration OR we need to save restart
    if (  (global->dis_iter == global->dis_nav) || global->write_restart || global->write_end_restart ) {
	//if (  (global->dis_iter == global->dis_nav)) {
	sim_log( "Ohms averaging:stop accumulation---dis_iter="<<global->dis_iter );
	sim_log( "Ohms averaging:stop accumulation---r="<<r );

      
	if (global->dis_iter == global->dis_nav) { // check if this is a regular write
	  for_all_cells {
	    Bx(i,j,k) = Bx(i,j,k)/global->dis_iter;
	    By(i,j,k) = By(i,j,k)/global->dis_iter;
	    Bz(i,j,k) = Bz(i,j,k)/global->dis_iter;

	    Ex(i,j,k) = Ex(i,j,k)/global->dis_iter;
	    Ey(i,j,k) = Ey(i,j,k)/global->dis_iter;
	    Ez(i,j,k) = Ez(i,j,k)/global->dis_iter;
	    
	    //	  Jfx(i,j,k) = Jfx(i,j,k)/global->dis_iter;
	    //	  Jfy(i,j,k) = Jfy(i,j,k)/global->dis_iter;
	    //	  Jfz(i,j,k) = Jfz(i,j,k)/global->dis_iter;
	    
	    Jex(i,j,k) = Jex(i,j,k)/global->dis_iter;
	    Jey(i,j,k) = Jey(i,j,k)/global->dis_iter;
	    Jez(i,j,k) = Jez(i,j,k)/global->dis_iter;

	    Jix(i,j,k) = Jix(i,j,k)/global->dis_iter;
	    Jiy(i,j,k) = Jiy(i,j,k)/global->dis_iter;
	    Jiz(i,j,k) = Jiz(i,j,k)/global->dis_iter;
	    
	    Pe_xx(i,j,k) = Pe_xx(i,j,k)/global->dis_iter;
	    Pe_yy(i,j,k) = Pe_yy(i,j,k)/global->dis_iter;
	    Pe_zz(i,j,k) = Pe_zz(i,j,k)/global->dis_iter;
	    Pe_xy(i,j,k) = Pe_xy(i,j,k)/global->dis_iter;
	    Pe_xz(i,j,k) = Pe_xz(i,j,k)/global->dis_iter;
	    Pe_yz(i,j,k) = Pe_yz(i,j,k)/global->dis_iter;
	    
	    Pi_xx(i,j,k) = Pi_xx(i,j,k)/global->dis_iter;
	    Pi_yy(i,j,k) = Pi_yy(i,j,k)/global->dis_iter;
	    Pi_zz(i,j,k) = Pi_zz(i,j,k)/global->dis_iter;
	    Pi_xy(i,j,k) = Pi_xy(i,j,k)/global->dis_iter;
	    Pi_xz(i,j,k) = Pi_xz(i,j,k)/global->dis_iter;
	    Pi_yz(i,j,k) = Pi_yz(i,j,k)/global->dis_iter;
	    
	    
	    //	  mix_e(i,j,k) = mix_e(i,j,k)/ne(i,j,k);	  
	    //mix_i(i,j,k) = mix_i(i,j,k)/ni(i,j,k);	  
	    
	    ne(i,j,k) = ne(i,j,k)/global->dis_iter;
	    ni(i,j,k) = ni(i,j,k)/global->dis_iter;
	    
	    //Ohm's law terms
	    neEx(i,j,k) /= global->dis_iter; //n_e E
	    neEy(i,j,k) /= global->dis_iter;
	    neEz(i,j,k) /= global->dis_iter;	  
	    
	    edotj(i,j,k) /= global->dis_iter;	  
	    
	    jexbx(i,j,k) /= global->dis_iter; //j_e x B
	    jexby(i,j,k) /= global->dis_iter;
	    jexbz(i,j,k) /= global->dis_iter;
	    
	    neUUe_xx(i,j,k) /= global->dis_iter; // n_e U_e U_e
	    neUUe_yy(i,j,k) /= global->dis_iter;	  
	    neUUe_zz(i,j,k) /= global->dis_iter;	  
	    neUUe_xy(i,j,k) /= global->dis_iter;	  
	    neUUe_xz(i,j,k) /= global->dis_iter;	  
	    neUUe_yz(i,j,k) /= global->dis_iter;	  
	    
	    // d(n_e u_e)/dt
	    /*Note that since the differences are centered between the timesteps
	      there are 1 fewer of them */
	    djedtx(i,j,k) =  djedtx(i,j,k)/(( global->dis_iter-1 )* dis_sample * grid->dt ); // ADAM: * dis_sample in denominator.
	    djedty(i,j,k) =  djedty(i,j,k)/(( global->dis_iter-1 )* dis_sample * grid->dt );
	    djedtz(i,j,k) =  djedtz(i,j,k)/(( global->dis_iter-1 )* dis_sample * grid->dt );
	    
	    // dB/dt - similar to d(n_e u_e)/dt
	    dBdtx(i,j,k) =  dBdtx(i,j,k)/(( global->dis_iter-1 )* dis_sample * grid->dt );  // ADAM: * dis_sample in denominator.
	    dBdty(i,j,k) =  dBdty(i,j,k)/(( global->dis_iter-1 )* dis_sample * grid->dt );
	    dBdtz(i,j,k) =  dBdtz(i,j,k)/(( global->dis_iter-1 )* dis_sample * grid->dt );
	    
	    //Ion Terms
	    niEx(i,j,k) /= global->dis_iter; //n_e E
	    niEy(i,j,k) /= global->dis_iter;
	    niEz(i,j,k) /= global->dis_iter;	    
	    
	    jixbx(i,j,k) /= global->dis_iter; //j_e x B
	    jixby(i,j,k) /= global->dis_iter;
	    jixbz(i,j,k) /= global->dis_iter;
	    
	    niUUi_xx(i,j,k) /= global->dis_iter; // n_e U_e U_e
	    niUUi_yy(i,j,k) /= global->dis_iter;	  
	    niUUi_zz(i,j,k) /= global->dis_iter;	  
	    niUUi_xy(i,j,k) /= global->dis_iter;	  
	    niUUi_xz(i,j,k) /= global->dis_iter;	  
	    niUUi_yz(i,j,k) /= global->dis_iter;
	    
	    // d(n_i u_i)/dt
	    /*Note that since the differences are centered between the timesteps
	      there are 1 fewer of them */
	    djidtx(i,j,k) =  djidtx(i,j,k)/(( global->dis_iter-1 )* dis_sample * grid->dt );  // ADAM: * dis_sample in denominator.
	    djidty(i,j,k) =  djidty(i,j,k)/(( global->dis_iter-1 )* dis_sample * grid->dt );
	    djidtz(i,j,k) =  djidtz(i,j,k)/(( global->dis_iter-1 )* dis_sample * grid->dt );
	    
	    
	  }; // normalize  edotj
	  global->dis_iter = 0;                                       // reset iteration count
	  sim_log("Time-averaging diagnostic: writing regular output");
	}  else {
	  sim_log("Time-averaging diagnostic: writing restart files");
	}	// if this is regular save

      if (rank() == 0) MESSAGE(("Time Averaging ---> start = %d   end=%d  steps=%d\n",global->dis_begin_int,step(),1+step()-global->dis_begin_int)); 

      // write the file
      sprintf(fname,"hydro-ave/%d/T.%d",NUMFOLD,global->dis_begin_int);
      dump_mkdir(fname);
     
      sprintf(fname,"hydro-ave/%d/T.%d/%s.%d.%d",NUMFOLD,global->dis_begin_int,"data-ave",global->dis_begin_int,(int) rank()); \
      if ( ! (fp.open(fname, io_write)==ok) ) ERROR(("Cannot open file."));
      sim_log( "saving the in-place averaging diagnostics Bx" );
      fp.write(Bx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics By" );
      fp.write(By, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Bz" );
      fp.write(Bz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Ex" );
      fp.write(Ex, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Ey" );
      fp.write(Ey, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Ez" );
      fp.write(Ez, nx*ny*nz );
      //      sim_log( "saving the in-place averaging diagnostics Jfx" );
      //      fp.write(Jfx, nx*ny*nz );
      //      sim_log( "saving the in-place averaging diagnostics Jfy" );
      //      fp.write(Jfy, nx*ny*nz );
      //      sim_log( "saving the in-place averaging diagnostics Jfz" );
      //      fp.write(Jfz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Jex" );
      fp.write(Jex, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Jey" );
      fp.write(Jey, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Jez" );
      fp.write(Jez, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Jix" );
      fp.write(Jix, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Jiy" );
      fp.write(Jiy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics Jiz" );
      fp.write(Jiz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pe-xx" );
      fp.write(Pe_xx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pe-yy" );
      fp.write(Pe_yy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pe-zz" );
      fp.write(Pe_zz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pe-xy" );
      fp.write(Pe_xy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pe-xz" );
      fp.write(Pe_xz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pe-yz" );
      fp.write(Pe_yz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pi-xx" );
      fp.write(Pi_xx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pi-yy" );
      fp.write(Pi_yy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pi-zz" );
      fp.write(Pi_zz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pi-xy" );
      fp.write(Pi_xy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pi-xz" );
      fp.write(Pi_xz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :Pi-yz" );
      fp.write(Pi_yz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :ne" );
      fp.write(ne, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics :ni" );
      fp.write(ni, nx*ny*nz );
      //      sim_log( "saving the in-place averaging diagnostics :mix_e" );
      //fp.write(mix_e, nx*ny*nz );
      //sim_log( "saving the in-place averaging diagnostics :mix_i" );
      //fp.write(mix_i, nx*ny*nz );
      //Ohm's law terms
      sim_log( "saving the in-place averaging diagnostics ne*Ex" );
      fp.write(neEx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics ne*Ey" );
      fp.write(neEy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics ne*Ez" );
      fp.write(neEz, nx*ny*nz );	    

      sim_log( "saving the in-place averaging diagnostics E*J" );
      fp.write(edotj, nx*ny*nz );	    

      sim_log( "saving the in-place averaging diagnostics (j_e x B)_x" );
      fp.write(jexbx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (j_e x B)_y" );
      fp.write(jexby, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (j_e x B)_z" );
      fp.write(jexbz, nx*ny*nz );
 
      sim_log( "saving the in-place averaging diagnostics (n_e U_e U_e)_xx" );
      fp.write(neUUe_xx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_e U_e U_e)_yy" );
      fp.write(neUUe_yy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_e U_e U_e)_zz" );
      fp.write(neUUe_zz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_e U_e U_e)_xy" );
      fp.write(neUUe_xy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_e U_e U_e)_xz" );
      fp.write(neUUe_xz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_e U_e U_e)_yz" );
      fp.write(neUUe_yz, nx*ny*nz );

      sim_log( "saving the in-place averaging diagnostics d(n_e u_e)/dt_x" );
      fp.write(djedtx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics d(n_e u_e)/dt_y" );
      fp.write(djedty, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics d(n_e u_e)/dt_z" );
      fp.write(djedtz, nx*ny*nz );

      sim_log( "saving the in-place averaging diagnostics dB/dt_x" );
      fp.write(dBdtx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics dB/dt_y" );
      fp.write(dBdty, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics dB/dt_z" );
      fp.write(dBdtz, nx*ny*nz );

      //Ion Terms
      sim_log( "saving the in-place averaging diagnostics ni*Ex" );
      fp.write(niEx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics ni*Ey" );
      fp.write(niEy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics ni*Ez" );
      fp.write(niEz, nx*ny*nz );	    	    

      sim_log( "saving the in-place averaging diagnostics (j_i x B)_x" );
      fp.write(jixbx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (j_i x B)_y" );
      fp.write(jixby, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (j_i x B)_z" );
      fp.write(jixbz, nx*ny*nz );
 
      sim_log( "saving the in-place averaging diagnostics (n_i U_i U_i)_xx" );
      fp.write(niUUi_xx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_i U_i U_i)_yy" );
      fp.write(niUUi_yy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_i U_i U_i)_zz" );
      fp.write(niUUi_zz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_i U_i U_i)_xy" );
      fp.write(niUUi_xy, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_i U_i U_i)_xz" );
      fp.write(niUUi_xz, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics (n_i U_i U_i)_yz" );
      fp.write(niUUi_yz, nx*ny*nz );

      sim_log( "saving the in-place averaging diagnostics d(n_i u_i)/dt_x" );
      fp.write(djidtx, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics d(n_i u_)/dt_y" );
      fp.write(djidty, nx*ny*nz );
      sim_log( "saving the in-place averaging diagnostics d(n_i u_i)/dt_z" );
      fp.write(djidtz, nx*ny*nz );

      fp.close();
    } // end of file write

  } // end of inside the everaging interval 

} // end of time-averaging diagnostic


