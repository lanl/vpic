/* 
   Utility to post-process field and hydro dumps from VPIC into a form
   easily read by various postprocessing facilities such as IDL.  

   The utility uses an input deck to specify problem topology and then 
   creates arrays from the data files that have the specified data and
   topology.  The user is permitted to specify strides for the various 
   array dimensions in order to create more manageable file sizes.  
   (This is particularly important on very large problems).  Strides do
   not need to divide evenly into the simulation array dimensions. 

   Usage: data_join [fname base] [input deck] <optional tag>
          data_join template

   The output data consist of a set of integer array limits followed by 
   an array of binary single precision data of specified dimensions 
   containing the data fields requested.

   The optional tag is an identifier string appended to each output file. 
   For example, to generate a data files "ex.0.bin", one would provide a
   tag "0".  

   ------------------------------------------------------------------------
   Last modified: Brian Albright, X-1, 1/09/2005. 
*/ 


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

typedef unsigned short int material_id;

/* Various structs */ 

struct v0_grid_data {
  int step, nx, ny, nz;
  float dt, dx, dy, dz, x0, y0, z0, cvac, eps0, damp;
  int mp_rank, mp_num_proc;
  int species_id;
  float q_m;
};

struct array_header {
  int size; 
  int ndim;
  int dim;
};

typedef struct field {
  float ex, ey, ez;                /* Electric field */
  float cbx, cby, cbz;             /* Magnetic field */
  float tcax, tcay, tcaz;          /* TCA fields */
  float jfx, jfy, jfz;             /* Free current */
  material_id ematx, ematy, ematz; /* Material at edge centers */
  material_id fmatx, fmaty, fmatz; /* Material at face centers */
  material_id nmat, cmat;          /* Material at cell centers and nodes */
  float rhof, rhob;                /* Free and bound charge density */
  float div_e_err, div_b_err;      /* Divergence errors */
} field_t;

typedef struct hydro {
  float rho;           /* Charge density         => < q f > */
  float jx, jy, jz;    /* Current density        => < q v_i f > */
  float ke;            /* Kinetic energy density => < m c^2 (gamma-1) f > */
  float px, py, pz;    /* Momentum density       => < p_i f > */
  float txx, tyy, tzz; /* Stress diagonal        => < p_i v_j f >, i==j */
  float tyz, tzx, txy; /* Stress off-diagonal    => < p_i v_j f >, i!=j */
  float pad0, pad1;    /* 16-byte align the structure */
} hydro_t;

void ierror( char *msg ); 
void *emalloc( size_t size ); 
void store_frame_data( FILE *fp, int frames );
void print_template( void );

#define READ(type,value,file) do {                \
    type __READ_tmp;                              \
    fread( &__READ_tmp, sizeof(type), 1, file );  \
    (value) = __READ_tmp;                         \
  } while(0); 

#define ABORT(cond) {if (cond) exit(0);}

#define ERROR(m) {     \
  char errmsg[512];    \
  sprintf(errmsg, m);  \
  ierror(errmsg);      \
}

/* Macros cribbed from VPIC partition.c */ 

#define RANK_TO_INDEX(rank,ix,iy,iz) {                  \
  int _ix, _iy, _iz;                                    \
  _ix  = (rank);  /* ix = ix + gpx*( iy + gpy*iz ) */   \
  _iy  = _ix/gpx; /* iy = iy + gpy*iz */                \
  _ix -= _iy*gpx; /* ix = ix */                         \
  _iz  = _iy/gpy; /* iz = iz */                         \
  _iy -= _iz*gpy; /* iy = iy */                         \
  (ix) = _ix;                                           \
  (iy) = _iy;                                           \
  (iz) = _iz;                                           \
} 

#define INDEX_TO_RANK(ix,iy,iz,rank) {                  \
  int _ix, _iy, _iz;                                    \
   /* Wrap processor index periodically */              \
   _ix = (ix) % gpx; while(_ix<0) _ix += gpx;           \
   _iy = (iy) % gpy; while(_iy<0) _iy += gpy;           \
   _iz = (iz) % gpz; while(_iz<0) _iz += gpz;           \
   /* Compute the rank */                               \
  (rank) = _ix + gpx*( _iy + gpy*_iz );                 \
} 

#define flush_line( fp ) { char c; do { c=fgetc( fp ); } while ( c!='\n' && c!=EOF ); }
#define INFILE_TYPE_NONE  0
#define INFILE_TYPE_FIELD 1
#define INFILE_TYPE_HYDRO 2

#define foreach_file    for ( ifile=0; ifile<nproc; ++ifile )

/* For Yee mesh offsets */ 
#define EOFF 1
#define COFF 0

/*************************************************************************
  Main routine                                   
**************************************************************************/

int main( int argc, char *argv[] ) {
  int nproc, ifile, nframes;
  FILE **fparr, *ofp[38], *ologp, *input_deck; 
  int *nxarr; 
  int nx, nx_output, nx_total, nx_ave; 
  int nvx, nvx_output=200; 
  float dv; 
  float *data_vec, *odata_vec; 
  union idatum  {
    field_t fd;
    hydro_t hd;
  };
  union idatum idat; 
  char fname[256], ofname[256], errmsg[512], ideck_name[256], dummy[256];
  int not_eof, infile_type=INFILE_TYPE_NONE; 
  int i, xstride=1, ystride=1, zstride=1, gpx=0, gpy=0, gpz=0;
  int xsize, ysize, zsize, osize[38], ilow=0, ihigh=0;
  int iflag[38] = { 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0,
                    0, 0, 0 };
  float *odata[38] = { NULL, NULL, NULL, NULL, NULL, 
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL };
  char *fname_list[] = { "ex", "ey", "ez", "cbx", "cby", "cbz", 
                         "tcax", "tcay", "tcaz", "jfx", "jfy", "jfz",
                         "ematx", "ematy", "ematz", "fmatx", "fmaty", "fmatz", 
                         "nmat", "cmat",
                         "rhof", "rhob", "dive_err", "divb_err", 
                         "rho", "jx", "jy", "jz", "ke", "px", "py", "pz", 
                         "txx", "tyy", "tzz", "tyz", "tzx", "txy" }; 
  int offsets[38][3] = { { COFF, EOFF, EOFF }, /* Ex */ 
                         { EOFF, COFF, EOFF }, /* Ey */ 
                         { EOFF, EOFF, COFF }, /* Ez */ 
                         { EOFF, COFF, COFF }, /* cBx */ 
                         { COFF, EOFF, COFF }, /* cBy */ 
                         { COFF, COFF, EOFF }, /* cBz */ 
                         { COFF, EOFF, EOFF }, /* TCAx */ 
                         { EOFF, COFF, EOFF }, /* TCAy */ 
                         { EOFF, EOFF, COFF }, /* TCAz */ 
                         { COFF, EOFF, EOFF }, /* jfx */ 
                         { EOFF, COFF, EOFF }, /* jfy */ 
                         { EOFF, EOFF, COFF }, /* jfz */ 
                         { COFF, EOFF, EOFF }, /* ematx */ 
                         { EOFF, COFF, EOFF }, /* ematy */ 
                         { EOFF, EOFF, COFF }, /* ematz */ 
                         { EOFF, COFF, COFF }, /* fmatx */ 
                         { COFF, EOFF, COFF }, /* fmaty */ 
                         { COFF, COFF, EOFF }, /* fmatz */ 
                         { EOFF, EOFF, EOFF }, /* nmatz */ 
                         { COFF, COFF, COFF }, /* cmatz */ 
                         { EOFF, EOFF, EOFF }, /* rhof */ 
                         { EOFF, EOFF, EOFF }, /* rhob */ 
                         { EOFF, EOFF, EOFF }, /* div e err */ 
                         { COFF, COFF, COFF }, /* div b err */ 
                         { COFF, COFF, COFF }, /* rho */ /* I think hydro accum to cell centers */
                                                         /* (if nodes, then this is off by half a cell) */ 
                         { COFF, COFF, COFF }, /* jx */
                         { COFF, COFF, COFF }, /* jy */
                         { COFF, COFF, COFF }, /* jz */
                         { COFF, COFF, COFF }, /* ke */
                         { COFF, COFF, COFF }, /* px */
                         { COFF, COFF, COFF }, /* py */
                         { COFF, COFF, COFF }, /* pz */
                         { COFF, COFF, COFF }, /* txx */
                         { COFF, COFF, COFF }, /* tyy */
                         { COFF, COFF, COFF }, /* tzz */
                         { COFF, COFF, COFF }, /* tyz */
                         { COFF, COFF, COFF }, /* tzx */
                         { COFF, COFF, COFF }  /* txy */  };

  if ( argc!=3 && argc!=4 ) {
    if  ( argc==2 && !strcmp(argv[1], "template") ) print_template(); 
    fprintf( stderr, 
             "Usage: %s [fname base] [input deck] <tag>\n       %s template\n", 
	     argv[0], argv[0] ); 
    return 0;
  } 

  /* Open the input deck file */ 
  sprintf( ideck_name, "%s", argv[2] ); 
  if ( !(input_deck = fopen( ideck_name, "r")) ) 
    ERROR(("Cannot open file: %s\n", ideck_name));         

  /* Parse input deck */ 
  printf("Parsing input deck....\n"); 
  fscanf( input_deck, "%d %d %d", &gpx, &gpy, &gpz ); flush_line( input_deck ); 
  ABORT( gpx<=0 || gpy<=0 || gpz<=0 ); 
  nproc=gpx*gpy*gpz;
  fscanf( input_deck, "%d %d %d", &xstride, &ystride, &zstride ); flush_line( input_deck ); 
  ABORT( xstride<=0 || ystride<=0 || zstride<=0 ); 
  for ( i=0; i<38; ++i ) { 
    fscanf( input_deck, "%d", &iflag[i] ); flush_line( input_deck );  
  }
  printf("Done.\n"); 

  if ( !(ologp=fopen("data_join.log", "w")) ) ERROR(("Cannot open file: data_join.log\n")); 
  for ( i=0; i<38; ++i ) 
    if ( iflag[i] ) printf("Data set %s is requested.\n", fname_list[i] );  

  /* Initialize output file handles */ 
  for ( i=0; i<38; ++i ) {
    if ( argc==3 ) sprintf( ofname, "%s.bin",    fname_list[i] ); 
    else           sprintf( ofname, "%s.%s.bin", fname_list[i], argv[3] ); 
    if ( iflag[i] && !(ofp[i]=fopen(ofname, "wb")) ) 
      ERROR(("Cannot open file: %s\n", ofname));         
  }

  /* Allocate space for input file handles. */ 
  fparr=emalloc((size_t)nproc*sizeof(*fparr));
  nxarr=emalloc((size_t)nproc*sizeof(*nxarr)); 
  foreach_file fparr[ifile]=emalloc(sizeof(**fparr)); 

  /* Check #1 for file consistency */ 
  for ( i=0; i<24; ++i )  if ( iflag[i] ) infile_type=INFILE_TYPE_FIELD; 
  for ( i=24; i<38; ++i ) 
    if ( iflag[i] ) {
      if ( infile_type!=INFILE_TYPE_NONE ) {
        ERROR(("Multiple infile types. Need to process field and hydro data separately!")); 
      } else { 
        infile_type=INFILE_TYPE_HYDRO; 
        break;
      }
    }
  if ( infile_type==INFILE_TYPE_NONE ) {
    fprintf( stderr, "Nothing to process.  Terminating....\n" );
    exit(0);
  }

  /************************************************************************* 
     Read data 
  **************************************************************************/
 
  /* Loop over all of the input files. */ 
  foreach_file {
    struct v0_grid_data gd;
    struct array_header ah; 
    int ip, px, py, pz, ix, iy, iz;
    int nparticles_this_file=0, nparticles_output_from_this_file=0; 

    sprintf(fname, "%s.%d", argv[1], ifile); 
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) ERROR( "Cannot open input data files.\n" ); 
    printf( "Processing data from file: %s\n", fname ); 

    /* Read the WRITE_HEADER_V0 data */ 
    {
      char c;
      short s;
      int n;
      float f;
      double d; 

      /* Check consistency of machine representation of data */ 
      READ(char,     c,fparr[ifile]);  ABORT(c!=8                );
      READ(char,     c,fparr[ifile]);  ABORT(c!=2                );
      READ(char,     c,fparr[ifile]);  ABORT(c!=4                );
      READ(char,     c,fparr[ifile]);  ABORT(c!=4                );
      READ(char,     c,fparr[ifile]);  ABORT(c!=8                );
      READ(short int,s,fparr[ifile]);  ABORT(s!=(short int)0xcafe);
      READ(int,      n,fparr[ifile]);  ABORT(n!=(int)0xdeadbeef  );
      READ(float,    f,fparr[ifile]);  ABORT(f!=1.0              );
      READ(double,   d,fparr[ifile]);  ABORT(d!=1.0              );

      /* Version and dump type */
      READ(int,n,fparr[ifile]); /* Version */ 
      READ(int,n,fparr[ifile]); 
      /* Check #2 for file consistency */ 
      if (    n==1 /* field */ && infile_type==INFILE_TYPE_HYDRO 
           || n==2 /* hydro */ && infile_type==INFILE_TYPE_FIELD )
        ERROR(("Incompatable data file type for the requested fields."));  

      /* High level information */
      READ(int,  gd.step,      fparr[ifile]);
      READ(int,  gd.nx,  fparr[ifile]);
      READ(int,  gd.ny,  fparr[ifile]);
      READ(int,  gd.nz,  fparr[ifile]);
      READ(float,gd.dt,  fparr[ifile]);
      READ(float,gd.dx,  fparr[ifile]);
      READ(float,gd.dy,  fparr[ifile]);
      READ(float,gd.dz,  fparr[ifile]);
      READ(float,gd.x0,  fparr[ifile]);
      READ(float,gd.y0,  fparr[ifile]);
      READ(float,gd.z0,  fparr[ifile]);
      READ(float,gd.cvac,fparr[ifile]);
      READ(float,gd.eps0,fparr[ifile]);
      READ(float,gd.damp,fparr[ifile]);
      READ(int,  n,      fparr[ifile]);
      READ(int,  n,      fparr[ifile]);
      /* Species parameters */
      READ(int,  n,fparr[ifile]);
      READ(float,f,fparr[ifile]);
      
      /* Array size/dimension information */ 
      READ(int,ah.size,  fparr[ifile]); 
      READ(int,ah.ndim,  fparr[ifile]);
      READ(int,ah.dim,   fparr[ifile]);  

      /* Check #3 for file consistency; probably overkill */ 
      if ( ah.size!=(infile_type==INFILE_TYPE_HYDRO ? sizeof(hydro_t) : sizeof(field_t)) )
        ERROR(("Incompatible data struct sizes.")); 
    }

    RANK_TO_INDEX( ifile, px, py, pz ); 

    /* Allocate space for output arrays; write limits to log file. */ 
    if ( !ifile ) 
      for ( i=0; i<38; ++i )
        if ( iflag[i] ) {
          xsize=(gpx*gd.nx+offsets[i][0])/xstride;
          ysize=(gpy*gd.ny+offsets[i][1])/ystride;
          zsize=(gpz*gd.nz+offsets[i][2])/zstride;
          osize[i]=xsize*ysize*zsize; 
          odata[i]=emalloc( osize[i] ); 

          /* DEBUG */ 
          printf("px, py, pz = %d %d %d\n", gpx, gpy, gpz );
          printf("osize=%d\n", osize[i]);
          printf("gd.nx, gd.ny, gd.nz = %d %d %d\n", gd.nx, gd.ny, gd.nz );

          fprintf(ologp, "%s %d %d %d\n", fname_list[i], xsize, ysize, zsize ); 
        }

    /* DEBUG */ 
    exit(0);

    for ( ix=0; ix<=gd.nx; ++ix ) {
      for ( iy=0; iy<=gd.ny; ++iy ) {
        for ( iz=0; iz<=gd.nz; ++iz ) {
          int sz, iox, ioy, ioz;
          ilow =(infile_type==INFILE_TYPE_FIELD ? 0 : 24);
          ihigh=(infile_type==INFILE_TYPE_FIELD ? 24 : 38);
          sz   =(infile_type==INFILE_TYPE_FIELD ? sizeof(field_t) : sizeof(hydro_t)); 
          fread( &idat, sz, 1, fparr[ifile] );
#         define COMPUTE_INDEX_1D(PX,NX,IX) ((PX)*(NX)+(IX))
          if ( !(COMPUTE_INDEX_1D(px,gd.nx,ix)%xstride) && 
               !(COMPUTE_INDEX_1D(py,gd.ny,iy)%ystride) &&
               !(COMPUTE_INDEX_1D(pz,gd.nz,iz)%zstride) && 
               COMPUTE_INDEX_1D(px,gd.nx,ix)<px*gd.nx+offsets[i][0] &&
               COMPUTE_INDEX_1D(py,gd.ny,iy)<py*gd.ny+offsets[i][1] &&
               COMPUTE_INDEX_1D(pz,gd.nz,iz)<pz*gd.nz+offsets[i][2] )  {
            for ( i=ilow; i<ihigh; ++i ) {
              if ( iflag[i] ) {
                int nox, noy, noz;
                iox=COMPUTE_INDEX_1D(px,gd.nx,ix)/xstride; 
                ioy=COMPUTE_INDEX_1D(py,gd.ny,iy)/ystride;
                ioz=COMPUTE_INDEX_1D(pz,gd.nz,iz)/zstride;
                nox=(px*gd.nx+offsets[i][0])/xstride;   
                noy=(py*gd.ny+offsets[i][1])/ystride;   
                noz=(pz*gd.nz+offsets[i][2])/zstride;   
#               define INDEX(IX,IY,IZ) (IX*noy+IY)*noz+IZ
                if ( ix<gd.nx+offsets[i][0] && 
                     iy<gd.ny+offsets[i][1] && 
                     iz<gd.nz+offsets[i][2] ) {
                  switch ( i ) {
                  case 0:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.ex;
                  case 1:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.ey;
                  case 2:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.ez;
                  case 3:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.cbx;
                  case 4:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.cby;
                  case 5:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.cbz;
                  case 6:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.tcax;
                  case 7:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.tcay;
                  case 8:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.tcaz;
                  case 9:  odata[i][INDEX(iox,ioy,ioz)]=idat.fd.jfx;
                  case 10: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.jfy;
                  case 11: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.jfz;
                  case 12: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.ematx;
                  case 13: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.ematy;
                  case 14: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.ematz;
                  case 15: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.fmatx;
                  case 16: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.fmaty;
                  case 17: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.fmatz;
                  case 18: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.nmat;
                  case 19: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.cmat;
                  case 20: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.rhof;
                  case 21: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.rhob;
                  case 22: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.div_e_err;
                  case 23: odata[i][INDEX(iox,ioy,ioz)]=idat.fd.div_b_err;
                  case 24: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.rho;
                  case 25: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.jx;
                  case 26: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.jy;
                  case 27: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.jz;
                  case 28: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.ke;
                  case 29: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.px;
                  case 30: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.py;
                  case 31: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.pz;
                  case 32: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.txx;
                  case 33: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.tyy;
                  case 34: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.tzz;
                  case 35: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.tyz;
                  case 36: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.tzx;
                  case 37: odata[i][INDEX(iox,ioy,ioz)]=idat.hd.txy;
                  }     
                }
              }
            }
          }

        }
      }
    }

  }

  /* Data are stored.  Now write to file */ 
  for ( i=ilow; i<ihigh; ++i ) 
    if ( iflag[i] ) fwrite( odata[i], sizeof(float), osize[i], ofp[i] );  

  /* Close input and output files */ 
  foreach_file fclose(fparr[ifile]); 
  for ( i=ilow; i<ihigh; ++i ) if ( iflag[i] ) fclose(ofp[i]);
  fclose( ologp ); 

  return 0;
}


/************************************************************************ 
   Print an error message to stderr and exit program. 
************************************************************************/ 

void ierror( char *msg ) {
  fprintf( stderr, msg ); 
  exit(1); 
}


/*********************************************************************** 
   Allocate space with an error check that aborts if malloc is
   unsuccessful.  Returns a pointer to allocated space. 
***********************************************************************/ 

void *emalloc( size_t size ) {
  void *ptr = NULL; 

  if ( size<0 )
    ierror( "emalloc: request for array of negative size.\n" ); 
  if ( !(ptr=malloc(size)) )
    ierror( "emalloc: cannot allocate space.\n" ); 
  return ptr; 
}


/*********************************************************************** 
   Print a template input deck for this utility to stdout.  
***********************************************************************/ 

void print_template( void ) {
  printf( "1 1 1 ! number of processors in x, y, z directions\n"  
          "1 1 1 ! stride to use in printing arrays in x, y, z\n"
          "0     ! create ex mesh      (E field)  ----------------------------v in field data \n"
          "0     ! create ey mesh\n"
          "0     ! create ez mesh\n"
          "0     ! create cbx mesh     (B field)\n"
          "0     ! create cby mesh\n"
          "0     ! create cbz mesh\n"
          "0     ! create tcax mesh    (TCA fields)\n"
          "0     ! create tcay mesh\n"
          "0     ! create tcaz mesh\n"
          "0     ! create jfx mesh     (Free current)\n"
          "0     ! create jfy mesh\n"
          "0     ! create jfz mesh\n"
          "0     ! create ematx mesh   (Material at edge center) \n"
          "0     ! create ematy mesh\n"
          "0     ! create ematz mesh\n"
          "0     ! create fmatx mesh   (Material at face center) \n"
          "0     ! create fmaty mesh\n"
          "0     ! create fmatz mesh\n"
          "0     ! create nmat mesh    (Material at node)\n"
          "0     ! create cmat mesh    (Material at cell center)\n"
          "0     ! create rhof mesh    (Free charge density)\n"
          "0     ! create rhob mesh    (Bound charge density)\n"
          "0     ! create diveerr mesh (Div E error)\n"
          "0     ! create divberr mesh (Div B error)---------------------------^ in field data \n"
          "0     ! create rho mesh     (Charge density         => < q f >)-----v in hydro data \n"
          "0     ! create jx mesh      (Current density        => < q v_i f >)\n"
          "0     ! create jy mesh\n"
          "0     ! create jz mesh\n"
          "0     ! create ke mesh      (Kinetic energy density => < m c^2 (gamma-1) f >) \n"
          "0     ! create px mesh      (Momentum density       => < p_i f >) \n"
          "0     ! create py mesh\n"
          "0     ! create pz mesh\n"
          "0     ! create txx mesh     (Stress diagonal        => < p_i v_j f >, i==j)\n"
          "0     ! create tyy mesh\n"
          "0     ! create tzz mesh\n"
          "0     ! create tyz mesh     (Stress off-diagonal    => < p_i v_j f >, i!=j) \n"
          "0     ! create tzx mesh\n"
          "0     ! create txy mesh     ----------------------------------------^ in hydro data \n" );
  exit(0);
}
