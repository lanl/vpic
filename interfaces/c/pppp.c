/* 
   Utility to process vpic particle dump prior to postprocessing in IDL. 
   The utility name "pppp" stands for "particle pre-postprocessor". 

   This utility will convert the VPIC particle data (which are of 
   necessity written with local position data) into positions with
   global data by operating on the particle index.  

   Usage:  pppp 
   [num processors] [input filename prefix] [interval]

   The input filename prefix works as follows:  Suppose one has 32 files
   with names eparticle.1000.0 ... eparticle.1000.32 ; to process these
   data, one specifies a prefix of eparticle.1000

   The interval argument is an optional particle interval.  Set it equal
   to ten, e.g., in order to write only every 10th particle in the data. 

   The output data consist of binary single precision data containing the 
   following seven fields:

   { x, y, z, ux, uy, uz, q }

   x, y, z:	Position of particle in global coordinates
   ux, uy, uz:	Normalized momenta of particles (p/c)
   q:		Charge of particle (proportional to its weight)

   Modify the FILTER macro and recompile in order to select only certain 
   particles, e.g., those within a restricted spatial range or those having 
   energies above a threshold.

   Sample IDL file ppp.pro to use these data:

   ;==========================================================================
   ; Particle postprocessor.  Works off of particle data file prepared using
   ; pppp.c particle pre-postprocessor.  
   ;
   pro ppp, pdata
     fname = dialog_pickfile(/read, filter = '*particle*.bin')
     openr, unit1, fname, /get_lun, /swap_if_big_endian
     particlestruct = {  x:0.0,  y:0.0,  z:0.0, $    ; particle spatial coords
                        ux:0.0, uy:0.0, uz:0.0, $    ; particle normalized momenta
                         q:0.0 }                     ; particle charge
     pdata=assoc(unit1,particlestruct)
   end
   ;=========================================================================

   IDL usage: 
     .run ppp
     ppp, pdata
     print, pdata(1)
   
   ------------------------------------------------------------------------
   Last modified: Brian Albright, X-1, 4/24/2005. 
*/ 


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>


/* To handle more than 2 billion particles, we need an int64 type */ 

#ifndef INT64_TYPE
#if LONG_MAX==2147483647L
/* Long int is 32-bits ... assume long long int is 64-bits */
/* Warning: This is not ANSI-C89 compliant */
__extension__
typedef long long int int64;
#define INT64_TYPE int64
#else
/* Long int is not 32-bits ... assume long int is a 64-bits */
#define INT64_TYPE long int
#endif
#endif

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

struct ipdata {
  float dx, dy, dz; 
  int i;
  float ux, uy, uz; 
  float q; 
};

struct opdata {
  float x, y, z;
  float ux, uy, uz;
  float q;
};

void ierror( char *msg ); 
void *emalloc( size_t size ); 
void store_frame_data( FILE *fp, int frames );

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

#define foreach_file    for ( ifile=0; ifile<nproc; ++ifile )

/* Speed up I/0 by reading blocks of up to PBUF_SIZE particles at a time */ 
#define PBUF_SIZE  32768 /* 1MB of particles */ 

/* Set filter to something else and recompile to filter particles */ 
#define FILTER 1

#if 0
/* Example:  Filter only particles with positive x momentum and position 
             from 1 to 1.5 in y : */ 
	     
#define FILTER ( po.ux>0 && po.x>=1 && po.x<=1.5 )

#endif 

/*************************************************************************
  Main routine                                   
**************************************************************************/

int main( int argc, char *argv[] ) {
  int nproc, ifile, nframes;
  FILE **fparr, *ofp, *ologp; 
  int *nxarr; 
  int nx, nx_output, nx_total, nx_ave; 
  int nvx, nvx_output=200; 
  float dv; 
  float *data_vec, *odata_vec; 
  char fname[256], errmsg[512], ofname[256];
  int freq=1;
  int not_eof; 
  INT64_TYPE nparticles_total=0, nparticles_outfile=0;

  if ( argc!=3 && argc!=4 ) {
    fprintf( stderr, "Usage:  %s "
             "[num processors] [input filename prefix] [interval]\n\n", 
             argv[0] ); 
    return 0;
  } 
  
  /* Obtain runtime data from comand line. */ 
  nproc=atoi(argv[1]); 
  if ( nproc<=0 ) ERROR( "Bad number of input files requested.\n" );
  if ( argc==4 ) {
    freq=atoi(argv[3]);
    if ( freq<=0 ) ERROR( "Bad particle dump interval specified.\n" ); 
  }

  /* Allocate space for input file handles. */ 
  fparr=emalloc((size_t)nproc*sizeof(*fparr));
  nxarr=emalloc((size_t)nproc*sizeof(*nxarr)); 
  foreach_file fparr[ifile]=emalloc(sizeof(**fparr)); 

  /* Open the output file */ 
  sprintf( ofname, "%s.bin", argv[2] ); 
  if ( !(ofp = fopen( ofname, "wb")) ) ERROR(("Cannot open file: %s\n", ofname));         

  /************************************************************************* 
     Write electron phase space movie data. 
  **************************************************************************/
 
  foreach_file {
    struct v0_grid_data   gd;
    struct array_header   ah; 
    int ip, nxp2, nyp2, nzp2;
    int nparticles_this_file=0, nparticles_output_from_this_file=0; 

    sprintf(fname, "%s.%d", argv[2], ifile); 
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) ERROR( "Cannot open particle data files.\n" ); 
    printf( "Processing particles in file: %s  ", fname ); 

    /* Read the WRITE_HEADER_V0 data */ 
    {
      char c;
      short s;
      int n;
      float f;
      double d; 

      READ(char,     c,fparr[ifile]);  ABORT(c!=8                );
      READ(char,     c,fparr[ifile]);  ABORT(c!=2                );
      READ(char,     c,fparr[ifile]);  ABORT(c!=4                );
      READ(char,     c,fparr[ifile]);  ABORT(c!=4                );
      READ(char,     c,fparr[ifile]);  ABORT(c!=8                );
      READ(short int,s,fparr[ifile]);  ABORT(s!=(short int)0xcafe);
      READ(int,      n,fparr[ifile]);  ABORT(n!=(int)0xdeadbeef  );
      READ(float,    f,fparr[ifile]);  ABORT(f!=1.0              );
      READ(double,   d,fparr[ifile]);  ABORT(d!=1.0              );

      /* Dump type and header version */
      READ(int,n,fparr[ifile]); 
      READ(int,n,fparr[ifile]); 

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
      
      /* Particle array information */ 
      READ(int,ah.size,  fparr[ifile]); 
      READ(int,ah.ndim,  fparr[ifile]);
      READ(int,ah.dim,   fparr[ifile]);  
    }

    /* Now loop over particles and write data to ofp stream */ 

    /* Particles are read in blocks of PBUF_SIZE and are written 
       to a buffer before output to speed up I/O */ 

    nxp2=gd.nx+2;  
    nyp2=gd.ny+2;  
    nzp2=gd.nz+2; 
    for ( ip=0; ip<ah.dim; ip+=PBUF_SIZE ) {
      struct ipdata pi[PBUF_SIZE];
      struct opdata po, po_buf[PBUF_SIZE];
      int n_buf, ibuf, obuf; 

      n_buf=( ip+PBUF_SIZE>ah.dim ? ah.dim-ip : PBUF_SIZE ); 
      fread( &pi, sizeof(struct ipdata), n_buf, fparr[ifile] ); 
      nparticles_total     += n_buf;
      nparticles_this_file += n_buf; 
      obuf=0; 
      for ( ibuf=0; ibuf<n_buf; ++ibuf ) {
        if ( (ip+ibuf)%freq==0 ) {
          int ix, iy, iz, i; 
          float x, y, z, dx, dy, dz;

          i = pi[ibuf].i;
	  /* Turn index i into separate ix, iy, iz indices */ 
          iz = i/(nxp2*nyp2);                          
          iy = (i - iz*nxp2*nyp2)/nxp2;          
          ix = i - nxp2*(iy+nyp2*iz);                  
          /* Compute real particle position from relative coords and grid data */ 
          /* Store data first in tmp particle struct po so we can apply filter */ 
          po.x = gd.x0+((ix-1)+(pi[ibuf].dx+1)*0.5)*gd.dx;
          po.y = gd.y0+((iy-1)+(pi[ibuf].dy+1)*0.5)*gd.dy;
          po.z = gd.z0+((iz-1)+(pi[ibuf].dz+1)*0.5)*gd.dz;
          po.ux = pi[ibuf].ux;  
          po.uy = pi[ibuf].uy;  
          po.uz = pi[ibuf].uz; 
          po.q  = pi[ibuf].q; 
          if ( FILTER ) {
            nparticles_outfile++; 
            nparticles_output_from_this_file++; 
            memcpy( &po_buf[obuf], &po, sizeof(po) ); 
#if 0     /* Commented this out since memcpy should be more efficient */        
            po_buf[obuf].x=po.x;  
            po_buf[obuf].y=po.y;  
            po_buf[obuf].z=po.z; 
            po_buf[obuf].ux=po.ux;
            po_buf[obuf].uy=po.uy; 
            po_buf[obuf].uz=po.uz; 
            po_buf[obuf].q=po.q; 
#endif 
            obuf++; 
            if ( obuf==PBUF_SIZE ) {
              fwrite( &po_buf, sizeof(struct opdata), PBUF_SIZE, ofp);
              obuf=0; 
            }
          }
        }
      }
      /* Write remaining output buffer to file */ 
      if ( obuf!=0 ) fwrite( &po_buf, sizeof(struct opdata), obuf, ofp);
    }
    printf( "Local particles: %d  Local particles written: %d\n", 
            nparticles_this_file, nparticles_output_from_this_file); 
  }
  printf( "Done processing particles.\n"); 
  printf( "Number of particles processed: %d\nNumber of particles written to file: %d\n", 
          nparticles_total, nparticles_outfile);
  fclose(ofp); 

  /* Close input and output data files */ 
  foreach_file fclose(fparr[ifile]); 

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



