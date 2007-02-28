/* 
   Utility to join together the poynting flux data files (written one file 
   per processor) for the 2D diagnostic output from vpic-3 as applied to 
   the 2D LPI problem.  

   Usage:  poynting2d [num_processors] [binary write flag]

   num_processors    = number of processors used to make the data.
   binary_write_flag = "b" for binary files, "a" for ascii files.

   The output files are unformatted binary "direct access" files.

   Format:  for top, bottom: sequence of (length nx arrays of float)
            for left, right: sequence of (length nz arrays of float) 
   ------------------------------------------------------------------------
   Brian Albright, X-1, 6/18/2005. 
*/ 


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

void ierror( char *msg ); 
void *emalloc( size_t size ); 
void store_frame_data( FILE *fp, int frames );

#define ERROR(m) {     \
  char errmsg[512];    \
  sprintf(errmsg, m);  \
  ierror(errmsg);      \
}

#define foreach_file    for ( ifile=0; ifile<nproc; ++ifile )

#define WRITE_MODE(flag)  ( (flag) ? "wb" : "w" )
#define WRITE_VEC(stream,data,length,size,flag)  \
    if ( flag )                                  \
      fwrite( data, size, length, stream );      \
    else  { /* ascii write */                    \
      int itmp;                                  \
      for ( itmp=0; itmp<length; ++itmp )        \
        fprintf( stream, "%e ", data[itmp] );    \
     /* fprintf( stream, "\n" ); */              \
    }
#define SUFFIX(flag) ((flag) ? "bin" : "dat")


/*************************************************************************
  Main routine                                   
**************************************************************************/

int main( int argc, char *argv[] ) {
  int nproc, ifile, nframes;
  FILE **fparr, *olog, *oftop, *ofbottom, *ofleft, *ofright; 
  int nx, nz; 
  float *data_vec_x, *data_vec_z; 
  float *left, *right, *top, *bottom; 
  char fname[256], errmsg[512];
  char ofname_top[256], ofname_bottom[256], ofname_left[256], ofname_right[256];
  int write_binary=1; 
  int not_eof; 

  if ( argc!=3 ) {
    fprintf( stderr, "Usage:  %s [num_processors] "
	     " [a or b to write ascii or binary data]\n\n",
             argv[0] ); 
    return 0;
  } 
  
  /* Obtain runtime data from comand line. */ 
  nproc=atoi(argv[1]); 
  if ( nproc<=0 ) ERROR( "Bad number of input files requested.\n" );
  if ( argv[2][0]=='a' || argv[2][0]=='A' ) write_binary=0;

  /* Allocate space for input file handles. */ 
  fparr=emalloc((size_t)nproc*sizeof(*fparr));
  foreach_file fparr[ifile]=emalloc(sizeof(**fparr)); 

  /* Prepare output files */ 
# define PREPARE_OUTPUT_FILE(DIR, DIRSTR)                                 \
  sprintf(ofname_##DIR, "poynting_%s.%s", DIRSTR, SUFFIX(write_binary));  \
  if ( !(of##DIR = fopen( ofname_##DIR, WRITE_MODE(write_binary) )) )     \
    ERROR(("Cannot open file: %s\n", ofname_##DIR));           
  PREPARE_OUTPUT_FILE(top,    "top");
  PREPARE_OUTPUT_FILE(bottom, "bottom");
  PREPARE_OUTPUT_FILE(left,   "left");
  PREPARE_OUTPUT_FILE(right,  "right");

  /* Open all of the input files. */ 
  foreach_file {
    sprintf(fname, "poynting.%d", ifile); 
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) {
      fprintf(stderr, "Cannot open file: %s\n", fname); 
      exit(1);
    }
    /* First two elements in each file are the number of nx, nz points */
    fread( &nx, sizeof(int), 1, fparr[ifile] ); 
    fread( &nz, sizeof(int), 1, fparr[ifile] ); 
  }
  /* Assume 1D domain decomposition in x */ 
  printf("Total x mesh points in phase source data: %d\n", nx*nproc); 
  printf("Total z mesh points in phase source data: %d\n", nz); 
  data_vec_x = (float *)emalloc(2*(size_t)nx*sizeof(float)); 
  data_vec_z = (float *)emalloc(2*(size_t)nz*sizeof(float)); 
  left = data_vec_z;
  right = data_vec_z+nz;
  bottom = data_vec_x;
  top = data_vec_x+nx; 

  /* Prepare to write data files. */ 
  not_eof=1;
  nframes=0;
  while ( not_eof ) {
    int ix=0;
    /* For each input file, read one poynting data block.  Then write the 
       data to the appropriate places in the output files. */  
    printf("Processing frame %i... ", ++nframes); 
    foreach_file {
      fread( data_vec_z, sizeof(float), 2*nz, fparr[ifile] ); 
      fread( data_vec_x, sizeof(float), 2*nx, fparr[ifile] ); 
      if ( feof(fparr[ifile]) ) not_eof=0; 
      if ( ifile==0 ) {
        WRITE_VEC(ofleft,  left,   nz, sizeof(float), write_binary);     
        if ( !write_binary ) fprintf(ofleft, "\n"); 
      }
      if ( ifile==nproc-1 ) {     
        WRITE_VEC(ofright, right,  nz, sizeof(float), write_binary);     
	if ( !write_binary ) fprintf(ofright, "\n"); 
      }
      WRITE_VEC(ofbottom,  bottom, nx, sizeof(float), write_binary);     
      WRITE_VEC(oftop,     top,    nx, sizeof(float), write_binary);     
    }
    if ( !write_binary ) fprintf(ofbottom, "\n"); 
    if ( !write_binary ) fprintf(oftop, "\n"); 
    printf("done.\n"); 
  }  
  ++nframes; 
  free(data_vec_x); 
  free(data_vec_z); 

  /* Close input and output data files */ 
  foreach_file fclose(fparr[ifile]); 
  fclose( oftop );
  fclose( ofbottom );
  fclose( ofleft );
  fclose( ofright );

  /* Write various parameters to log file poynting.dat */
  if ( !(olog=fopen( "poynting.dat", "w" ))  )  ERROR(("Cannot open log file."));
  fprintf( olog, "%d\n", nx*nproc ); 
  fprintf( olog, "%d\n", nz );
  fprintf( olog, "%d\n", nframes );
  fclose(olog);
 
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


