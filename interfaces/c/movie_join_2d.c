/* 
   Utility to join together movie data files (written one file per 
   processor) for 2D diagnostic output from vpic-3.  Currently, 
   the output is for ex and ey (normalized to incident intensity). 
   
   Output data files:  movie_ex.bin, movie_ey.bin

   Usage:  movie_join [num_processors] [binary write flag] [toplogy]

   num_processors    = number of processors used to make the data.
   binary_write_flag = "b" for binary files, "a" for ascii files.
   topology = 0 for decomp in x direction 1 for decomp in z direction

   Example.  Suppose one has 8 files of each type named "movie_ex.0"...
   "movie_ex.7", "movie_ey.0"..."movie_ey.7"... from a run with 
   decomposition in the x direction. 
   
   These are joined using the command from the command line:

     movie_join 8 movie_ex 0
     movie_join 8 movie_ey 0

   or, optionally (since x decomposition is assumed if one leaves off
   the toplogy field),

     movie_join 8 movie_ex
     movie_join 8 movie_ey

   The Ex, Ey output are unformatted binary files of the format: 

   Ex[0,0,0], Ex[0,0,1], ...  Ex[NX-1,NY-1,NZ-1]  <--frame 1 data
   
   ...

   Ex[0,0,0], Ex[0,0,1], ...  Ex[NX-1,NY-1,NZ-1]  <--frame NF data

   The movie_log.dat file is ascii formatted with entries: 
   NX                  total number of x cells in phase frames
   NY                  total number of y cells in phase frames 
   NZ                  total number of z cells in phase frames
   NF                  total number of frames in movie files

   Note:  These data can be read in using the IDL facility movie_2d_data.pro
   ------------------------------------------------------------------------
   Last modified: Brian Albright, X-1, 4/25/2005. 
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
      fprintf( stream, "\n" );                   \
    }
#define SUFFIX(flag) ((flag) ? "bin" : "dat")


/*************************************************************************
  Main routine                                   
**************************************************************************/

int main( int argc, char *argv[] ) {
  int nproc, ifile, nframes;
  FILE **fparr, *ofp, *ologp; 
  int nx_total, ny_total, nz_total, irecord_length, orecord_length; 
  float *idata_vec, *odata_vec; 
  char fname[256], errmsg[512], ofname[256];
  int write_binary=1, decomp=0; 
  int nx, ny, nz; 

  if ( argc!=3 && argc!=4 ) {
    fprintf( stderr, "Usage:  %s [num_processors] [filename]"
	     " [0, 1, 2 for x, y, z decomp]\n\n",
             argv[0] ); 
    return 0;
  } 
  
  /* Obtain runtime data from comand line. */ 
  nproc=atoi(argv[1]); 
  if ( nproc<=0 ) ERROR( "Bad number of input files requested.\n" );
  if ( argc==4 ) {
    decomp = atoi(argv[3]); 
    if ( decomp<0 || decomp>2 ) {
      fprintf(stderr, "Warning: unknown topology; assuming x decomposition."); 
      decomp=0; 
    }
  }

  /* Allocate space for input file handles. */ 
  fparr=emalloc((size_t)nproc*sizeof(*fparr));
  foreach_file fparr[ifile]=emalloc(sizeof(**fparr)); 

  /* Open all of the input files.  */
  foreach_file {
    sprintf(fname, "%s.%d", argv[2], ifile);
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) ERROR(("Cannot open file: %s\n", fname));
    fread( &nx, sizeof(int), 1, fparr[ifile] );
    fread( &ny, sizeof(int), 1, fparr[ifile] ); 
    fread( &nz, sizeof(int), 1, fparr[ifile] ); 
    irecord_length=nx*ny*nz; 
  }

  if      ( decomp==0 ) { nx_total=nx*nproc; ny_total=ny; nz_total=nz; }
  else if ( decomp==1 ) { ny_total=ny*nproc; nz_total=nz; nx_total=nx; }
  else                  { nz_total=nz*nproc; nx_total=nx; ny_total=ny; }
  orecord_length = nx_total*ny_total*nz_total; 
  printf("Total x, y, z points in source data: %d %d %d\n", nx_total, ny_total, nz_total);
  printf("Input, output record lengths: %d %d\n", irecord_length, orecord_length); 
  if ( irecord_length<=0 ) ierror("Cannot proceed. Input record length <= 0"); 
  if ( orecord_length<=0 ) ierror("Cannot proceed. Output record length <= 0"); 
  idata_vec = (float *)emalloc((size_t)(irecord_length)*sizeof(float));
  odata_vec = (float *)emalloc((size_t)(orecord_length)*sizeof(float)); 
  
  /* Prepare and write movie data file. */
  sprintf(ofname, "%s.%s", argv[2], SUFFIX(write_binary));
  if ( !(ofp=fopen(ofname,WRITE_MODE(write_binary) )) ) ERROR(("Cannot open file: %s\n", ofname));
  nframes=0;
  while ( !feof(fparr[0]) ) {
    printf("Reading frame %d... ", ++nframes); 
    foreach_file {
      int ix, ixmin, iy, iymin, iz, izmin; 
      fread( idata_vec, sizeof(float), irecord_length, fparr[ifile] );
      if ( feof(fparr[ifile]) ) {
        printf("end of file reached.  Frame %d not written.\n", nframes--); 
	goto DoneReading; 
      }
      if      ( decomp==0 ) { ixmin=ifile*nx; iymin=0; izmin=0; }
      else if ( decomp==1 ) { iymin=ifile*ny; izmin=0; ixmin=0; }
      else                  { izmin=ifile*nz; ixmin=0; iymin=0; }
      for ( ix=0; ix<nx; ++ix ) 
        for ( iy=0; iy<ny; ++iy )
	  for ( iz=0; iz<nz; ++iz ) 
	    odata_vec[(ix+ixmin)*ny_total*nz_total+(iy+iymin)*nz_total+(iz+izmin)]
	      = idata_vec[ix*ny*nz+iy*nz+iz];
    }
    WRITE_VEC(ofp, odata_vec, orecord_length, sizeof(float), write_binary);
    printf("done.\n");
  } 

 DoneReading: 
  fclose(ofp); 

  /* Close input data files */
  foreach_file fclose(fparr[ifile]);
  free(idata_vec);
  free(odata_vec);

  /* Write movie_log.dat file with frame data */ 
  if ( !(ologp=fopen( "movie_log.dat", "w" )) )  ERROR(("Cannot open log file."));
  fprintf( ologp, "%d\n", nx_total ); 
  fprintf( ologp, "%d\n", ny_total );
  fprintf( ologp, "%d\n", nz_total ); 
  fprintf( ologp, "%d\n", nframes );
  fclose(ologp);
  
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



