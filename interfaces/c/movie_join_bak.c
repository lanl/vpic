/* 
   Utility to join together the movie data files (written one file per 
   processor) for the 1D diagnostic output from vpic-3 as applied to 
   the 1D LPI problem.  This utility is specifically for the phase space
   arrays:  vx and x binned particle data.  Ex and Ey field arrays can 
   use fft_join to bind them together. 

   Usage:  movie_join [num_processors] [file header] [nx] [binary write flag]

   Example.  Suppose one has 8 files named "movie_ex.0"..."movie_ex.7". 
   These are joined using the command from the command line:

     movie_join 8 movie_ex 200 1

   The binary file "movie_ex.bin" will be created with the output.  

   Notes: 

   1) If we run into troubles with excessively large files, we need to rewrite
   in order to use the ISO standard library functions fgetpos() and fsetpos() 
   in stdio.h because with very large movie files, one may not be guaranteed 
   that all file addresses are representable by a long int.

   2) The output file is an unformatted binary "direct access" file of the 
   following format:  Let fe[index] represent the electron distribution with 
   NX points in x, NVX=200 points in v.  The file is written as:

   [number NF of movie frames]  [system length in lambda_de]
   fe[0], fe[1], ... fe[ix*NVX + ivx] ... fe[(NX-1)*NVX + NVX-1]    <--frame 1 data

   ...
   
   fe[0], fe[1], ... fe[ix*NVX + ivx] ... fe[(NX-1)*NVX + NVX-1]    <--frame NF data


   ------------------------------------------------------------------------
   Last modified: Brian Albright, X-1, 3/29/2005. 
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

int main( int argc, char *argv[] ) {
  int nproc, ifile, nframes=0;
  FILE **fparr, *ofp, *ologp; 
  int *nxarr; 
  int nx, nx_output, nx_total, nx_ave; 
  int nvx, nvx_output=200; 
  float dv; 
  float *data_vec, *odata_vec; 
  char fname[256], errmsg[512], ofname[256];
  int write_binary=1; 
  int not_eof=1; 

  if ( argc!=5 ) {
    fprintf( stderr, "Usage:  %s [num_processors] "
	     " [NX] [a or b to write ascii or binary data]\n\n"
             "Note: ascii data write is intended for debugging one's diagnostics. The files"
             "produced can be huge.\n", argv[0] ); 
    return 0;
  } 
  
  /* Obtain runtime data from comand line. */ 
  nproc=atoi(argv[1]); 
  if ( nproc<=0 ) ERROR( "Bad number of input files requested.\n" );
  nx_output=atoi(argv[2]); 
  if ( nx_output<=0 ) ERROR( "Bad nx requested.\n" ); 
  if ( argv[3][0]=='a' || argv[3][0]=='A' ) write_binary=0;

  /* Allocate space for input file handles. */ 
  fparr=emalloc((size_t)nproc*sizeof(*fparr));
  nxarr=emalloc((size_t)nproc*sizeof(*nxarr)); 
  foreach_file fparr[ifile]=emalloc(sizeof(**fparr)); 

  /************************************************************************* 
     Write phase space movie data. 
  **************************************************************************/
 
  /* Open all of the input files. */ 
  foreach_file {
    sprintf(fname, "movie_phase.%d", ifile); 
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) ERROR(("Cannot open file: %s\n", fname)); 
    /* First element in each file is the num of cells in local x mesh */
    fread( &nxarr[ifile], sizeof(int), 1, fparr[ifile] ); 
    fread( &nvx, sizeof(int), 1, fparr[ifile] ); 
    fread( &dv, sizeof(float), 1, fparr[ifile] ); 
  }
  
  nx_total = 0; 
  foreach_file nx_total += nxarr[ifile]; 
  printf("Total x mesh points: %d\n", nx_total); 
  nx_ave = nx_total/nx_output; 
  if ( nx_total%nx_output ) 
    fprintf(stderr, "Warning:  requested number of x points doesn't "
	    "evenly divide nx_total.\n" ); 
  data_vec = (float *)emalloc((size_t)(nx_total*nvx_output)*sizeof(*data_vec)); 
  odata_vec = (float *)emalloc((size_t)(nx_output*nvx_output)*sizeof(*odata_vec)); 

  /* Write various movie parameters to log file movie_log.dat */ 
  if ( !(ologp=fopen( "movie_log.dat", "w" ))  )  ERROR(("Cannot open log file."));
  fprintf( ologp, "%d\n", nframes );   
  fprintf( ologp, "%d\n", nx_output );    
  fprintf( ologp, "%d\n", nvx_output ); 
  fprintf( ologp, "%e\n", dv ); 

  /* Prepare to write movie data file. */ 
  sprintf(ofname, (write_binary ? "%s.bin" : "%s.dat"), argv[2]);
  if ( !(ofp = fopen( ofname, WRITE_MODE(write_binary) )) ) 
    ERROR(("Cannot open file: %s\n", ofname));         
  while ( not_eof ) {
    int ix=0;
    /* Read one line of data from each input file and join together into 
       data_vec array.  We have assumed that each output file has the same 
       number of lines in it. */ 
    printf("Reading frame %d... ", ++nframes); 
    foreach_file {
      fread( &data_vec[ix], sizeof(float), nxarr[ifile]*nvx_output, fparr[ifile] ); 
      ix += nxarr[ifile]*nvx_output; 
      if ( feof(fparr[ifile]) ) not_eof=0; 
    }

    /* Now create the output data array by averaging over values in data_vec */ 
    for ( ix=0; ix<nx_output; ++ix ) {
      int ivx;
      for ( ivx=0; ivx<nvx_output; ++ivx ) {
	int ixa;
	odata_vec[ix*nvx_output+ivx] = 0;
	for ( ixa=0; ixa<nx_ave; ++ixa )
	  if ( ix*nx_ave+ixa<nx_total )
	    odata_vec[ix*nvx_output+ivx] 
	      += data_vec[(ix*nx_ave+ixa)*nvx_output+ivx]; 
      }
    }

    /* Write odata_vec into output stream */ 
    WRITE_VEC(ofp, odata_vec,nx_output*nvx_output,sizeof(float),write_binary);
#if 0
    if ( write_binary )
      fwrite( odata_vec, sizeof(float), nx_output*nvx_output, ofp );
    else  /* ascii write */ 
      for ( ix=0; ix<nx_output; ++ix ) {
        int ivx;
        for ( ivx=0; ivx<nvx_output; ++ivx ) 
          fprintf( ofp, "%e ", odata_vec[ix*nvx_output+ivx] );
        fprintf( ofp, "\n" );  
      }
#endif
    printf("done.\n"); 
  }  
  fprintf( ologp, "%d\n", nframes );   
  fclose( ologp ); 

  
  /* Close input and output data files */ 
  foreach_file fclose(fparr[ifile]); 
  fclose(ofp); 


  /*************************************************************************  
     Write field movie data.  
  **************************************************************************/

  
  return 0;
}


/* 
   Print an error message to stderr and exit program. 
*/ 

void ierror( char *msg ) {
  fprintf( stderr, msg ); 
  exit(1); 
}


/* 
   Allocate space with an error check that aborts if malloc is
   unsuccessful.  Returns a pointer to allocated space. 
*/ 

void *emalloc( size_t size ) {
  void *ptr = NULL; 

  if ( size<0 )
    ierror( "emalloc: request for array of negative size.\n" ); 
  if ( !(ptr=malloc(size)) )
    ierror( "emalloc: cannot allocate space.\n" ); 
  return ptr; 
}


