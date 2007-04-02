/* 
   Utility to join together the fft data files (written one file per 
   processor) for the 1D diagnostic output from vpic-3 as applied to 
   the 1D LPI problem.  

   Usage:  fft_join [num_processors] [file header]

   Example.  Suppose one has 8 files named "fft_ex.0"..."fft_ex.7". 
   These are joined using the command from the command line:

     fft_join 8 fft_ex

   The file "fft_ex.bin" will be created with the output. 

   ------------------------------------------------------------------------
   Last modified: Brian Albright, X-1, 3/29/2005. 
*/ 


#include <stdio.h>
#include <stdlib.h>

void ierror( char *msg ); 
void *emalloc( size_t size ); 

#define ERROR(m) {     \
  char errmsg[512];    \
  sprintf(errmsg, m);  \
  ierror(errmsg);      \
}

#define foreach_file    for ( ifile=0; ifile<nproc; ++ifile )

int main( int argc, char *argv[] ) {
  int nproc, ifile, nx, nx_total=0, nlines=0;
  FILE **fparr, *ofp; 
  int *nxarr; 
  float *data_vec; 
  char fname[256], errmsg[512], ofname[256];
  int not_eof=1; 

  if ( argc!=3 ) {
    fprintf( stderr, "Usage:  %s [num_processors] [filename_header]\n", argv[0] ); 
    return 0;
  } 
  
  /* Obtain number of data files from command line. */ 
  nproc=atoi(argv[1]); 
  if ( nproc<=0 ) ERROR( "Bad number of input files requested.\n" );

  /* Allocate space for input file handles. */ 
  fparr=(FILE **)emalloc((size_t)nproc*sizeof(*fparr));
  nxarr=(int *)emalloc((size_t)nproc*sizeof(*nxarr)); 
  foreach_file fparr[ifile]=emalloc(sizeof(**fparr)); 
  
  /* Open all of the input files. */ 
  foreach_file {
    sprintf(fname, "%s.%d", argv[2], ifile); 
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) ERROR(("Cannot open file: %s\n", fname)); 
    /* First element in each file is the num of cells in local x mesh */
    fread( &nxarr[ifile], sizeof(int), 1, fparr[ifile] ); 
  }
  
  foreach_file nx_total += nxarr[ifile]; 
  printf("Total x mesh points: %d\n", nx_total); 
  data_vec = (float *)emalloc((size_t)nx_total*sizeof(*data_vec)); 
  sprintf(ofname, "%s.bin", argv[2]); 
  if ( !(ofp=fopen(ofname, "wb")) ) ERROR(("Cannot open file: %s\n", ofname)); 
  while ( not_eof ) {
    int ix=0;

    /* Read one line of data from each input file and join together into 
       data_vec array.  We have assumed that each output file will have the 
       same number of lines in it. */ 
    printf("Reading line %d... ", ++nlines); 
    foreach_file {
      fread( &data_vec[ix], sizeof(float), nxarr[ifile], fparr[ifile] ); 
      ix += nxarr[ifile]; 
      if ( feof(fparr[ifile]) ) not_eof=0; 
    }

    /* Write data_vec into output stream */ 
    fwrite( data_vec, sizeof(float), nx_total, ofp );
    printf("done.\n"); 
  }  
  
  /* Close all files */ 
  foreach_file fclose(fparr[ifile]); 
  fclose(ofp); 
  
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
