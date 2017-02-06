/*
   Utility to combine data files for making fig 10 from Wilks et al
   PoP 2001.

   Usage:  fig10  [num processors] [frame] [a or b for ascii/binary write]

   ------------------------------------------------------------------------
   Last modified: Brian Albright, X-1, 4/22/2005.
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
  int *nxarr;
  int nux, nuz, nframe;
  int nx, nx_output, nx_total, nx_ave;
  int nvx, nvx_output=200;
  float dv;
  double *data_vec, *odata_vec;
  char fname[256], errmsg[512], ofname[256];
  int write_binary=1;
  int not_eof;

  if ( argc!=4 ) {
    fprintf( stderr, "Usage:  %s [num_processors] "
	     " [frame] [a or b to write ascii or binary data]\n\n",
             argv[0] );
    return 0;
  }

  /* Obtain runtime data from comand line. */
  nproc=atoi(argv[1]);
  if ( nproc<=0 ) ERROR( "Bad number of input files requested.\n" );
  nframe=atoi(argv[2]);
  if ( nframe<=0 ) ERROR( "Bad nx requested.\n" );
  if ( argv[3][0]=='a' || argv[3][0]=='A' ) write_binary=0;

  /* Allocate space for input file handles. */
  fparr=emalloc((size_t)nproc*sizeof(*fparr));
  nxarr=emalloc((size_t)nproc*sizeof(*nxarr));
  foreach_file fparr[ifile]=emalloc(sizeof(**fparr));

  /*************************************************************************
     Write electron phase space movie data.
  **************************************************************************/

  /* Open all of the input files. */
  foreach_file {
    sprintf(fname, "fig10_H.%d", ifile);
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) ERROR(("Cannot open phase space files.\n"));
    /* First element in each file is nux, second is nuz */
    fread( &nux, sizeof(int), 1, fparr[ifile] );
    fread( &nuz, sizeof(int), 1, fparr[ifile] );
  }

  printf("Total nux, nuz points in source data: %d %d\n", nux, nuz);
  printf("Frame requested: %d\n", nframe);
  data_vec  = (double *)emalloc((size_t)(nux*nuz)*sizeof(*data_vec));
  odata_vec = (double *)emalloc((size_t)(nux*nuz)*sizeof(*odata_vec));

  /* Prepare to write movie data file. */
  sprintf(ofname, "fig10.%s", SUFFIX(write_binary));
  if ( !(ofp = fopen( ofname, WRITE_MODE(write_binary) )) ) ERROR(("Cannot open file: %s\n", ofname));
  {
    int iux, iuz;

    for ( iux=0; iux<nux; ++iux )
      for ( iuz=0; iuz<nuz; ++iuz )
        odata_vec[iux*nuz+iuz]=0;
    foreach_file {
      fseek( fparr[ifile], 2*sizeof(int)+nframe*nux*nuz*sizeof(double), SEEK_SET );
      fread( data_vec, sizeof(double), nux*nuz, fparr[ifile] );
      for ( iux=0; iux<nux; ++iux )
        for ( iuz=0; iuz<nuz; ++iuz )
	  odata_vec[iux*nuz+iuz]+=data_vec[iux*nuz+iuz];
    }
  }

  /* Write odata_vec into output stream */
  WRITE_VEC(ofp, odata_vec, nux*nuz, sizeof(double), write_binary);
  printf("done.\n");
  fclose(ofp);

  /* Close input and output data files */
  foreach_file fclose(fparr[ifile]);
  free(data_vec);
  free(odata_vec);

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

  if ( size<0 ) ierror( "emalloc: request for array of negative size.\n" );
  if ( !(ptr=malloc(size)) ) ierror( "emalloc: cannot allocate space.\n" );
  return ptr;
}
