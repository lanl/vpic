/* 
   Utility to join together the movie data files (written one file per 
   processor) for the 1D diagnostic output from vpic-3 as applied to 
   the 1D LPI problem.  This utility is for both the phase space
   arrays of vx and x binned particle data as well as Ex and Ey field 
   arrays; all are joined in this single utility into output movie files: 

    movie_phase.bin, movie_ex.bin, movie_ey.bin

   Usage:  movie_join [num_processors] [nx] [binary write flag]

   num_processors    = number of processors used to make the data.
   nx                = number of desired x-mesh cells in the phase space  
                       movie files (bins are averaged to this value). 
   binary_write_flag = "b" for binary files, "a" for ascii files.

   Example.  Suppose one has 8 files of each type named "movie_ex.0"...
   "movie_ex.7", "movie_ey.0"..."movie_ey.7", "movie_phase.0"...
   "movie_phase.7".   
   
   These are joined using the command from the command line:

     movie_join 8 200 b

   The phase output file is an unformatted binary "direct access" file of the 
   following format:  Let fe[index] represent the electron distribution with 
   NX points in x, NVX=200 points in v.  The file is written as:

   fe[0], fe[1], ... fe[ix*NVX + ivx] ... fe[(NX-1)*NVX + NVX-1]    <--frame 1 data

   ...
   
   fe[0], fe[1], ... fe[ix*NVX + ivx] ... fe[(NX-1)*NVX + NVX-1]    <--frame NF data

   The Ex, Ey field files are unformatted binary of the format: 

   Ex[0], Ex[1], ...  Ex[NX_TOTAL]  <--frame 1 data

   ...
 
   Ex[0], Ex[1], ...  EX[NX_TOTAL]  <--frame NF data


   The movie_log.dat file is ascii formatted with entries: 
   NX_TOTAL            total number of cells in Ex, Ey frames
   NX                  total number of x cells in phase frames
   NVX                 total number of vx cells in phase frames
   DV                  velocity space interval 
   NF                  total number of frames in movie files

   ------------------------------------------------------------------------
   Last modified: Brian Albright, X-1, 4/3/2005. 
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
  int nx, nx_output, nx_total, nx_ave; 
  int nvx, nvx_output=200; 
  float dv; 
  float *data_vec, *odata_vec; 
  char fname[256], errmsg[512], ofname[256];
  int write_binary=1; 
  int not_eof; 

  if ( argc!=4 ) {
    fprintf( stderr, "Usage:  %s [num_processors] "
	     " [NX] [a or b to write ascii or binary data]\n\n",
             argv[0] ); 
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
     Write electron phase space movie data. 
  **************************************************************************/
 
  /* Open all of the input files. */ 
  foreach_file {
    sprintf(fname, "movie_phase_e.%d", ifile); 
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) {
      printf("Cannot open electron phase space files. Skipping ahead to H ions.\n"); 
      goto Process_Hydrogen_ions;
    }
    /* First element in each file is the num of cells in local x mesh */
    fread( &nxarr[ifile], sizeof(int), 1, fparr[ifile] ); 
    fread( &nvx, sizeof(int), 1, fparr[ifile] ); 
    fread( &dv, sizeof(float), 1, fparr[ifile] ); 
  }
  
  nx_total = 0; 
  foreach_file nx_total += nxarr[ifile]; 
  printf("Total x mesh points in phase source data: %d\n", nx_total); 
  nx_ave = nx_total/nx_output; 
  if ( nx_total%nx_output ) 
    fprintf(stderr, "Warning:  requested number of x points doesn't "
	    "evenly divide nx_total.\n" ); 
  data_vec = (float *)emalloc((size_t)(nx_total*nvx_output)*sizeof(*data_vec)); 
  odata_vec = (float *)emalloc((size_t)(nx_output*nvx_output)*sizeof(*odata_vec)); 

  /* Prepare to write movie data file. */ 
  sprintf(ofname, "movie_phase_e.%s", SUFFIX(write_binary));
  if ( !(ofp = fopen( ofname, WRITE_MODE(write_binary) )) ) 
    ERROR(("Cannot open file: %s\n", ofname));         
  not_eof=1;
  nframes=0;
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
	for ( ixa=0; ixa<nx_ave && ix*nx_ave+ixa<nx_total; ++ixa )
	  odata_vec[ix*nvx_output+ivx] += data_vec[(ix*nx_ave+ixa)*nvx_output+ivx];
        /* divide by nx_ave unless interval straddles nx_total; then, divide by
           actual number of cells summed. */ 
        odata_vec[ix*nvx_output+ivx] /= ixa; 
      }
    }

    /* Write odata_vec into output stream */ 
    WRITE_VEC(ofp, odata_vec, nx_output*nvx_output, sizeof(float), write_binary);
    printf("done.\n"); 
  }  
  fclose(ofp); 

  /* Close input and output data files */ 
  foreach_file fclose(fparr[ifile]); 
  free(data_vec);
  free(odata_vec);


  /*************************************************************************
     Write H phase space movie data.
  **************************************************************************/

 Process_Hydrogen_ions: 
  /* Open all of the input files. */
  foreach_file {
    sprintf(fname, "movie_phase_H.%d", ifile);
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) {
      printf("Cannot open H phase space files. Skipping ahead to He ions.\n");       
      goto Process_Helium_ions; 
    }
    /* First element in each file is the num of cells in local x mesh */
    fread( &nxarr[ifile], sizeof(int), 1, fparr[ifile] );
    fread( &nvx, sizeof(int), 1, fparr[ifile] );
    fread( &dv, sizeof(float), 1, fparr[ifile] );
  }

  nx_total = 0;
  foreach_file nx_total += nxarr[ifile];
  printf("Total x mesh points in phase source data: %d\n", nx_total);
  nx_ave = nx_total/nx_output;
  if ( nx_total%nx_output )
    fprintf(stderr, "Warning:  requested number of x points doesn't "
            "evenly divide nx_total.\n" );
  data_vec = (float *)emalloc((size_t)(nx_total*nvx_output)*sizeof(*data_vec));
  odata_vec = (float *)emalloc((size_t)(nx_output*nvx_output)*sizeof(*odata_vec));

  /* Prepare to write movie data file. */
  sprintf(ofname, "movie_phase_H.%s", SUFFIX(write_binary));
  if ( !(ofp = fopen( ofname, WRITE_MODE(write_binary) )) )
    ERROR(("Cannot open file: %s\n", ofname));
  not_eof=1;
  nframes=0;
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
        for ( ixa=0; ixa<nx_ave && ix*nx_ave+ixa<nx_total; ++ixa )
          odata_vec[ix*nvx_output+ivx] += data_vec[(ix*nx_ave+ixa)*nvx_output+ivx];
        /* divide by nx_ave unless interval straddles nx_total; then, divide by
           actual number of cells summed. */
        odata_vec[ix*nvx_output+ivx] /= ixa;
      }
    }

    /* Write odata_vec into output stream */
    WRITE_VEC(ofp, odata_vec, nx_output*nvx_output, sizeof(float), write_binary);
    printf("done.\n");
  }
  fclose(ofp);

  /* Close input and output data files */
  foreach_file fclose(fparr[ifile]);
  free(data_vec);
  free(odata_vec);


  /*************************************************************************
     Write He phase space movie data.
  **************************************************************************/

 Process_Helium_ions: 
  /* Open all of the input files. */
  foreach_file {
    sprintf(fname, "movie_phase_He.%d", ifile);
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) {
      printf("Cannot open He phase space files. Skipping ahead to fields.\n"); 
      goto Process_fields; 
    }
    /* First element in each file is the num of cells in local x mesh */
    fread( &nxarr[ifile], sizeof(int), 1, fparr[ifile] );
    fread( &nvx, sizeof(int), 1, fparr[ifile] );
    fread( &dv, sizeof(float), 1, fparr[ifile] );
  }

  nx_total = 0;
  foreach_file nx_total += nxarr[ifile];
  printf("Total x mesh points in phase source data: %d\n", nx_total);
  nx_ave = nx_total/nx_output;
  if ( nx_total%nx_output )
    fprintf(stderr, "Warning:  requested number of x points doesn't "
            "evenly divide nx_total.\n" );
  data_vec = (float *)emalloc((size_t)(nx_total*nvx_output)*sizeof(*data_vec));
  odata_vec = (float *)emalloc((size_t)(nx_output*nvx_output)*sizeof(*odata_vec));

  /* Prepare to write movie data file. */
  sprintf(ofname, "movie_phase_He.%s", SUFFIX(write_binary));
  if ( !(ofp = fopen( ofname, WRITE_MODE(write_binary) )) )
    ERROR(("Cannot open file: %s\n", ofname));
  not_eof=1;
  nframes=0;
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
        for ( ixa=0; ixa<nx_ave && ix*nx_ave+ixa<nx_total; ++ixa )
          odata_vec[ix*nvx_output+ivx] += data_vec[(ix*nx_ave+ixa)*nvx_output+ivx];
        /* divide by nx_ave unless interval straddles nx_total; then, divide by
           actual number of cells summed. */
        odata_vec[ix*nvx_output+ivx] /= ixa;
      }
    }

    /* Write odata_vec into output stream */
    WRITE_VEC(ofp, odata_vec, nx_output*nvx_output, sizeof(float), write_binary);
    printf("done.\n");
  }
  fclose(ofp);

  /* Close input and output data files */
  foreach_file fclose(fparr[ifile]);
  free(data_vec);
  free(odata_vec);


  /*************************************************************************  
     Write Ex movie data.  
  **************************************************************************/
 
 Process_fields:
  /* Open all of the input files. */
  foreach_file {
    sprintf(fname, "movie_ex.%d", ifile);
    if ( !(fparr[ifile]=fopen(fname, "rb")) ) 
      ERROR(("Cannot open file: %s\n", fname));
    /* First element in each file is the num of cells in local x mesh */
    fread( &nxarr[ifile], sizeof(int), 1, fparr[ifile] );
  }
 
  nx_total = 0;
  foreach_file nx_total += nxarr[ifile];
  printf("Total x mesh points in Ex source data: %d\n", nx_total);
  data_vec = (float *)emalloc((size_t)(nx_total)*sizeof(*data_vec));
  
  /* Prepare to write movie data file. */
  sprintf(ofname, "movie_ex.%s", SUFFIX(write_binary));
  if ( !(ofp = fopen( ofname, WRITE_MODE(write_binary) )) )
    ERROR(("Cannot open file: %s\n", ofname));
  not_eof=1;
  nframes=0;
  while ( not_eof ) {
    int ix=0;
    /* Read one line of data from each input file and join together into
       data_vec array.  We have assumed that each output file has the same
       number of lines in it. */
    printf("Reading frame %d... ", ++nframes);
    foreach_file {
      fread( &data_vec[ix], sizeof(float), nxarr[ifile], fparr[ifile] );
      ix += nxarr[ifile];
      if ( feof(fparr[ifile]) ) not_eof=0;
    }

    /* Write data_vec into output stream */
    WRITE_VEC(ofp, data_vec, nx_total, sizeof(float), write_binary);
    printf("done.\n");
  } 
  fclose(ofp); 

  /* Close input data files */
  foreach_file fclose(fparr[ifile]);
  free(data_vec);


  /*************************************************************************
     Write Ey movie data.
  **************************************************************************/

  /* Open all of the input files. */
  foreach_file {
    sprintf(fname, "movie_ey.%d", ifile);
    if ( !(fparr[ifile]=fopen(fname, "rb")) )
      ERROR(("Cannot open file: %s\n", fname));
    /* First element in each file is the num of cells in local x mesh */
    fread( &nxarr[ifile], sizeof(int), 1, fparr[ifile] );
  }

  nx_total = 0;
  foreach_file nx_total += nxarr[ifile];
  printf("Total x mesh points in Ey source data: %d\n", nx_total);
  data_vec = (float *)emalloc((size_t)(nx_total)*sizeof(*data_vec));
 
  /* Prepare to write movie data file. */
  sprintf(ofname, "movie_ey.%s", SUFFIX(write_binary));
  if ( !(ofp = fopen( ofname, WRITE_MODE(write_binary) )) )
    ERROR(("Cannot open file: %s\n", ofname));
  not_eof=1;
  nframes=0;
  while ( not_eof ) {
    int ix=0;
    /* Read one line of data from each input file and join together into
       data_vec array.  We have assumed that each output file has the same
       number of lines in it. */
    printf("Reading frame %d... ", ++nframes);
    foreach_file {
      fread( &data_vec[ix], sizeof(float), nxarr[ifile], fparr[ifile] );
      ix += nxarr[ifile];
      if ( feof(fparr[ifile]) ) not_eof=0;
    }

    /* Write data_vec into output stream */
    WRITE_VEC(ofp, data_vec, nx_total, sizeof(float), write_binary);
    printf("done.\n");
  }
  fclose(ofp);

  /* Close input data files */
  foreach_file fclose(fparr[ifile]);
  free(data_vec);

  /************************************************************************
     Write various movie parameters to log file movie_log.dat 
  ************************************************************************/
  if ( !(ologp=fopen( "movie_log.dat", "w" ))  )  ERROR(("Cannot open log file."));
  fprintf( ologp, "%d\n", nx_total ); 
  fprintf( ologp, "%d\n", nx_output );
  fprintf( ologp, "%d\n", nvx_output );
  fprintf( ologp, "%e\n", dv );
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



