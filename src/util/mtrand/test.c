#include <stdio.h>
#include <stdlib.h>
#include "mtrand.h"

int
main( int argc,
      char **argv ) {

  const char * fname = argv[1];
  const int N = atoi( argv[2] );
  float * x;
  mt_rng_t * rng;
  FILE * file;

  x = malloc( N*sizeof(x[0]) );

  rng = new_mt_rng( atoi( argv[3] ) );
  mt_frandn_fill( rng, x, N );
  delete_mt_rng( rng );

  file = fopen( fname, "w" );
  fwrite( x, sizeof(x[0]), N, file );
  fclose( file );

  free( x );

  return 0;
}
