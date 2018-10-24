#include "pipelines.h"

#if defined(VPIC_USE_OPENMP)

#include <omp.h>

// omp_boot - mimics thread_boot - see 'pipelines_thread.c'

static void
omp_boot( int * pargc,
	  char *** pargv )
{
  int n_pipeline;

  if ( omp_helper.n_pipeline != 0 ) ERROR(( "OMP container has already booted" ));

  n_pipeline = strip_cmdline_int( pargc, pargv, "--tpp", 1 );

  omp_set_num_threads( n_pipeline );

  if ( n_pipeline < 1 || n_pipeline > MAX_PIPELINE )
    ERROR(( "Invalid number of pipelines requested (%i)", n_pipeline ));

  //initialize dispatch_to_host
  int dispatch_to_host = strip_cmdline_int( pargc, pargv, "--dispatch_to_host", 1 );

  //assign our helper values
  omp_helper.n_pipeline       = n_pipeline;
  omp_helper.dispatch_to_host = dispatch_to_host;
}

/*
static void
omp_util( int * pargc,
          char *** pargv )
{
  const char * f_dump = strip_cmdline_string( pargc, pargv, "--dump-fields",   "none.txt" );
  const char * e_dump = strip_cmdline_string( pargc, pargv, "--dump-energies", "none.txt" );
}
*/

omp_container_t omp_helper = {
  0,
  0,
  omp_boot
};

#endif
