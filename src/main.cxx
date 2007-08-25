/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised from earlier V4PIC versions
 *
 */

#include <vpic.hxx>

# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
#include <fenv.h>

#ifndef PPE_ROUNDING_MODE
#define PPE_ROUNDING_MODE FE_TONEAREST
#endif

# endif

int
main( int argc,
      char **argv ) {
  int m, n;

# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)

  // set PPU rounding mode
  fesetround(PPE_ROUNDING_MODE);

  // Allow processing of SPU-accelerated pipeline workloads on the 8 SPUs

  // Strip threads-per-process arguments from the argument list

  int spp = 8;
  for( m=n=0; n<argc; n++ )
    if( strncmp( argv[n], "-spp=", 5 )==0 ) spp = atoi( argv[n]+5 );
    else                                    argv[m++] = argv[n];
  argv[m] = NULL; // ANSI - argv is NULL terminated
  argc = m;

  spu.boot( spp, // Total number of SPUs for processing pipeline workloads
            0 ); // This PPU thread physically cannot process SPU workloads!

# endif

  // Strip threads-per-process arguments from the argument list

  int tpp = 1;
  for( m=n=0; n<argc; n++ )
    if( strncmp( argv[n], "-tpp=", 5 )==0 ) tpp = atoi( argv[n]+5 );
    else                                    argv[m++] = argv[n];
  argv[m] = NULL; // ANSI - argv is NULL terminated
  argc = m;

  thread.boot( tpp, 1 );
  serial.boot( tpp, 1 );

  // Note: Some MPIs will bind threads to cores if threads are booted
  // after MPI is initialized.  So we start up the pipeline
  // dispatchers _before_ starting up MPI.

  mp_init(argc, argv);

  vpic_simulation simulation; 
  if( argc>=3 && strcmp(argv[1],"restart")==0 ) simulation.restart(argv[2]);
  else simulation.initialize(argc,argv);
  
  // Allow us to change a few run variables such as quota, num_step
  // "on the fly"

  if ( argc==4 ) simulation.modify_runparams( argv[3] );  

  while( simulation.advance() ); 

  // Let all processors finish up
  
  simulation.barrier();
  
  // Issue a termination message when we exit cleanly.
  
  if( simulation.rank()==0 ) 
    MESSAGE(( "Maximum number of time steps reached.  Job has completed." )); 

  mp_finalize( simulation.grid_mp() );

  serial.halt();
  thread.halt();

# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  spu.halt();
# endif

  return 0;
}

