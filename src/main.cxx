/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised from earlier V4PIC versions
 *
 */

#include <mpi.h>    /* For MPI_Init and MPI_Finalize ... WRAP THIS! */
#include <vpic.hxx>

#ifndef PIPELINE_BOOT_REQUEST
#define PIPELINE_BOOT_REQUEST 1
#endif

int
main( int argc,
      char **argv ) {

  /* Note: Some MPIs will bind threads to cores if threads are booted
     after MPI is initialized.  So we start up the pipeline dispatchers
     _before_ starting up MPI. */

  PSTYLE.boot( PIPELINE_BOOT_REQUEST );

# ifdef USE_CELL_SPUS
  /* FIXME: n_pipeline rationalization throughout code.
     Also, need to discuss MPI processor model for RoadRunner with Ben
     to determine if each there is a one-to-one correspondence with
     an MPI processor and PPU (hopefully there is) or if there are
     multiple PPU's per MPI process.  It may make some subtle
     differences in PPU threading model and SPU dispatcher. */
  spu.boot( 8 );
# endif

  mp_init(argc, argv);

  vpic_simulation simulation; 
  if( argc>=3 && strcmp(argv[1],"restart")==0 ) simulation.restart(argv[2]);
  else simulation.initialize(argc,argv);
  
  /* Allow us to change a few run variables such as quota, num_step "on the fly" */ 
  if ( argc==4 ) simulation.modify_runparams( argv[3] );  

  while(simulation.advance()); 

  /* Let all processors finish up */

  simulation.barrier();

  /* Issue a termination message when we exit cleanly. */ 

  {
    int rank = (int)simulation.rank();
    if ( rank==0 ) 
      MESSAGE(("Maximum number of time steps reached.  Job has completed.")); 
  }

  mp_finalize();

  PSTYLE.halt();

# ifdef USE_CELL_SPUS
  spu.halt();
# endif

  return 0;
}

