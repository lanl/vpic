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

int main( int argc, char **argv ) {

  /* Note: Some MPIs will bind threads to cores if threads are booted
     after MPI is initialized.  So we start up the pipeline dispatchers
     _before_ starting up MPI. */

  serial.boot(1);
  MPI_Init(&argc,&argv); 

  vpic_simulation simulation; 
  if( argc>=3 && strcmp(argv[1],"restart")==0 ) simulation.restart(argv[2]);
  else simulation.initialize(argc,argv);
  
  /* Allow us to change a few run variables such as quota, num_step "on the fly" */ 
  if ( argc==4 ) simulation.modify_runparams( argv[3] );  

  while(simulation.advance()); 

  /* Issue a termination message when we exit cleanly. */ 
  {
    int rank; 
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if ( rank==0 ) 
      MESSAGE(("Maximum number of time steps reached.  Job has completed.")); 
  }

  MPI_Finalize();
  serial.halt();

  return 0;
}

