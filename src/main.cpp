/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised from earlier V4PIC versions
 *
 */

#if defined ENABLE_DMP
#include <mpi.h>    /* For MPI_Init and MPI_Finalize */
#else
#include <mpi_stubs.h>
#endif // ENABLE_DMP

#include <string.h> /* For strcmp */
#include <vpic.hpp>

int main( int argc, char **argv ) {
  MPI_Init(&argc,&argv); 

  /* Do some hardware compatibility checking */
  if( sizeof(INT32_TYPE)!=4 ||
      sizeof(INT64_TYPE)!=8 ) {
    ERROR(( "INT32_TYPE and/or INT64_TYPE have turned this code into a house of lies." ));
    MPI_Finalize();
    return 1;
  }

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
  return 0;
}


