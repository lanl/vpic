#include "util.h"

double _boot_timestamp = 0;

double
uptime( void ) {
  double local_time = wallclock(), time_sum;
  mp_allsum_d( &local_time, &time_sum, 1 );
  return time_sum/(double)world_size - _boot_timestamp;
}

void
boot_services( int * pargc,
               char *** pargv ) {

  // Start up the checkpointing service.  This should be first.

  boot_checkpt( pargc, pargv );

  // Start up the threads.  Note that some MPIs will bind threads to
  // cores if threads are booted _after_ MPI is initialized.  So we
  // start up the pipeline dispatchers _before_ starting up MPI.

  // FIXME: The thread utilities should take responsibility for
  // thread-core affinity instead of leaving this to chance.

  serial.boot( pargc, pargv );
  thread.boot( pargc, pargv );

  // Boot up the communications layer
  // See note above about thread-core-affinity

  boot_mp( pargc, pargv );

  // Set the boot_timestamp
  
  mp_barrier();
  _boot_timestamp = 0;
  _boot_timestamp = uptime();
}

// This operates in reverse order from boot_services

void
halt_services( void ) {
  _boot_timestamp = 0;
  halt_mp();
  thread.halt();
  serial.halt();
  halt_checkpt();
}

