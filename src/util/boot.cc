#include "util.h"

#include "stdio.h"

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

  // Boot up the communications layer

#if defined(VPIC_USE_PTHREADS)
  #if defined(VPIC_SWAP_MPI_PTHREAD_INIT)
    boot_mp( pargc, pargv );      // Boot communication layer first.
    serial.boot( pargc, pargv );
    thread.boot( pargc, pargv );
  #else
    serial.boot( pargc, pargv );
    thread.boot( pargc, pargv );
    boot_mp( pargc, pargv );      // Boot communication layer last.
  #endif

#elif defined(VPIC_USE_OPENMP)
  boot_mp( pargc, pargv );        // Boot communication layer first.
  omp_helper.boot( pargc, pargv );

#endif

  // Set the boot_timestamp

  mp_barrier();
  _boot_timestamp = 0;
  _boot_timestamp = uptime();

  if (_world_rank == 0)
  {
      printf("Booting with %d threads and %d (MPI) ranks \n", thread.n_pipeline, _world_size);
  }
}

// This operates in reverse order from boot_services

void
halt_services( void )
{
  _boot_timestamp = 0;

#if defined(VPIC_USE_PTHREADS)
  #if defined(VPIC_SWAP_MPI_PTHREAD_INIT)
    thread.halt();
    serial.halt();
    halt_mp();
  #else
    halt_mp();
    thread.halt();
    serial.halt();
  #endif

#elif defined(VPIC_USE_OPENMP)
  halt_mp();

#endif

  halt_checkpt();
}

