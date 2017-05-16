/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised from earlier V4PIC versions
 *
 */

#include "vpic/vpic.h"

#ifdef USE_CATALYST
#include "VPICAdaptor.h"
#endif

/* The simulation variable is set up this way so both the checkpt
   service and main can see it.  This allows main to find where
   the restored objects are after a restore. */

vpic_simulation * simulation = NULL;

void
checkpt_main( vpic_simulation ** _simulation ) {
  CHECKPT_PTR( simulation );
}

vpic_simulation **
restore_main( void ) {
  RESTORE_PTR( simulation );
  return &simulation;
}

void
checkpt( const char * fbase,
         int tag ) {
  char fname[256];
  if( !fbase ) ERROR(( "NULL filename base" ));
  sprintf( fname, "%s.%i.%i", fbase, tag, world_rank );
  if( world_rank==0 ) log_printf( "*** Checkpointing to \"%s\"\n", fbase );
  checkpt_objects( fname );
}

int
main( int argc,
      char **argv ) {
  boot_services( &argc, &argv );
 
  const char * fbase = strip_cmdline_string(&argc, &argv, "--restore", NULL);
  if( fbase ) {

    // We are restoring from a checkpoint.  Determine checkpt file
    // for this process, restore all the objects in that file,
    // wait for all other processes to finishing restoring (such
    // that communication within reanimate functions is safe),
    // reanimate all the objects and issue a final barrier to
    // so that all processes come of a restore together.

    if( world_rank==0 ) log_printf( "*** Restoring from \"%s\"\n", fbase );
    char fname[256];
    sprintf( fname, "%s.%i", fbase, world_rank );
    restore_objects( fname );
    mp_barrier();
    reanimate_objects();
    mp_barrier();

  } else {

    // We are initializing from scratch.

    if( world_rank==0 ) log_printf( "*** Initializing\n" );
    simulation = new vpic_simulation;
    simulation->initialize( argc, argv );
    REGISTER_OBJECT( &simulation, checkpt_main, restore_main, NULL );

  }
 
  // Do any post init/restore simulation modifications
  // FIXME-KJB: STRIP_CMDLINE COULD MAKE THIS CLEANER AND MORE POWERFUL.
 
  fbase = strip_cmdline_string( &argc, &argv, "--modify", NULL );
  if( fbase ) {
    if( world_rank==0 ) log_printf( "*** Modifying from \"%s\"\n", fbase );
    simulation->modify( fbase );  
  }
 
  // Advance the simulation

  if( world_rank==0 ) log_printf( "*** Advancing\n" );
  double elapsed = wallclock();
  while( simulation->advance() ); 
  elapsed = wallclock() - elapsed;
  if( world_rank==0 ) {
    int  s = (int)elapsed, m  = s/60, h  = m/60, d  = h/24, w = d/ 7;
    /**/ s -= m*60,        m -= h*60, h -= d*24, d -= w*7;
    log_printf( "*** Done (%gs / %iw:%id:%ih:%im:%is elapsed)\n",
                elapsed, w, d, h, m, s );
  }

  // Cleaning up
#ifdef USE_CATALYST
  coprocessorfinalize();
#endif

  if( world_rank==0 ) log_printf( "*** Cleaning up\n" );
  UNREGISTER_OBJECT( &simulation );
  simulation->finalize();
  delete simulation;
  if( world_rank==0 ) log_printf( "normal exit\n" ); 

  halt_services();
  return 0;
}

