/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <vpic.hxx>

void vpic_simulation::finalize( void ) {

  // FIXME!!!
  // stop pipelines
  // This is a hack to get us through the RR assessment.  At some
  // point this will have to be re-worked to use overlays to allow
  // for multple accelerated implementations
  advance_p_finalize();

  // block for all message passing processes to finish
  mp_barrier( grid->mp );
}
