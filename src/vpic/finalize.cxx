/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "vpic.hxx"

void vpic_simulation::finalize( void ) {
  // block for all message passing processes to finish
  mp_barrier( grid->mp );
  update_profile( mp_rank( grid->mp )==0 );
}
