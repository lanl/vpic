/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "sf_interface.h"

static int
aa_n_pipeline(void) {
  int                       n = serial.n_pipeline;
  if( n<thread.n_pipeline ) n = thread.n_pipeline;
# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  if( n<spu.n_pipeline    ) n = spu.n_pipeline;
# endif
  return n; /* max( {serial,thread,spu}.n_pipeline ) */
}

void
checkpt_accumulator_array( const accumulator_array_t * aa ) {
  CHECKPT( aa, 1 );
  CHECKPT_ALIGNED( aa->a, (size_t)(aa->n_pipeline+1)*(size_t)aa->stride, 128 );
  CHECKPT_PTR( aa->g );
}

accumulator_array_t *
restore_accumulator_array( void ) {
  accumulator_array_t * aa;
  RESTORE( aa );
  RESTORE_ALIGNED( aa->a );
  RESTORE_PTR( aa->g );
  if( aa->n_pipeline!=aa_n_pipeline() )
    ERROR(( "Number of accumulators restored is not the same as the number of "
            "accumulators checkpointed.  Did you change the number of threads "
            "per process between checkpt and restore?" ));
  return aa;
}

accumulator_array_t *
new_accumulator_array( grid_t * g ) {
  accumulator_array_t * aa;
  if( !g ) ERROR(( "Bad grid."));
  MALLOC( aa, 1 );
  aa->n_pipeline = aa_n_pipeline();
  aa->stride     = POW2_CEIL(g->nv,2);
  aa->g          = g;
  MALLOC_ALIGNED( aa->a, (size_t)(aa->n_pipeline+1)*(size_t)aa->stride, 128 );
  CLEAR( aa->a, (size_t)(aa->n_pipeline+1)*(size_t)aa->stride );
  REGISTER_OBJECT( aa, checkpt_accumulator_array, restore_accumulator_array,
                  NULL );
  return aa;
}

void
delete_accumulator_array( accumulator_array_t * aa ) {
  if( !aa ) return;
  UNREGISTER_OBJECT( aa );
  FREE_ALIGNED( aa->a );
  FREE( aa );
}

