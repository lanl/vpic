#define IN_spa

#include "spa_private.h"

//----------------------------------------------------------------------------//
// Include various programming model implementation files. For now, this is
// just the pipeline model. When there are more models, probably want these
// to be conditionally included.
//----------------------------------------------------------------------------//

#include "advance_p_pipeline.cc"

//----------------------------------------------------------------------------//
// Top level function to select and call particle advance function using the
// desired particle advance abstraction.  Currently, the only abstraction
// available is the pipeline abstraction.
//----------------------------------------------------------------------------//

void
advance_p( /**/  species_t            * RESTRICT sp,
           /**/  accumulator_array_t  * RESTRICT aa,
           const interpolator_array_t * RESTRICT ia )
{
  // Once more options are available, this should be conditionally executed
  // based on user choice.
  advance_p_pipeline( sp, aa, ia );
}
