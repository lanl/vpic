#define IN_spa

#include "spa_private.h"

//----------------------------------------------------------------------------//
// Include various programming model implementation files. For now, this is
// just the pipeline model. When there are more models, probably want these
// to be conditionally included.
//----------------------------------------------------------------------------//

#include "center_p_pipeline.cc"

//----------------------------------------------------------------------------//
// Top level function to select and call particle center function using the
// desired particle center abstraction.  Currently, the only abstraction
// available is the pipeline abstraction.
//----------------------------------------------------------------------------//

void
center_p( /**/  species_t            * RESTRICT sp,
          const interpolator_array_t * RESTRICT ia )
{
  // Once more options are available, this should be conditionally executed
  // based on user choice.
  center_p_pipeline( sp, ia );
}
