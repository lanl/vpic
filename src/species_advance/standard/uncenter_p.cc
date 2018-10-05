#define IN_spa

#include "spa_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call particle uncenter function using the
// desired particle center abstraction.  Currently, the only abstraction
// available is the pipeline abstraction.
//----------------------------------------------------------------------------//

void
uncenter_p( /**/  species_t            * RESTRICT sp,
            const interpolator_array_t * RESTRICT ia )
{
  // Once more options are available, this should be conditionally executed
  // based on user choice.
  uncenter_p_pipeline( sp, ia );
}
