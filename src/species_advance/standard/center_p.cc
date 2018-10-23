#define IN_spa

#include "../species_advance.h"

//----------------------------------------------------------------------------//
// Top level function to select and call particle center function using the
// desired particle center abstraction.  Currently, the only abstraction
// available is the pipeline abstraction.
//----------------------------------------------------------------------------//

void
center_p( species_t * RESTRICT sp,
          const interpolator_array_t * RESTRICT ia )
{
  // Once more options are available, this should be conditionally executed
  // based on user choice.
  center_p_pipeline( sp, ia );
}
