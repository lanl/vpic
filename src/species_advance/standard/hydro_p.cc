#define IN_spa

#include "../species_advance.h"

//----------------------------------------------------------------------------//
// Top level function to select and call particle hydro accumulation function
// using the desired particle advance abstraction.  Currently, the only
// abstraction available is the pipeline abstraction.
//----------------------------------------------------------------------------//

void
accumulate_hydro_p( hydro_array_t * RESTRICT ha,
                    const species_t * RESTRICT sp,
                    const interpolator_array_t * RESTRICT ia )
{
  // Once more options are available, this should be conditionally executed
  // based on user choice.
  accumulate_hydro_p_pipeline( ha, sp, ia );
}
