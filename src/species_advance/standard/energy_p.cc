#define IN_spa

#include "spa_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call particle energy function using the
// desired particle energy abstraction.  Currently, the only abstraction
// available is the pipeline abstraction.
//----------------------------------------------------------------------------//

double
energy_p( const species_t * RESTRICT sp,
          const interpolator_array_t * RESTRICT ia )
{
  double energy_particles;

  // Once more options are available, this should be conditionally executed
  // based on user choice.
  energy_particles = energy_p_pipeline( sp, ia );

  return energy_particles;
}
