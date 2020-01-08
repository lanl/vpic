#ifndef _sf_interface_private_h_
#define _sf_interface_private_h_

#ifndef IN_sf_interface
#error "Do not include sf_interface_private.h; include sf_interface.h"
#endif

#include "sf_interface.h"

///////////////////////////////////////////////////////////////////////////////
// load_interpolator_pipeline interface

void load_interpolator_array_pipeline( interpolator_array_t *RESTRICT ia,
                                       const field_array_t *RESTRICT fa );

///////////////////////////////////////////////////////////////////////////////
// clear_accumulators_pipeline interface

void clear_accumulator_array_pipeline( accumulator_array_t *RESTRICT aa );

void reduce_accumulator_array_pipeline( accumulator_array_t *RESTRICT aa );

///////////////////////////////////////////////////////////////////////////////

void unload_accumulator_array_pipeline(
    field_array_t *RESTRICT fa, const accumulator_array_t *RESTRICT aa );

#endif // _sf_interface_private_h_
