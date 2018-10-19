#define IN_sf_interface

#include "sf_interface_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper reduce_accumulator_array
// function.
//----------------------------------------------------------------------------//

void
reduce_accumulator_array( accumulator_array_t * RESTRICT aa )
{
  if ( !aa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  reduce_accumulator_array_pipeline( aa );
}
