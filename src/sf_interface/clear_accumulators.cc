#define IN_sf_interface

#include "sf_interface_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper clear_accumulator_array
// function.
//----------------------------------------------------------------------------//

void
clear_accumulator_array( accumulator_array_t * RESTRICT aa )
{
  if ( !aa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  clear_accumulator_array_pipeline( aa );
}
