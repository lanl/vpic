#define IN_sfa

#include "sfa_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper advance_b function.
//----------------------------------------------------------------------------//

void
advance_b( field_array_t * RESTRICT fa,
           float _frac )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  advance_b_pipeline( fa, _frac );
}
