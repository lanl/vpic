#define IN_sfa

#include "sfa_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper clean_div_e function.
//----------------------------------------------------------------------------//

void
clean_div_e( field_array_t * fa )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  clean_div_e_pipeline( fa );
}

