#define IN_sfa

#include "sfa_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper vacuum_clean_div_e
// function.
//----------------------------------------------------------------------------//

void
vacuum_clean_div_e( field_array_t * fa )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  vacuum_clean_div_e_pipeline( fa );
}
