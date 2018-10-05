#define IN_sfa

#include "sfa_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper vacuum_compute_curl_b
// function.
//----------------------------------------------------------------------------//

void
vacuum_compute_curl_b( field_array_t * RESTRICT fa )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  vacuum_compute_curl_b_pipeline( fa );
}
