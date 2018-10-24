// Note: This is virtually identical to vacuum_compute_rhob

#define IN_sfa

#include "sfa_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper vacuum_compute_div_e_err
// function.
//----------------------------------------------------------------------------//

void
vacuum_compute_div_e_err( field_array_t * RESTRICT fa )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  vacuum_compute_div_e_err_pipeline( fa );
}
