#define IN_sfa

#include "sfa_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper compute_rms_div_b_err
// function.
//----------------------------------------------------------------------------//

double
compute_rms_div_b_err( const field_array_t * fa )
{
  double rms_div_b_err;

  if ( !fa )
  {
    ERROR( ( "Bad args") );
  }

  // Conditionally execute this when more abstractions are available.
  rms_div_b_err = compute_rms_div_b_err_pipeline( fa );

  return rms_div_b_err;
}
