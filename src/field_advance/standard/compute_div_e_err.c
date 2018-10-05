// Note: This is virtually identical to compute_rhob

#define IN_sfa

#include "sfa_private.h"

//----------------------------------------------------------------------------//
// Include various programming model implementation files. For now, this is
// just the pipeline model. When there are more models, probably want these
// to be conditionally included.
//----------------------------------------------------------------------------//

#include "compute_div_e_err_pipeline.c"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper compute_div_e_err function.
//----------------------------------------------------------------------------//

void
compute_div_e_err( field_array_t * RESTRICT fa )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  compute_div_e_err_pipeline( fa );
}
