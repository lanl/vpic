// FIXME: USE THE DISCRETIZED VARIATIONAL DEFINITION OF ENERGY

#define IN_sfa

#include "sfa_private.h"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper vacuum_energy_f function.
//----------------------------------------------------------------------------//

void
vacuum_energy_f( double * global,
                 const field_array_t * RESTRICT fa )
{
  if ( !global || !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  vacuum_energy_f_pipeline( global, fa );
}

