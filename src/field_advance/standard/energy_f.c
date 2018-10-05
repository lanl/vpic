// FIXME: USE THE DISCRETIZED VARIATIONAL PRINCIPLE DEFINITION OF ENERGY

#define IN_sfa

#include "sfa_private.h"

//----------------------------------------------------------------------------//
// Include various programming model implementation files. For now, this is
// just the pipeline model. When there are more models, probably want these
// to be conditionally included.
//----------------------------------------------------------------------------//

#include "energy_f_pipeline.c"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper energy_f function.
//----------------------------------------------------------------------------//

void
energy_f( double * global,
          const field_array_t * RESTRICT fa )
{
  if ( !global || !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  energy_f_pipeline( global, fa );
}

