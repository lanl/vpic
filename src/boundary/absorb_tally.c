#include "boundary.h"

// Special absorbing boundary condition that tallies the number of
// particles absorbed of each species.
//
// Written by:  Brian J. Albright, X-1, LANL   July, 2005

int
absorb_tally( void * params,
              particle_t *r, 
              particle_mover_t *pm,
              field_t *f,
              accumulator_t *a, 
              const grid_t *g,
              species_t *s, 
              particle_injector_t *pi,
              mt_rng_t *rng,
              int face ) {
  absorb_tally_t * absorb_tally_data = (absorb_tally_t *)params;
  int ispec; 

  // obtain reflux data for particle according to particle species
  for ( ispec=0; 
        ispec<absorb_tally_data->nspec && s->id!=absorb_tally_data->id[ispec]; 
	++ispec )
    ; 
  if ( ispec==absorb_tally_data->nspec ) 
    ERROR(("Unknown species passed to boundary handler."));
  ++(absorb_tally_data->nabs[ispec]);

  accumulate_rhob( f, r, g );
  return 0;
}
