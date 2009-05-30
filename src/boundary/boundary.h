#ifndef _boundary_h_
#define _boundary_h_

#include "../species_advance/standard/spa.h"

BEGIN_C_DECLS

// In maxwellian_reflux.c

// FIXME: MAXWELLIAN_REFLUX SILENTLY ASSUMES LESS THAN 32 SPECIES

typedef struct maxwellian_reflux {
  float ut_perp[32]; // array of perp thermal velocities 
  float ut_para[32]; // array of parallel thermal velocities
} maxwellian_reflux_t; 

int
maxwellian_reflux( void * params,
                   particle_t *r,
                   particle_mover_t *pm, 
                   field_t *f,
                   accumulator_t *a,
                   const grid_t *g,
                   species_t *s,
                   particle_injector_t *pi,
                   mt_rng_t *rng,
                   int face ); 

// In absorb_tally.c

typedef struct absorb_tally {
  int nspec;         // number of species interacting with boundary
  species_id id[32]; // array of species ids
  int nabs[32];      // array of numbers of particles absorbed
} absorb_tally_t; 

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
              int face ); 

// In absorb_tally.c

typedef struct link_boundary {
  char fbase[256]; // base of file name to contain link info
  double n_out;    // number of writes so far on this node (double to
                   // accomodate long long runs)
} link_boundary_t;

int
link_boundary( void * params,
               particle_t *r,
               particle_mover_t *pm, 
               field_t *f,
               accumulator_t *a,
               const grid_t *g,
               species_t *s,
               particle_injector_t *pi,
               mt_rng_t *rng,
               int face ); 

END_C_DECLS

#endif // _boundary_h_

