#ifndef _chemical_h_
#define _chemical_h_

#include "collision_private.h"

typedef struct chemical_collision_model {
  char * name;
  chemical_rate_constant_func_t rate_constant;
  chemical_collision_func_t collision;
  void * params;
  rng_pool_t * rp;
  field_array_t * fa;
  int interval;

  // Species information.
  int n_reactants, n_products, *consumable;
  species_t  **reactants, **products;

  // Sampling parameters.
  double sample;
  int dynamic_sampling;

  // Information about the number of collisions and their probabilities.
  int n_large_pr[ MAX_PIPELINE ];
  int n_tested[ MAX_PIPELINE ];
  int unaccumulated[ MAX_PIPELINE ];
  float pr_max[ MAX_PIPELINE ];

  // Information on the number of particles produced and modified.
  int *n_modified;                     // size = N * MAX_PIPELINE
  particle_mover_t **reactant_movers;  // size = N * MAX_PIPELINE
  int *n_produced;                     // size =     MAX_PIPELINE
  particle_t **product_particles;      // size = M * MAX_PIPELINE

} chemical_collision_model_t;

void
apply_chemical_collision_model_pipeline( chemical_collision_model_t * cm );

// Helper macros for looping.
#define FOR_PRODUCTS  for( sp=products[0],  m=0 ; m<M ; sp=products[++m]  )
#define FOR_REACTANTS for( sp=reactants[0], n=0 ; n<N ; sp=reactants[++n] )

#endif /* _chemical_h_ */
