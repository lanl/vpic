#ifndef _langevin_h_
#define _langevin_h_

#include "collision_private.h"

typedef struct langevin
{
  species_t  * sp;
  rng_pool_t * rp;
  float kT;
  float nu;
  int interval;
} langevin_t;

void
apply_langevin_pipeline( langevin_t * l );

void
langevin_pipeline_scalar( langevin_pipeline_args_t * RESTRICT args,
                          int pipeline_rank,
                          int n_pipeline );

#endif /* _langevin_h_ */
