#ifndef _binary_h_
#define _binary_h_

#include "collision_private.h"

typedef struct binary_collision_model
{
  char * name;
  binary_rate_constant_func_t rate_constant;
  binary_collision_func_t collision;
  void * params;
  species_t  * spi;
  species_t  * spj;
  rng_pool_t * rp;
  double sample;
  int interval;
  int n_large_pr[ MAX_PIPELINE ];
  int strategy;
} binary_collision_model_t;

void
apply_binary_collision_model_pipeline( binary_collision_model_t * cm );

#endif /* _binary_h_ */
