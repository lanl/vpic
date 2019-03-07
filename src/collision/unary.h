#ifndef _unary_h_
#define _unary_h_

#include "collision_private.h"

typedef struct unary_collision_model
{
  char * name;
  unary_rate_constant_func_t rate_constant;
  unary_collision_func_t collision;
  void * params;
  species_t * sp;
  rng_pool_t * rp;
  int interval;
  int n_large_pr[ MAX_PIPELINE ];
} unary_collision_model_t;

void
apply_unary_collision_model_pipeline( unary_collision_model_t * cm );

#endif /* _unary_h_ */
