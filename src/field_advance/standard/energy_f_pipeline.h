#ifndef _energy_f_pipeline_h_
#define _energy_f_pipeline_h_

#include "sfa_private.h"

typedef struct pipeline_args
{
  const field_t      * ALIGNED(128) f;
  const sfa_params_t *              p;
  const grid_t       *              g;
  double en[ MAX_PIPELINE+1 ][ 6 ];
} pipeline_args_t;

void
energy_f_pipeline_scalar( pipeline_args_t * args,
                          int pipeline_rank,
                          int n_pipeline );

#endif // _energy_f_pipeline_h_
