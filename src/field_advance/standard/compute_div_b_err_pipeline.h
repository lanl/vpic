#ifndef _compute_div_b_err_pipeline_h_
#define _compute_div_b_err_pipeline_h_

#include "../field_advance.h"

typedef struct pipeline_args
{
  field_t      * ALIGNED(128) f;
  const grid_t *              g;
} pipeline_args_t;

void
compute_div_b_err_pipeline_scalar( pipeline_args_t * args,
                                   int pipeline_rank,
                                   int n_pipeline );

#endif // _compute_div_b_err_pipeline_h_
