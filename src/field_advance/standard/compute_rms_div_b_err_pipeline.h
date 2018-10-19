#ifndef _compute_rms_div_b_err_pipeline_h_
#define _compute_rms_div_b_err_pipeline_h_

#include "../field_advance.h"

typedef struct pipeline_args
{
  const field_t * ALIGNED(128) f;
  const grid_t  *              g;
  double err[MAX_PIPELINE+1];
} pipeline_args_t;

static void
compute_rms_div_b_err_pipeline_scalar( pipeline_args_t * args,
                                       int pipeline_rank,
                                       int n_pipeline );

#endif // _compute_rms_div_b_err_pipeline_h_
