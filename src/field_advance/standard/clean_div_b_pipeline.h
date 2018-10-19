#ifndef _clean_div_b_pipeline_h_
#define _clean_div_b_pipeline_h_

#include "../field_advance.h"

typedef struct pipeline_args
{
  field_t      * ALIGNED(128) f;
  const grid_t *              g;
} pipeline_args_t;

void
clean_div_b_pipeline_scalar( pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline );

void
clean_div_b_pipeline_v4( pipeline_args_t * args,
                         int pipeline_rank,
                         int n_pipeline );

void
clean_div_b_pipeline_v8( pipeline_args_t * args,
                         int pipeline_rank,
                         int n_pipeline );

void
clean_div_b_pipeline_v16( pipeline_args_t * args,
                          int pipeline_rank,
                          int n_pipeline );

#endif // _clean_div_b_pipeline_h_
