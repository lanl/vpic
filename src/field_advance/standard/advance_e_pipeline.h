#ifndef _advance_e_pipeline_h_
#define _advance_e_pipeline_h_

#include "sfa_private.h"

typedef struct pipeline_args
{
  field_t            * ALIGNED(128) f;
  const sfa_params_t *              p;
  const grid_t       *              g;
} pipeline_args_t;

void
advance_e_pipeline_scalar( pipeline_args_t * args,
                           int pipeline_rank,
                           int n_pipeline );

void
advance_e_pipeline_v4( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline );

void
advance_e_pipeline_v8( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline );

void
advance_e_pipeline_v16( pipeline_args_t * args,
                        int pipeline_rank,
                        int n_pipeline );

#endif // _advance_e_pipeline_h_
