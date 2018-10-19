#ifndef _vacuum_clean_div_e_pipeline_h_
#define _vacuum_clean_div_e_pipeline_h_

#include "sfa_private.h"

typedef struct pipeline_args
{
  field_t            * ALIGNED(128) f;
  const sfa_params_t *              p;
  const grid_t       *              g;
} pipeline_args_t;

static void
vacuum_clean_div_e_pipeline_scalar( pipeline_args_t * args,
                                    int pipeline_rank,
                                    int n_pipeline );

#endif // _vacuum_clean_div_e_pipeline_h_
