#ifndef _compute_rms_div_e_err_pipeline_h_
#define _compute_rms_div_e_err_pipeline_h_

#ifndef IN_compute_rms_div_e_err_pipeline
#error "Only include compute_rms_div_e_err_pipeline.h in compute_rms_div_e_err_pipeline source files."
#endif

#include "../field_advance.h"

typedef struct pipeline_args
{
  const field_t * ALIGNED(128) f;
  const grid_t  *              g;
  double err[MAX_PIPELINE+1];
} pipeline_args_t;

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

static void
compute_rms_div_e_err_pipeline_scalar( pipeline_args_t * args,
                                       int pipeline_rank,
                                       int n_pipeline );

#endif // _compute_rms_div_e_err_pipeline_h_
