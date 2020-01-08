#ifndef _clean_div_b_pipeline_h_
#define _clean_div_b_pipeline_h_

#ifndef IN_clean_div_b_pipeline
#error                                                                         \
    "Only include clean_div_b_pipeline.h in clean_div_b_pipeline source files."
#endif

#include "../../field_advance.h"

typedef struct pipeline_args
{
    field_t *ALIGNED( 128 ) f;
    const grid_t *g;
} pipeline_args_t;

#define f( x, y, z ) f[VOXEL( x, y, z, nx, ny, nz )]

#define MARDER_CBX() f0->cbx += px * ( f0->div_b_err - fx->div_b_err )
#define MARDER_CBY() f0->cby += py * ( f0->div_b_err - fy->div_b_err )
#define MARDER_CBZ() f0->cbz += pz * ( f0->div_b_err - fz->div_b_err )

void clean_div_b_pipeline_scalar( pipeline_args_t *args, int pipeline_rank,
                                  int n_pipeline );

void clean_div_b_pipeline_v4( pipeline_args_t *args, int pipeline_rank,
                              int n_pipeline );

void clean_div_b_pipeline_v8( pipeline_args_t *args, int pipeline_rank,
                              int n_pipeline );

void clean_div_b_pipeline_v16( pipeline_args_t *args, int pipeline_rank,
                               int n_pipeline );

#endif // _clean_div_b_pipeline_h_
