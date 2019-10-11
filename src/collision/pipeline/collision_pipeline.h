#ifndef _collision_pipeline_h_
#define _collision_pipeline_h_

#include "../binary.h"
#include "../langevin.h"
#include "../unary.h"
#include "../takizuka_abe.h"

void
binary_pipeline_scalar( binary_collision_model_t * RESTRICT cm,
                        int pipeline_rank,
                        int n_pipeline );

void
langevin_pipeline_scalar( langevin_pipeline_args_t * RESTRICT args,
                          int pipeline_rank,
                          int n_pipeline );

void
unary_pipeline_scalar( unary_collision_model_t * RESTRICT cm,
                       int pipeline_rank,
                       int n_pipeline );

void
takizuka_abe_pipeline_scalar( takizuka_abe_t * RESTRICT cm,
                              int pipeline_rank,
                              int n_pipeline );

#endif /* _collision_pipeline_h_ */
