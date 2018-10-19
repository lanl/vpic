#ifndef _collision_private_h_
#define _collision_private_h_

#ifndef IN_collision
#error "Do not include collision_private.h; use collision.h"
#endif

#include "collision.h"

typedef void
(*collision_op_func_t)( void * params );

typedef void
(*delete_collision_op_func_t) ( struct collision_op * cop );

struct collision_op {
  void * params;
  collision_op_func_t apply;
  delete_collision_op_func_t delete_cop;
  collision_op_t * next;
};

BEGIN_C_DECLS

void
checkpt_collision_op_internal( const collision_op_t * cop );

collision_op_t *
restore_collision_op_internal( void * params );

collision_op_t *
new_collision_op_internal( void * params,
                           collision_op_func_t apply,
                           delete_collision_op_func_t delete_cop,
                           checkpt_func_t checkpt,
                           restore_func_t restore,
                           reanimate_func_t reanimate );

void
delete_collision_op_internal( collision_op_t * cop );

END_C_DECLS

///////////////////////////////////////////////////////////////////////////////
// Langevin pipeline interface

typedef struct langevin_pipeline_args {
  MEM_PTR( particle_t, 128 ) p;
  MEM_PTR( rng_t,      128 ) rng[ MAX_PIPELINE ];
  float decay; 
  float drive;
  int np;
  PAD_STRUCT( (1+MAX_PIPELINE)*SIZEOF_MEM_PTR+2*sizeof(float)+sizeof(int) )
} langevin_pipeline_args_t;

// PROTOTYPE_PIPELINE( langevin, langevin_pipeline_args_t );

void
apply_langevin_pipeline( langevin_t * l );

void
langevin_pipeline_scalar( langevin_pipeline_args_t * RESTRICT args,
                          int pipeline_rank,
                          int n_pipeline );

#endif /* _collision_h_ */
