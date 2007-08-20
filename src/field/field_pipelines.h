#ifndef _field_pipelines_h_
#define _field_pipelines_h_

#ifndef IN_field_pipeline
#error "Do not include field_pipelines.h; include field.h"
#endif

#include <field.h>

#define FOR_SPU ( defined(CELL_SPU_BUILD)       || \
                  ( defined(CELL_PPU_BUILD)   &&   \
                    defined(USE_CELL_SPUS)    &&   \
                    defined(HAS_SPU_PIPELINE) ) )

#if FOR_SPU

# if defined(CELL_PPU_BUILD)

    // Use the SPU dispatcher on the SPU pipeline

#   define EXEC_PIPELINES(name,args,sz_args)                          \
    spu.dispatch( SPU_PIPELINE(name##_pipeline_spu), args, sz_args ); \
    name##_pipeline( args, spu.n_pipeline, spu.n_pipeline )

#   define WAIT_PIPELINES() spu.wait()

#   define N_PIPELINE       spu.n_pipeline

# else

    // SPUs cannot dispatch pipelines

# endif

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  // Use the thread dispatcher on the v4 pipeline

# define EXEC_PIPELINES(name,args,sz_args)                               \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v4, args, sz_args ); \
  name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

# define WAIT_PIPELINES() thread.wait()

# define N_PIPELINE thread.n_pipeline

#else

  // Use the thread dispatcher on the scalar pipeline

# define EXEC_PIPELINES(name,args,sz_args)                            \
  thread.dispatch( (pipeline_func_t)name##_pipeline, args, sz_args ); \
  name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

# define WAIT_PIPELINES() thread.wait()

# define N_PIPELINE       thread.n_pipeline

#endif

///////////////////////////////////////////////////////////////////////////////
// Field advance pipelines interfaces

typedef struct advance_b_pipeline_args {
  field_t      * ALIGNED(16) f;
  const grid_t *             g;
  float frac;
} advance_b_pipeline_args_t;

typedef struct advance_e_pipeline_args {
  field_t                      * ALIGNED(16) f;
  const material_coefficient_t * ALIGNED(16) m;
  const grid_t                 *             g;
} advance_e_pipeline_args_t;

///////////////////////////////////////////////////////////////////////////////
// Field-particle coupling interfaces

typedef struct load_interpolator_pipeline_args {
  interpolator_t * ALIGNED(128) fi;
  const field_t  * ALIGNED(16)  f;
  const grid_t   *              g;
} load_interpolator_pipeline_args_t;

typedef struct clear_accumulators_pipeline_args {
  accumulator_t * ALIGNED(128) a;       // Base of all the accumulators
  int                          n_voxel; // Total number of voxels
} clear_accumulators_pipeline_args_t;

typedef struct reduce_accumulators_pipeline_args {
  accumulator_t * ALIGNED(128) a;
  const grid_t  *              g;
  int na;                               // Number of accumulators
} reduce_accumulators_pipeline_args_t;

typedef struct unload_accumulator_pipeline_args {
  field_t             * ALIGNED(16)  f;
  const accumulator_t * ALIGNED(128) a;
  const grid_t        *              g;
} unload_accumulator_pipeline_args_t;

//////////////////////////////////////////////////////////////////////////////
// Field diagnostics and initialization interfaces

typedef struct energy_f_pipeline_args {
  const field_t                * ALIGNED(16) f;
  const material_coefficient_t * ALIGNED(16) m;
  const grid_t                 *             g;
  double en[MAX_PIPELINE+1][6];
} energy_f_pipeline_args_t;

typedef struct compute_curl_b_pipeline_args {
  field_t                      * ALIGNED(16) f;
  const material_coefficient_t * ALIGNED(16) m;
  const grid_t                 *             g;
} compute_curl_b_pipeline_args_t;

///////////////////////////////////////////////////////////////////////////////
// Electric field divergence cleaning interfaces

typedef struct compute_rhob_pipeline_args {
  field_t                      * ALIGNED(16) f;
  const material_coefficient_t * ALIGNED(16) m;
  const grid_t                 *             g;
} compute_rhob_pipeline_args_t;

typedef struct compute_div_e_err_pipeline_args {
  field_t                      * ALIGNED(16) f;
  const material_coefficient_t * ALIGNED(16) m;
  const grid_t                 *             g;
} compute_div_e_err_pipeline_args_t;

typedef struct compute_rms_div_e_err_pipeline_args {
  field_t      * ALIGNED(16) f;
  const grid_t *             g;
  double err[MAX_PIPELINE+1];
} compute_rms_div_e_err_pipeline_args_t;

typedef struct clean_div_e_pipeline_args {
  field_t                      * ALIGNED(16) f;
  const material_coefficient_t * ALIGNED(16) m;
  const grid_t                 *             g;
} clean_div_e_pipeline_args_t;

///////////////////////////////////////////////////////////////////////////////
// Magnetic field divergence cleaning interfaces

typedef struct compute_div_b_err_pipeline_args {
  field_t      * ALIGNED(16) f;
  const grid_t *             g;
} compute_div_b_err_pipeline_args_t;

typedef struct compute_rms_div_b_err_pipeline_args {
  field_t      * ALIGNED(16) f;
  const grid_t *             g;
  double err[MAX_PIPELINE+1];
} compute_rms_div_b_err_pipeline_args_t;

typedef struct clean_div_b_pipeline_args {
  field_t      * ALIGNED(16) f;
  const grid_t *             g;
} clean_div_b_pipeline_args_t;

#undef FOR_SPU

#endif // _field_pipelines_h_
