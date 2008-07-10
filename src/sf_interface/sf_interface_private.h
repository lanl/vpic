#ifndef _sf_interface_private_h_
#define _sf_interface_private_h_

#ifndef IN_sf_interface
#error "Do not include sf_interface_private.h; include sf_interface.h"
#endif

#include "sf_interface.h"

// FIXME: THIS IS COMMON ACROSS MULTIPLE MODULES NOW AND SHOULD BE
// ABSTRACTED SOMEWHERE

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
// Field-particle coupling interfaces

typedef struct load_interpolator_pipeline_args {
  interpolator_t * ALIGNED(128) fi;
  const field_t  * ALIGNED(128) f;
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
  field_t             * ALIGNED(128) f;
  const accumulator_t * ALIGNED(128) a;
  const grid_t        *              g;
} unload_accumulator_pipeline_args_t;

#undef FOR_SPU

#endif // _sf_interface_private_h_
