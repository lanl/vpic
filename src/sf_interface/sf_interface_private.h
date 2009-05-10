#ifndef _sf_interface_private_h_
#define _sf_interface_private_h_

#ifndef IN_sf_interface
#error "Do not include sf_interface_private.h; include sf_interface.h"
#endif

#include "sf_interface.h"

// FIXME: THIS SHOULD BE ABSTRACTED SOMEWHERE ELSE

#define FOR_SPU ( defined(CELL_SPU_BUILD)        || \
                  ( defined(CELL_PPU_BUILD)    &&   \
                     defined(USE_CELL_SPUS)    &&   \
                     defined(HAS_SPU_PIPELINE) ) )

#if FOR_SPU

# if defined(CELL_PPU_BUILD) 

    // Use SPU dispatcher on the SPU pipeline
    // PPU will do straggler cleanup with scalar pipeline

#   define EXEC_PIPELINES(name,args,sz_args)                      \
    spu.dispatch( (pipeline_func_t)( (size_t)                     \
        ( root_segment_##name##_pipeline_spu) ), args, sz_args ); \
    name##_pipeline( args, spu.n_pipeline, spu.n_pipeline )

#   define WAIT_PIPELINES() spu.wait()

#   define N_PIPELINE spu.n_pipeline

#   define PROTOTYPE_PIPELINE( name, args_t )                      \
    extern uint32_t root_segment_##name##_pipeline_spu;            \
                                                                   \
    void                                                           \
    name##_pipeline( args_t * args,                                \
                     int pipeline_rank,                            \
                     int n_pipeline )

#   define PAD_STRUCT( sz ) char _pad[ PAD( (sz), 16 ) ];

# else

    // SPUs cannot dispatch pipelines

#   define PROTOTYPE_PIPELINE( name, args_t )                   \
    void                                                        \
    _SPUEAR_##name##_pipeline_spu( MEM_PTR( args_t, 128 ) argp, \
                                   int pipeline_rank,           \
                                   int n_pipeline )

#   define PAD_STRUCT( sz ) char _pad[ PAD( (sz), 16 ) ];

# endif

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  // Use thread dispatcher on the v4 pipeline
  // Caller will do straggler cleanup with scalar pipeline

# define EXEC_PIPELINES(name,args,sz_args)                               \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v4, args, sz_args ); \
  name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

# define WAIT_PIPELINES() thread.wait()

# define N_PIPELINE thread.n_pipeline

# define PROTOTYPE_PIPELINE( name, args_t ) \
  void                                      \
  name##_pipeline_v4( args_t * args,        \
                      int pipeline_rank,    \
                      int n_pipeline );     \
                                            \
  void                                      \
  name##_pipeline( args_t * args,           \
                   int pipeline_rank,       \
                   int n_pipeline )

# define PAD_STRUCT( sz )

#else

  // Use thread dispatcher on the scalar pipeline
  // Caller will do straggler cleanup with scalar pipeline

# define EXEC_PIPELINES(name,args,sz_args)                              \
  thread.dispatch( (pipeline_func_t)name##_pipeline, args, sz_args );   \
  name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

# define WAIT_PIPELINES() thread.wait()

# define N_PIPELINE thread.n_pipeline

# define PROTOTYPE_PIPELINE( name, args_t ) \
  void                                      \
  name##_pipeline( args_t * args,           \
                   int pipeline_rank,       \
                   int n_pipeline )

# define PAD_STRUCT( sz )

#endif

///////////////////////////////////////////////////////////////////////////////
// load_interpolator_pipeline interface

typedef struct load_interpolator_pipeline_args {

  interpolator_t * ALIGNED(128) fi;
  const field_t  * ALIGNED(128) f;
  const grid_t   *              g;

  PAD_STRUCT( 3*SIZEOF_MEM_PTR )

} load_interpolator_pipeline_args_t;

PROTOTYPE_PIPELINE( load_interpolator, load_interpolator_pipeline_args_t );

///////////////////////////////////////////////////////////////////////////////
// clear_accumulators_pipeline interface

// Pipelines are be assigned accumulator blocks in multiples of 256
// (16KB) which is particularly convenient on Cell.  The pipeline
// dispatcher will handle any stragglers.

enum { accumulators_n_block = 256 };

typedef struct accumulators_pipeline_args {

  MEM_PTR( accumulator_t, 128) a; // First accumulator to reduce
  int n;                          // Number of accumulators to reduce
  int n_array;                    // Number of accumulator arrays
  int s_array;                    // Stride between each array

  PAD_STRUCT( SIZEOF_MEM_PTR + 3*sizeof(int) )

} accumulators_pipeline_args_t;

PROTOTYPE_PIPELINE( clear_accumulators,  accumulators_pipeline_args_t );
PROTOTYPE_PIPELINE( reduce_accumulators, accumulators_pipeline_args_t );

///////////////////////////////////////////////////////////////////////////////

typedef struct unload_accumulator_pipeline_args {

  field_t             * ALIGNED(128) f;
  const accumulator_t * ALIGNED(128) a;
  const grid_t        *              g;

  PAD_STRUCT( 3*SIZEOF_MEM_PTR )

} unload_accumulator_pipeline_args_t;

PROTOTYPE_PIPELINE( unload_accumulator, unload_accumulator_pipeline_args_t );

#undef FOR_SPU

#endif // _sf_interface_private_h_
