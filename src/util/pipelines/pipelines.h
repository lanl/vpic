#ifndef _pipelines_h_
#define _pipelines_h_

#include <util_base.h>

enum { MAX_PIPELINE = 16 };

// A pipeline function takes a pointer to arguments for the pipeline
// and a integer which gives the rank of the pipeline and the total
// number of pipelines dispatched.

typedef void
(*pipeline_func_t)( void * args,
                    int pipeline_rank,
                    int n_pipeline );

///////////////////////////////////////////////////////////////////////////////

#define FOR_SPU ( defined(CELL_SPU_BUILD)        || \
                  ( defined(CELL_PPU_BUILD)    &&   \
                     defined(USE_CELL_SPUS)    &&   \
                     defined(HAS_SPU_PIPELINE) ) )

#if FOR_SPU

# if defined(CELL_PPU_BUILD) 

    // Use SPU dispatcher on the SPU pipeline
    // PPU will do straggler cleanup with scalar pipeline

#   define N_PIPELINE spu.n_pipeline
#   define EXEC_PIPELINES(name,args,str)                               \
    spu.dispatch( (pipeline_func_t)                                    \
                  ((size_t)(root_segment_##name##_pipeline_spu)),      \
                  args, sizeof(*args), str );                          \
    name##_pipeline( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )
#   define WAIT_PIPELINES() spu.wait()

#   define PROTOTYPE_PIPELINE( name, args_t )                          \
    extern uint32_t root_segment_##name##_pipeline_spu;                \
                                                                       \
    void                                                               \
    name##_pipeline( args_t * args,                                    \
                     int pipeline_rank,                                \
                     int n_pipeline )

#   define PAD_STRUCT( sz ) char _pad[ PAD( (sz), 16 ) ];

# else

    // SPUs cannot dispatch pipelines

#   define PROTOTYPE_PIPELINE( name, args_t )                          \
    void                                                               \
    _SPUEAR_##name##_pipeline_spu( args_t * args,                      \
                                   int pipeline_rank,                  \
                                   int n_pipeline )

#   define PAD_STRUCT( sz ) char _pad[ PAD( (sz), 16 ) ];
# endif

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  // Use thread dispatcher on the v4 pipeline
  // Caller will do straggler cleanup with scalar pipeline

# define N_PIPELINE thread.n_pipeline
# define EXEC_PIPELINES(name,args,str)                                 \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v4,                \
                   args, sizeof(*args), str );                         \
  name##_pipeline( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )
# define WAIT_PIPELINES() thread.wait()

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

# define N_PIPELINE thread.n_pipeline
# define EXEC_PIPELINES(name,args,str)                                  \
  thread.dispatch( (pipeline_func_t)name##_pipeline,                    \
                   args, sizeof(*args), str );                          \
  name##_pipeline( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )
# define WAIT_PIPELINES() thread.wait()

# define PROTOTYPE_PIPELINE( name, args_t ) \
  void                                      \
  name##_pipeline( args_t * args,           \
                   int pipeline_rank,       \
                   int n_pipeline )

# define PAD_STRUCT( sz )

#endif

///////////////////////////////////////////////////////////////////////////////

#if !defined(CELL_SPU_BUILD)

typedef struct pipeline_dispatcher {

  // n_pipelines indicates the number of pipelines currently running.
  // Technically, this should be read only for users!

  int n_pipeline;

  // boot creates the number of pipelines requested.  Generally, this
  // is number of cores on a node if using symmetric multiprocessing
  // or the number of pipeline processors if using heterogeneous
  // multiprocessing.  Dispatch to host is used to indicate whether
  // or not the host (e.g. the thread which call this) should be
  // used to process pipeline threads.  On some dispatchers, this
  // flag is meaningless (e.g. in the serial dispatcher, all pipelines
  // are executed on the host) or invalid to set (e.g. in the SPU
  // dispatcher ... pipelines physically cannot be executed on the
  // host).

  void (*boot)( int n_pipeline,
                int dispatch_to_host );

  // halt destroys all the resources used by the dispatcher created
  // in boot.

  void (*halt)( void );

  // dispatch begins executing the given pipeline function on all the
  // pipelines.
  //
  // pipeline is the pipeline function to execute on the pipelines.
  //
  // args is an array of arguments to pass to each pipeline.
  //
  // sz gives the byte size of an element of the argument
  // array.
  //
  // str gives the element stride between elements of the argument
  // array.  Pass 0 if you want all pipelines to get the same
  // arguments.
  //
  // If the pipeline functions do not take arguments, use NULL for
  // args and 0 for sz and str
                    
  void
  (*dispatch)( pipeline_func_t pipeline,
               void * args,
               int sz,
               int str );

  // wait waits for the previous dispatch to complete.

  void
  (*wait)( void );

} pipeline_dispatcher_t;

BEGIN_C_DECLS

extern pipeline_dispatcher_t serial; // For debugging purposes
extern pipeline_dispatcher_t thread;

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
extern pipeline_dispatcher_t spu;
#endif

END_C_DECLS

#endif // !CELL_SPU_BUILD

#endif // _pipelines_h_ 
