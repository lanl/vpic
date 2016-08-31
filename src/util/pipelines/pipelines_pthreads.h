#ifndef _pipelines_pthreads_h_
#define _pipelines_pthreads_h_

#ifndef THREAD_REROUTE
#error "Do not include pipelines_pthreads.h directly; use pipelines.h."
#endif

#include "../util_base.h"

enum { MAX_PIPELINE = 272 };

// A pipeline function takes a pointer to arguments for the pipeline
// and a integer which gives the rank of the pipeline and the total
// number of pipelines dispatched.

typedef void
(*pipeline_func_t)( void * args,
                    int pipeline_rank,
                    int n_pipeline );

//----------------------------------------------------------------------------//
// Generic macros that are used for all cases of vector acceleration as well
// as the standard case that does not use vector acceleration.
//----------------------------------------------------------------------------//

# define N_PIPELINE thread.n_pipeline

# define WAIT_PIPELINES() thread.wait()

# define PAD_STRUCT( sz )

//----------------------------------------------------------------------------//
// Macro defines to support v16 simd vector acceleration.  Uses thread
// dispatcher on the v16 pipeline and the caller does straggler cleanup with
// the scalar pipeline.
//----------------------------------------------------------------------------//

#if defined(V16_ACCELERATION) && defined(HAS_V16_PIPELINE)

# define EXEC_PIPELINES(name,args,str)                           \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v16,         \
                   args, sizeof(*args), str );                   \
  name##_pipeline( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )

# define PROTOTYPE_PIPELINE( name, args_t ) \
  void                                      \
  name##_pipeline_v16( args_t *args,        \
                       int pipeline_rank,   \
                       int n_pipeline );    \
                                            \
  void                                      \
  name##_pipeline( args_t *args,            \
                   int pipeline_rank,       \
                   int n_pipeline )

//----------------------------------------------------------------------------//
// Macro defines to support v8 simd vector acceleration.  Uses thread
// dispatcher on the v8 pipeline and the caller does straggler cleanup with
// the scalar pipeline.
//----------------------------------------------------------------------------//

#elif defined(V8_ACCELERATION) && defined(HAS_V8_PIPELINE)

# define EXEC_PIPELINES(name,args,str)                           \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v8,          \
                   args, sizeof(*args), str );                   \
  name##_pipeline( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )

# define PROTOTYPE_PIPELINE( name, args_t ) \
  void                                      \
  name##_pipeline_v8( args_t *args,         \
                      int pipeline_rank,    \
                      int n_pipeline );     \
                                            \
  void                                      \
  name##_pipeline( args_t *args,            \
                   int pipeline_rank,       \
                   int n_pipeline )

//----------------------------------------------------------------------------//
// Macro defines to support v4 simd vector acceleration.  Uses thread
// dispatcher on the v4 pipeline and the caller does straggler cleanup with
// the scalar pipeline.
//----------------------------------------------------------------------------//

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

# define EXEC_PIPELINES(name,args,str)                           \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v4,          \
                   args, sizeof(*args), str );                   \
  name##_pipeline( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )

# define PROTOTYPE_PIPELINE( name, args_t ) \
  void                                      \
  name##_pipeline_v4( args_t *args,         \
                      int pipeline_rank,    \
                      int n_pipeline );     \
                                            \
  void                                      \
  name##_pipeline( args_t *args,            \
                   int pipeline_rank,       \
                   int n_pipeline )

//----------------------------------------------------------------------------//
// Macro defines to support the standard implementation which does not use
// explicit simd vectorization.  Uses thread dispatcher on the scalar pipeline
// and the caller does straggler cleanup with the scalar pipeline.
//----------------------------------------------------------------------------//

#else

# define EXEC_PIPELINES(name,args,str)                           \
  thread.dispatch( (pipeline_func_t)name##_pipeline,             \
                   args, sizeof(*args), str );                   \
  name##_pipeline( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )

# define PROTOTYPE_PIPELINE( name, args_t ) \
  void                                      \
  name##_pipeline( args_t *args,            \
                   int pipeline_rank,       \
                   int n_pipeline )

#endif

////////////////////////////////////////////////////////////////////////////////

typedef struct pipeline_dispatcher
{
  // n_pipelines indicates the number of pipelines currently running.
  // Technically, this should be read only for users!

  int n_pipeline;

  // boot creates the number of pipelines requested (in the command
  // line args).  Generally, this is number of cores on a node if
  // using symmetric multiprocessing or the number of pipeline
  // processors if using heterogeneous multiprocessing.

  void
  (*boot)( int * pargc,
           char *** pargv );

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

END_C_DECLS

#endif // _pipelines_pthreads_h_ 
