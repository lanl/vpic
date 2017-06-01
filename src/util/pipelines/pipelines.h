#ifndef _pipelines_h_
#define _pipelines_h_

#include "../util_base.h"

enum { MAX_PIPELINE = 64 };

// A pipeline function takes a pointer to arguments for the pipeline
// and a integer which gives the rank of the pipeline and the total
// number of pipelines dispatched.

typedef void
(*pipeline_func_t)( void * args,
                    int pipeline_rank,
                    int n_pipeline );

///////////////////////////////////////////////////////////////////////////////

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

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

typedef struct pipeline_dispatcher {

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

#endif // _pipelines_h_
