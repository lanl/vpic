#ifndef _pipelines_h_
#define _pipelines_h_

#include <util_base.h>

enum { MAX_PIPELINE = 16 };

#if !defined(CELL_SPU_BUILD)
// FIXME: Should all these really be protected from a CELL_SPU_BUILD?

// A pipeline function takes a pointer to arguments for the pipeline
// and a integer which gives the rank of the pipeline and the total
// number of pipelines dispatched.

typedef void
(*pipeline_func_t)( void * args,
                    int pipeline_rank,
                    int n_pipeline );

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
  // size_args gives the byte size of an element of the argument
  // array.  Use 0 if you want to pass the exact same arguments to
  // each pipeline.
  //
  // If the pipeline functions do not take arguments, use NULL for
  // args and 0 for the size_args.
                    
  void
  (*dispatch)( pipeline_func_t pipeline,
               void * args,
               int size_args );

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
