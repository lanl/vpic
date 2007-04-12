/* NOTES: THESE HOST->PIPELINE COMMUNICATION INTERFACE STUBS
   ARE LOOSELY BASED ON MY THREADER LIBRARY.  IT SHOULD BE
   REPLACED BY A ROADRUNNER CELL SPECIFIC PPC/SPE INTERFACE.
   THIS ALLOWS FOR TESTING OF THE PIPELINES WITHOUT HAVING
   A CELL PROCESSOR. */

#ifndef _pipelines_h_
#define _pipelines_h_

#include <util_base.h>

enum { MAX_PIPELINE = 16 };

/* A pipeline function takes a pointer to arguments for the
   pipeline and a integer which gives the rank of the pipeline
   and the total number of pipelines dispatched. */

typedef void
(*pipeline_func_t)( void * args,
                    int pipeline_rank,
                    int n_pipeline );

typedef struct pipeline_dispatcher {

  /* n_pipelines indicates the number of pipelines currently
     running.  Technically, this should be read only for
     users! */

  int n_pipeline;

  /* boot creates the number of pipelines requested.
     Generally, this is number of cores on a node if using symmetric
     multiprocessing (or number of cores on a node if you want one core
     dedicated to the host tasks) or the number of pipeline processors
     if using asymmetric multiprocessing. */

  void (*boot)( int n_pipelines_requested );

  /* halt destroys all the pipelines created by boot. */

  void (*halt)( void );

  /* dispatch begins executing the given pipeline function on all
     the pipelines.

     pipeline is the pipeline function to execute on the pipelines.
   
     args is an array of arguments to pass to each pipeline.
   
     size_args gives the byte size of an element of the argument array.
     Use 0 if you want to pass the exact same arguments to each pipeline.
   
     request stores information for the host processor about completion of
     the request.  Use NULL if the host does not need to know about completion.
   
     If the pipeline functions do not take arguments, use NULL for
     args and 0 for the size_args. */
                    
  void
  (*dispatch)( pipeline_func_t pipeline,
               void * args,
               int size_args );

  /* wait waits for the previous dispatch to complete. */

  void
  (*wait)( void );

} pipeline_dispatcher_t;

/* Temporary hack until functions proper starting calling the
   appropriate kinds of dispatchers for their pipelines. */

#if 1

#define _n_pipeline ((int)serial.n_pipeline)
#define boot_pipelines serial.boot
#define halt_pipelines serial.halt
#define dispatch_pipelines(p,a,s) serial.dispatch((pipeline_func_t)(p),(a),(s))
#define wait_for_pipelines serial.wait

#else

#define _n_pipeline ((int)thread.n_pipeline)
#define boot_pipelines thread.boot
#define halt_pipelines thread.halt
#define dispatch_pipelines(p,a,s) thread.dispatch((pipeline_func_t)(p),(a),(s))
#define wait_for_pipelines thread.wait

#endif

BEGIN_C_DECLS

extern pipeline_dispatcher_t serial;
extern pipeline_dispatcher_t thread;
extern pipeline_dispatcher_t spu;

END_C_DECLS

#endif /* _pipelines_h_ */
