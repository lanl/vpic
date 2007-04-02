#ifndef _pipeline_h_
#define _pipeline_h_

#include <common.h>

/* NOTES: THESE HOST->PIPELINE COMMUNICATION INTERFACE STUBS
   ARE LOOSELY BASED ON MY THREADER LIBRARY.  IT SHOULD BE
   REPLACED BY A ROADRUNNER CELL SPECIFIC PPC/SPE INTERFACE.
   THIS ALLOWS FOR TESTING OF THE PIPELINES WITHOUT HAVING
   A CELL PROCESSOR. */

enum {
  MAX_PIPELINE = 1,
  n_pipeline   = MAX_PIPELINE
};

/* A pipeline function takes a pointer to arguments for the
   pipeline and a integer which gives the rank of the pipeline. */

typedef void
(*pipeline_func_t)( void * args,
                    int pipeline_rank );

/* A pipeline request is used to monitor the pipeline status */

typedef struct pipeline_request {
  volatile int complete[MAX_PIPELINE];
} pipeline_request_t;

BEGIN_C_DECLS

/* dispatch_pipelines begins executing the given pipeline function on all
   the pipelines.

   pipeline is the pieline function to execute on the pipelines.
   
   args is an array of arguments to pass to each pipeline.
   
   size_args gives the byte size of an element of the argument array.
   Use 0 if you want to pass the exact same arguments to each pipeline.
   
   request stores information for the host processor about completion of
   the request.  Use NULL if the host does not need to know about completion.
   
   If the pipeline functions do not take arguments, use NULL for
   args and 0 for the size_args. */
                    
void
_dispatch_pipelines( pipeline_func_t pipeline,
                     void * args,
                     int size_args,
                     pipeline_request_t * request );
                     
#define dispatch_pipelines(p,a,s,r) \
  _dispatch_pipelines((pipeline_func_t)(p),(a),(s),(r))

/* wait_for_pipelines waits for all the given pipeline_request to
   complete. */
   
void
wait_for_pipelines( pipeline_request_t * request );

END_C_DECLS

#endif /* _pipeline_h_ */
