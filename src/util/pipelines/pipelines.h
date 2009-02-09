#ifndef _pipelines_h_
#define _pipelines_h_

#include "../util_base.h"

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
// FIXME: Should CELL_SPU_BUILDs include this (and if so, under what
// conditions)?
#include <libspe2.h> // For spe_program_handle_t
#endif

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
                    
#define register_function(name) \
  void (*name)(void *) = thread_##name

  void
  (*dispatch)( pipeline_func_t pipeline,
               void * args,
               int size_args );

  // wait waits for the previous dispatch to complete.

  void
  (*wait)( void );

  // signal SPEs using mailboxes
  void
  (*signal)(uint32_t signal);

  // synchronize on signal from SPEs using mailboxes
  void
  (*sync)(uint32_t signal);

} pipeline_dispatcher_t;


BEGIN_C_DECLS

extern pipeline_dispatcher_t serial; // For debugging purposes
extern pipeline_dispatcher_t thread;

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)

extern pipeline_dispatcher_t spu;

// Due to the vagaries of Cell, spu.dispatch has some restrictions
// that the other pipelines do not have.
//
// Usage:
//
//   spu.dispatch( SPU_PIPELINE(pipeline), args, size );
//
// - The spu pipeline (a spe_program_handle_t) must be encapsulated in 
//   the SPU_PIPELINE macro for reasons of evil discussed below.
// - args should be 16 byte aligned (ideally 128)
// - size should be a multiple of 16 (ideally 128), 0 is permissible
//
// The spu pipeline itself sees:
//
//   int main( uint64_t spu_id, uint64_t args, uint64_t envp ) {
//     ... my spu pipeline ...
//   }
//
// where:
//   args = (uint64_t)user_args + job * size
// and:
//   envp = job in the low 32 bits and n_job in high 32-bit
//
// Note on the evil: The spu.dispatch pipeline wants a
// spe_program_handle_t * (an object pointer) instead of
// pipeline_func_t (a function pointer type) like the other
// dispatchers.  C forbids casting directly been object pointer and
// function pointer types.  The below evil macro:
//
// (1) Takes a "spe_program_handle_t".  (Not a pointer to one.)
//
// (2) Stops anybody from passing a typical function pointer to the
// spu dispatcher as it will not be directly castable to a
// "spe_program_handle_t *"
//
// (3) Casts indirectly a "spe_program_handle_t *" into a
// "pipeline_func_t."

#define register_function(name) \
	extern uint32_t spe_##name; \
	uint32_t name = spe_##name

#define SPU_PIPELINE(x) \
  ((pipeline_func_t)(uint64_t)(spe_program_handle_t *)&(x))

#endif

END_C_DECLS

#endif // !CELL_SPU_BUILD

#endif // _pipelines_h_ 
