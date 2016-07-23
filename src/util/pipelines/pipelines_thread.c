/****************************************************************************
 *
 * Copyright (c) 2000, Kevin James Bowers
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The names of its contributors may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ***************************************************************************/

// Modified extensively for VPIC-3P by K. Bowers 4/10/2007

// TO DO: 
// (?) Signal blocking in pipelines
// (?) Timeouts in thread_halt, thread_boot (spin wait)

#include <pthread.h>

#include "pipelines.h"

static void *
pipeline_mgr( void *_id );

enum {
  PIPELINE_SLEEP = 0,
  PIPELINE_EXECUTE = 1,
  PIPELINE_TERMINATE = 2,
  PIPELINE_ACK = 3
};

typedef struct pipeline_state {
  pthread_t handle;
  pthread_mutex_t mutex;
  pthread_cond_t wake;
  volatile int state;
  pipeline_func_t func;
  void *args;
  int job;
  int n_job;
  volatile int *flag;
} pipeline_state_t;

static pthread_t Host;
static pipeline_state_t Pipeline[ MAX_PIPELINE ];
static volatile int Done[ MAX_PIPELINE ];
static int Id = 0;
static int Busy = 0;
static int Dispatch_To_Host = 0;

/****************************************************************************/

#include "../checkpt/checkpt.h"

void
checkpt_thread( const pipeline_dispatcher_t * _thread ) {
  CHECKPT_VAL( int, thread.n_pipeline );
}

pipeline_dispatcher_t *
restore_thread( void ) {
  int n_pipeline;
  RESTORE_VAL( int, n_pipeline );
  if( thread.n_pipeline!=n_pipeline )
    ERROR(( "--tpp changed between checkpt (%i) and restore (%i)",
            n_pipeline, thread.n_pipeline ));
  return &thread;
}

/****************************************************************************
 *
 * Thread dispatcher notes:
 * - Only one instance of thread dispatcher can exist at time (i.e. the
 *   application can't spawn two threads and have each thread create its
 *   own thread dispatcher).  Otherwise bad things could happen.
 * - The dispatcher assumes pipelines always return (otherwise a
 *   possible deadlock exists in parallel_execute and thread_halt _will_
 *   deadlock).
 * - The dispatcher intends that only host will call parallel_execute.
 *   (Original version of the dispatcher was less strict and most of that
 *   infrastructure is still present though.)
 * - Calls to thread_boot should always be paired with calls to
 *   thread_halt (this insures good housekeeping). To reallocate the
 *   number of pipelines available to the program, the dispatcher should
 *   be halted and then re-booted with the new desired number of
 *   pipelines.
 * - Internal: thread.n_pipeline = 0 -> Global variables need to be
 *   initialized. No other global variables should be trusted when
 *   thread.n_pipeline = 0.
 *
 ***************************************************************************/

/****************************************************************************
 *
 * Function: thread_boot( n_pipeline, dispatch_to_host )
 *
 * Initialize the thread dispatcher to use n_pipeline threads.
 * The underlying thread library may create extra management threads.
 * This routine does not require you to have the requested number of
 * physical processors in your system (the host and/or pipeline threads
 * will be multitasked on physical processors in that case).
 *
 * If the boolean dispatch_to_host is true, then the host thread
 * (i.e. the caller of this function) will be treated as a pipeline
 * thread.  Thus, n_pipeline-1 pipeline threads will be created.
 * thread_dispatch will execute the workload given to pipelines
 * 0:n_pipeline-2 on the n_pipeline-1 pipeline threads and the workload
 * given to pipeline n_pipeline-1 on the host thread.
 *
 * If the boolean dispatch_to_host is false, then the host thread
 * will _not_ be treaded as a pipeline thread.  Thus, n_pipeline
 * pipeline threads will be created.  thread_dispatch will execute all
 * n_pipeline workloads on the n_pipeline threads (and the host thread
 * will presumably do other things while the pipeline threads are
 * executing).
 *
 * Arguments: n_pipeline (int) - The desired number of pipelines.
 *            dispatch_to_host (int) - Should the dispatch function
 *              dispatch tasks to the host thread 
 *
 * Returns: (void)
 *
 * Throws: Throws an error if
 *         * if n_pipeline is out of bounds. n_pipeline should satisfy
 *           1 <= n_pipeline <= MAX_PIPELINE
 *         * if the dispatcher was already initialized.
 *         * if something strange happened during global
 *           initialization or pipeline creation.
 *
 * Dependencies: pthread_self, pthread_cond_init, pthread_mutex_init
 *               pthread_create
 *
 * Globals read/altered: thread.n_pipeline (rw), Host (rw), Id (rw),
 *                       Pipeline (rw), Busy (rw), Dispatch_To_Host (rw)
 *
 * Notes:
 * - This should only be called once by at most one thread in a program.
 * - The calling thread becomes the Host thread.
 * - Calls to thread_boot should always be paired with calls to
 *   thread_halt for good housekeeping.
 *
 ***************************************************************************/

static void
thread_boot( int * pargc,
             char *** pargv ) {
  int i, n_pipeline;

  // Check if arguments are valid and dispatcher isn't already initialized

  if( thread.n_pipeline != 0 ) ERROR(( "Halt the thread dispatcher first!" ));

  // # if defined(CELL_PPU_BUILD)
  //   n_pipeline       = strip_cmdline_int( pargc, pargv, "--tpp",              2 );
  // # else
  //   n_pipeline       = strip_cmdline_int( pargc, pargv, "--tpp",              1 );
  // # endif

  n_pipeline       = strip_cmdline_int( pargc, pargv, "--tpp",              1 );
  Dispatch_To_Host = strip_cmdline_int( pargc, pargv, "--dispatch_to_host", 1 );

  if( n_pipeline<1 || n_pipeline>MAX_PIPELINE )
    ERROR(( "Invalid number of pipelines requested (%i)", n_pipeline ));

  // Initialize some global variables. Note: thread.n_pipeline = 0 here 

  Id   = 0;
  Busy = 0;
  Host = pthread_self();

  // Initialize all the pipelines

  for( i=0; i<n_pipeline-Dispatch_To_Host; i++ ) {

    // Initialize the pipeline state. Note: PIPELINE_ACK will be
    // cleared once the thread control function pipeline_mgr starts
    // executing and the thread is ready to execute pipelines.

    Pipeline[i].state = PIPELINE_ACK;

    // Initialize the pipeline mutex and signal condition variables
    // and spawn the pipeline.  Note: When mutexes are initialized,
    // they are initialized to the unlocked state

    if( pthread_cond_init( &Pipeline[i].wake, NULL ) )
      ERROR(( "pthread_cond_init failed" ));

    if( pthread_mutex_init( &Pipeline[i].mutex, NULL ) )
      ERROR(( "pthread_mutex_init failed" ));

    if( pthread_create( &Pipeline[i].handle,
                        NULL,
                        pipeline_mgr,
                        &Pipeline[i] ) )
      ERROR(( "pthread_create failed" ));

    // The thread spawn was successful.  Wait for the thread to
    // acknowledge and go to sleep.  This insures all created threads
    // will be ready to execute pipelines upon return from this
    // function.  Note: If the thread never acknowledges this function
    // will spin forever (maybe should add a time out) - however, if
    // this happens, it is not the library's fault (O/S or pthreads
    // screwed up)

    while( Pipeline[i].state != PIPELINE_SLEEP ) nanodelay(50);

  }

  thread.n_pipeline = n_pipeline;
  REGISTER_OBJECT( &thread, checkpt_thread, restore_thread, NULL );
}

/****************************************************************************
 *
 * Function: thread_halt
 *
 * Terminates all the pipeline threader and frees system resources
 * associate with the thread
 *
 * Arguments: None
 *
 * Returns: None
 * 
 * Throws: An error if:
 *         - The caller is not the host thread
 *         - Errors occurred terminating pipeline threads
 *         - Errors occurred freeing pipeline resources
 *         - Errors occurred freeing dispatcher resources
 *
 * Dependencies: pthread_mutex_lock, pthread_mutex_unlock, pthread_join
 *               pthread_cond_signal, pthread_mutex_destroy,
 *               pthread_cond_destroy
 *
 * Globals read/altered: thread.n_pipeline (r), Busy (r), Pipeline (rw)
 *
 * Notes:
 * - This function will spin forever if a thread is executing a
 *   pipeline which doesn't return.  Thus this library assumes pipelines
 *   _always_ return.
 *
 ***************************************************************************/

static void
thread_halt( void ) {
  int id;

  // Make sure the host is the caller and there are pipelines to
  // terminate

  if( !thread.n_pipeline ) ERROR(( "Boot the thread dispatcher first!" ));

  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may halt the thread dispatcher!" ));

  if( Busy ) ERROR(( "Pipelines are busy!" ));

  UNREGISTER_OBJECT( &thread );

  // Terminate the threads

  for( id=0; id<thread.n_pipeline-Dispatch_To_Host; id++ ) {

    // Signal the thread to terminate. Note: If the thread is
    // executing a pipeline which doesn't return, this function will
    // spin forever here (maybe should add a timeout, forced kill)

    pthread_mutex_lock( &Pipeline[id].mutex );
    Pipeline[id].state = PIPELINE_TERMINATE;
    pthread_cond_signal( &Pipeline[id].wake );
    pthread_mutex_unlock( &Pipeline[id].mutex );

    // Wait for the thread to terminate. Note: As Pipeline[id].state =
    // PIPELINE_TERMINATE now, non-terminated threads calling
    // parallel_execute will not be able to dispatch pipelines to
    // terminated pipelines (see the innermost switch in
    // parallel_execute).

    if( pthread_join( Pipeline[id].handle, NULL ) )
      ERROR(( "Unable to terminate pipeline!" ));

  }

  // Free resources associated with Threader as all pipelines are now
  // dead.  Note: This must be done in a separate loop because
  // non-terminated pipelines calling parallel_execute may try to
  // access mutexes for destroyed pipelines and Id and with unknown
  // results if the mutexes were destroyed before all the pipelines
  // were terminated.

  for( id=0; id<thread.n_pipeline-Dispatch_To_Host; id++ )
    if( pthread_mutex_destroy( &Pipeline[id].mutex ) ||
        pthread_cond_destroy( &Pipeline[id].wake ) )
      ERROR(( "Unable to destroy pipeline resources!" ));

  // Finish up

  thread.n_pipeline = 0;
}

/****************************************************************************
 *
 * Function: parallel_execute(pipeline,params,job,n_job,flag)
 *
 * Sends a thread to work on a given pipeline.
 *
 * Arguments: pipeline (pipeline_func_t) - Pointer to function to be executed
 *            params (void *) - Pointer to parameters to pass to the function
 *            job (int) - Identifier of the pipeline job
 *            n_job (int) - Identifier of total number of jobs
 *            flag (volatile int *) - Pointer to an int used to signal the
 *              host the pipeline is complete.  NULL if the host does not
 *              need notification.
 *
 * Returns: nothing
 *
 * Throws: An error if:
 *         - No threads exist to handle the pipeline
 *         - Parallel_execute detected thread_halt terminating pipelines.
 *         (The host DOES NOT complete the pipeline and the pipeline was
 *         NOT dispatched to any thread in error cases)
 *
 * Dependencies: pthread_mutex_unlock, pthread_self, pthread_mutex_trylock,
 *               pthread_equal, pthread_mutex_lock, pthread_cond_signal
 *
 * Globals read/altered: thread.n_pipeline (r), Id, Dispatch_To_Host,
 *                       Pipeline (rw)
 *
 * Notes:
 * - Only the host may call this.
 * - If all threads are busy, the function will wait until a thread
 *   becomes available to execute the pipeline.  There is a possibility
 *   of a deadlock here if the none of the other pipelines return.
 *   Thus, this library assumes pipelines _always_ eventually return.
 * - It is remotely possible for parallel_execute to reacquire the
 *   lock on a thread signalled to execute a pipeline before the
 *   thread wakes up. However, this is not a problem because the
 *   Pipeline[id].state will prevent the pipeline from being signaled a
 *   second time (see innermost switch).
 *
 ***************************************************************************/

static void
parallel_execute( pipeline_func_t func,
                  void *args,
                  int job,
                  int n_job,
                  volatile int *flag ) {
  int id;

  // Determine that the thread dispatcher is initialized and that
  // caller is the host thread.

  if( thread.n_pipeline==0 )
    ERROR(( "Boot the thread dispatcher first!" ));

  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may call parallel_execute" ));

  if( thread.n_pipeline-Dispatch_To_Host<=0 )
    ERROR(( "No threads available for parallel execution!" ));

  // If pipelines are potentially available to execute tasks, loop
  // until we dispatch the task to a pipeline (or we detect
  // thread_halt is running)

  for(;;) {

    // Get an id of a pipeline to query.

    id = Id; if( (++Id) >= thread.n_pipeline-Dispatch_To_Host ) Id = 0;

    if( !pthread_mutex_trylock( &Pipeline[id].mutex ) ) {
      
      // If the thread isn't executing, see if the job can be
      // dispatched to this pipeline.  Note: host now has control over
      // the thread's state.
      
      switch( Pipeline[id].state ) {
        
        // Note: all cases relinquish control over the pipeline's
        // state
        
      case PIPELINE_SLEEP:
        
        // The pipeline is available - assign the task, wake the
        // pipeline and return
        
        Pipeline[id].state = PIPELINE_EXECUTE;
        Pipeline[id].func  = func;
        Pipeline[id].args  = args;
        Pipeline[id].job   = job;
        Pipeline[id].n_job = n_job;
        Pipeline[id].flag  = flag;
        pthread_cond_signal( &Pipeline[id].wake );
        pthread_mutex_unlock( &Pipeline[id].mutex );
        return;
        
      case PIPELINE_TERMINATE:
        
        // The query pipeline was terminated.  Something is amiss!
        
        pthread_mutex_unlock( &Pipeline[id].mutex );
        ERROR(( "parallel_execute called while thread_halt is running - "
                "Should never happen in this restricted implementation." ));
        
      default:
        
        // Pipeline isn't ready to accept new tasks (pipeline is about
        // to wake-up to execute an already assigned task)
        
        pthread_mutex_unlock( &Pipeline[id].mutex );
        
      } // switch
    } // if pipeline might be ready
  } // for
  // Never get here
}

/****************************************************************************
 *
 * Function: pipeline_mgr(params) - internal use only
 *
 * The pipeline side of the parallel_execute
 *
 * Arguments: params (void *) - pointer to the pipeline to manage
 *
 * Returns: (void *) NULL always.
 *
 * Dependencies: pthread_mutex_lock, pthread_mutex_unlock, pthread_cond_wait
 *
 * Globals read/altered: Pipeline (rw)
 *
 * Notes:
 * - When pipeline is not sleeping (PIPELINE_SLEEP state), it has a
 *   mutex lock over over Pipeline[id].state to insure consistency in
 *   the loop flow control
 *
 ***************************************************************************/

static void *
pipeline_mgr( void *_pipeline ) {
  pipeline_state_t *pipeline = _pipeline;

  // Pipeline state is PIPELINE_ACK and the pipeline mutex is unlocked
  // when entering.  Since pthread_cond_wait unlockes the pipeline
  // mutex when the pipeline goes to sleep and the PIPELINE_ACK case
  // writes the pipeline state, the mutex needs to be locked first.
     
  pthread_mutex_lock( &pipeline->mutex );

  // Loop while the pipeline is still executing

  for(;;) {

    switch( pipeline->state ) {

    case PIPELINE_TERMINATE:

      // Terminate the pipeline. Unlock the mutex first (as
      // pthread_cond_wait locks it when the pipeline wakes up).

      pthread_mutex_unlock( &pipeline->mutex );
      return NULL;

    case PIPELINE_EXECUTE:

      // Execute the given task and set the completion flag if
      // necessary.  Note: the pipeline mutex is locked while the
      // pipeline is executing a task.

      if( pipeline->func )
        pipeline->func( pipeline->args, pipeline->job, pipeline->n_job );
      if( pipeline->flag ) *pipeline->flag = 1;

      // Pass through into the next case

    case PIPELINE_ACK:
    case PIPELINE_SLEEP:
    default:

      // Go to sleep. Note: pthread_cond_wait unlocks the pipeline
      // mutex while the pipeline is sleeping and locks it when the
      // pipeline wakes up

      pipeline->state = PIPELINE_SLEEP;
      pthread_cond_wait( &pipeline->wake, &pipeline->mutex );

      break;

    } // switch

  } // for

  return NULL; // Never get here - Avoid compiler warning
}

static void
thread_dispatch( pipeline_func_t func,
                 void * args,
                 int sz,
                 int str ) {
  int id;

  if( thread.n_pipeline==0 ) ERROR(( "Boot the thread dispatcher first!" ));

  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may call thread_dispatch" ));

  if( Busy ) ERROR(( "Pipelines are busy!" ));
  Busy = 1;

  for( id=0; id<thread.n_pipeline-Dispatch_To_Host; id++ ) {
    Done[id] = 0;
    parallel_execute( func,
                      ((char *)args) + id*sz*str,
                      id,
                      thread.n_pipeline,
                      &Done[id] );
  }

  if( Dispatch_To_Host ) {
    Done[id] = 0;
    if( func ) func( ((char *)args) + id*sz*str, id, thread.n_pipeline );
    Done[id] = 1;
  }
}
                 
static void
thread_wait( void ) {
  int id;

  if( thread.n_pipeline==0 ) ERROR(( "Boot the thread dispatcher first!" ));

  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may call thread_wait" ));

  if( !Busy ) ERROR(( "Pipelines are not busy!" ));

  id = 0;
  while( id<thread.n_pipeline ) {
    if( Done[id] ) id++;
    else           nanodelay(5);
  }

  Busy = 0;
}

pipeline_dispatcher_t thread = {
  0,               // n_pipeline
  thread_boot,     // boot
  thread_halt,     // halt
  thread_dispatch, // dispatch
  thread_wait      // wait
};

