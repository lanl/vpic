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

/* Modified extensively for VPIC-3P by K. Bowers 4/10/2007 */
/* Modified extensively for VPIC-3P SPU dispatch by K. Bowers 4/20/2007 */

/* TO DO: 
   (?) Signal blocking in spu_control_threads
   (?) Timeouts in spu_halt, spu_boot (spin wait) */

#include <pipelines.h>
#include <pthread.h>
#include <libspe2.h>

static void *
spu_control_thread( void *_id );

enum {
  SPU_CONTROL_THREAD_SLEEP = 0,
  SPU_CONTROL_THREAD_EXECUTE = 1,
  SPU_CONTROL_THREAD_TERMINATE = 2,
  SPU_CONTROL_THREAD_ACK = 3
};

typedef struct spu_control_state {
  pthread_t handle;
  pthread_mutex_t mutex;
  pthread_cond_t wake;
  volatile int state;
  spe_context_ptr_t context;
  spe_program_handle_t * pipeline;
  void *args;
  int job;
  int n_job;
  volatile int *flag;
} spu_control_state_t;

static pthread_t Host;
static spu_control_state_t SPU_Control_State[ MAX_PIPELINE ];
static volatile int Done[ MAX_PIPELINE ];
static int Id = 0;
static int Busy = 0;

/****************************************************************************
 *
 * SPU pipeline dispatcher notes:
 * - Only one instance of spu pipeline dispatcher can run at time (i.e. the
 *   application can't spawn two PPU threads and have each thread create its
 *   own dispatcher).  Otherwise bad things could happen.
 * - The dispatcher assumes pipelines always terminate (otherwise a
 *   possible deadlock exists in parallel_execute and spu_halt _will_
 *   deadlock).
 * - The dispatcher intends that only host PPU thread will call
 *   parallel_execute.  (Original version of threader was less strict and
 *   most of that infrastructure is still present though.)
 * - Calls to spu_boot should always be paired with calls to spu_halt
 *   (this insures good housekeeping). To reallocate the number of
 *   pipelines available to the host PPU thread, the dispatcher should
 *   be halted and then re-booted with the new desired number of
 *   pipelines.
 * - Internal: spu.n_pipeline = 0 -> Global variables need to be
 *   initialized. No other global variables should be trusted when
 *   spu.n_pipeline = 0.
 *
 ***************************************************************************/

/****************************************************************************
 *
 * Function: spu_boot(num_pipe)
 *
 * Initialize the dispatcher to have num_pipe spu control threads in
 * addition to the PPU host thread.  This routine does not require you
 * to have the requested number of physical SPUs in your system
 * (libspe2 will handle the multitasking of jobs in that case).
 *
 * Arguments: num_pipe (int) - The desired number of pipelines.
 *
 * Returns: (void)
 *
 * Throws: Throws an error if
 *         * if num_pipe is out of bounds. num_pipe should satisfy
 *           1 <= numproc <= MAX_PIPELINES
 *         * if the dispatcher was already initialized.
 *         * if something strange happened during global
 *           initialization or spu control thread creation.
 *
 * Dependencies: pthread_self, pthread_cond_init, pthread_mutex_init
 *               pthread_create
 *
 * Globals read/altered: spu.n_pipeline (rw), Host (rw),
 *                       Pipeline (rw)
 *
 * Notes:
 * - This should only be called once by at most one PPU thread in a program.
 * - The calling PPU thread becomes the PPU Host thread.
 * - Calls to spu_boot should always be paired with calls to
 *   spu_halt for good housekeeping.
 *
 ***************************************************************************/

static void
spu_boot( int num_pipe ) {
  int i;

  /* Check if arguments are valid and dispatcher isn't already initialized */

  if( spu.n_pipeline != 0 ) ERROR(( "Halt the spu dispatcher first!" ));
  if( num_pipe < 1 || num_pipe > MAX_PIPELINE )
    ERROR(( "Invalid number of pipelines requested" ));

  /* Initialize some global variables. Note: spu.n_pipeline = 0 here */ 

  Host = pthread_self();
  Id = 0;

  /* Initialize all the pipelines */

  for( i=0; i<num_pipe; i++ ) {

    /* Initialize the spu control thread state. Note:
       SPU_CONTROL_THREAD_ACK will be cleared once the spu control
       thread starts executing and the spu control thread is ready to
       execute pipelines. */

    SPU_Control_State[i].state = SPU_CONTROL_THREAD_ACK;

    /* Initialize the spu control thread mutex, signal condition
       variables and spe context and spawn the spu control thread.
       Note: When mutexes are initialized, they are initialized to the
       unlocked state */

    if( pthread_cond_init( &SPU_Control_State[i].wake, NULL ) )
      ERROR(( "pthread_cond_init failed" ));
    
    if( pthread_mutex_init( &SPU_Control_State[i].mutex, NULL ) )
      ERROR(( "pthread_mutex_init failed" ));

    if( ( SPU_Control_State[i].context=spe_context_create( 0, NULL ) )==NULL )
      ERROR(( "spe_context_create failed" ));

    if( pthread_create( &SPU_Control_State[i].handle,
                        NULL,
                        spu_control_thread,
                        &SPU_Control_State[i] ) )
      ERROR(( "pthread_create failed" ));

    /* The spu control thread spawn was successful.  Wait for the spu
       control thread to acknowledge and go to sleep.  This insures
       all created spu control threads will be ready to execute
       pipelines upon return from this function.  Note: If the spu
       control thread never acknowledges this function will spin
       forever (maybe should add a time out) - however, if this
       happens, it is not the library's fault (O/S or pthreads screwed
       up) */

    while( SPU_Control_State[i].state != SPU_CONTROL_THREAD_SLEEP )
      nanodelay(50);

  }

  spu.n_pipeline = num_pipe;
  Busy = 0;
}

/****************************************************************************
 *
 * Function: spu_halt
 *
 * Terminates all the spu control threads and frees all system resources
 * associated with the dispatcher.
 *
 * Arguments: None
 *
 * Returns: None
 * 
 * Throws: An error if:
 *         - The caller is not the host PPU thread
 *         - Errors occurred terminating SPU control threads
 *         - Errors occurred freeing dispatcher resources
 *
 * Dependencies: pthread_mutex_lock, pthread_mutex_unlock, pthread_join
 *               pthread_cond_signal, pthread_mutex_destroy,
 *               pthread_cond_destroy
 *
 * Globals read/altered: spu.n_pipeline (r), SPU_Control_State (rw)
 *
 * Notes:
 * - This function will spin forever if a spu control thread is
 *   executing a pipeline which doesn't terminate.  Thus this library
 *   assumes pipelines _always_ terminate.
 *
 ***************************************************************************/

static void
spu_halt( void ) {
  int id;

  /* Make sure the host is the caller and there are spu control threads to
     terminate */

  if( !spu.n_pipeline ) ERROR(( "Boot the spu dispatcher first!" ));
  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may halt the spu dispatcher!" ));

  /* Terminate the spu control threads */

  for( id=0; id<spu.n_pipeline; id++ ) {

    /* Signal the spu control thread to terminate. Note: If the spu
       control thread is executing a pipeline which doesn't return,
       this function will spin forever here (maybe should add a
       timeout, forced kill) */

    pthread_mutex_lock( &SPU_Control_State[id].mutex );
    SPU_Control_State[id].state = SPU_CONTROL_THREAD_TERMINATE;
    pthread_cond_signal( &SPU_Control_State[id].wake );
    pthread_mutex_unlock( &SPU_Control_State[id].mutex );

    /* Wait for the spu control thread to terminate. Note: As
       SPU_Control_State[id].state = SPU_CONTROL_THREAD_TERMINATE
       now, non-terminated spu control threads calling
       parallel_execute will not be able to dispatch pipelines to
       terminated spu_control_threads (see the innermost switch in
       parallel_execute) (This note does not apply for this
       restricted dispatcher.) */

    if( pthread_join( SPU_Control_State[id].handle, NULL ) )
      ERROR(( "Unable to terminate spu control thread!" ));

  }

  /* Free resources associated with dispatcher as all spu control
     threads are now dead.  Note: This must be done in a separate loop
     because non-terminated spu control threads calling
     parallel_execute may try to access mutexes for destroyed spu
     control threads and Id and with unknown results if the mutexes
     were destroyed before all the spu control threads were
     terminated.  (This note does not apply for this restricted
     dispatcher.) */

  for( id=0; id<spu.n_pipeline; id++ )
    if( spe_context_destroy( SPU_Control_State[id].context ) ||
        pthread_mutex_destroy( &SPU_Control_State[id].mutex ) ||
        pthread_cond_destroy( &SPU_Control_State[id].wake ) )
      ERROR(( "Unable to destroy spu control thread resources!" ));

  /* Finish up */

  spu.n_pipeline = Id = 0;
  Busy = 0;
}

/****************************************************************************
 *
 * Function: parallel_execute(pipeline,params,job,n_job,flag)
 *
 * Sends a spu to work on a given pipeline.
 *
 * Arguments: pipeline (spe_program_handle_t *) - Pointer to pipeline to be
 *              executed
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
 *         - spu dispatcher has not been booted or the caller is not
 *           the ppu host thread.
 *         - Parallel_execute detected spu_halt terminating spu control
 *           threads.  This should never happen in this implementation.
 *
 * Dependencies: pthread_mutex_unlock, pthread_self, pthread_mutex_trylock,
 *               pthread_equal, pthread_mutex_lock, pthread_cond_signal
 *
 * Globals read/altered: spu.n_pipeline (r), Pipeline (rw)
 *
 * Notes:
 * - Only the host may call this.
 * - If all spu control threads are busy, the function will wait until
 *   a spu control thread becomes available to execute the pipeline.
 *   There is a possibility of a deadlock here if the none of the other
 *   pipelines terminate.  Thus, this library assumes pipelines _always_
 *   eventually terminate.
 * - It is remotely possible for parallel_execute to reacquire the
 *   lock on a spu control thread signalled to execute a pipeline before
 *   the spu control thread wakes up. However, this is not a problem
 *   because the Pipeline[id].state will prevent the pipeline from being
 *   signaled a second time (see innermost switch).
 *
 ***************************************************************************/

static void
parallel_execute( spe_program_handle_t * pipeline,
                  void *args,
                  int size,
                  int job,
                  int n_job,
                  volatile int *flag ) {
  int id;

  /* Determine that the spu dispatcher is initialized and that
     caller is the host spu. */

  if( spu.n_pipeline==0 )
    ERROR(( "Boot the spu dispatcher first!" ));

  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may call parallel_execute" ));

  /* Loop until we hand off the task to a spu control thread (or we
     detect spu_halt is running) */

  for(;;) {

    /* Get an id of a spu control thread to query. */

    id = Id; if( (++Id) >= spu.n_pipeline ) Id = 0;

    if( !pthread_mutex_trylock( &SPU_Control_State[id].mutex ) ) {
      
      /* If the spu control thread isn't executing a pipeline, see if
         the task can be handed off to this spu control thread.  Note:
         host ppu thread now has control over the spu control thread's
         state. */
      
      switch( SPU_Control_State[id].state ) {
        
        /* Note: all cases relinquish control over the spu control
           thread's state */
        
      case SPU_CONTROL_THREAD_SLEEP:
        
        /* The spu control thread is available - assign the task, wake
           the spu control thread and return */
        
        SPU_Control_State[id].state    = SPU_CONTROL_THREAD_EXECUTE;
        SPU_Control_State[id].pipeline = pipeline;
        SPU_Control_State[id].args     = args;
        SPU_Control_State[id].job      = job;
        SPU_Control_State[id].n_job    = n_job;
        SPU_Control_State[id].flag     = flag;
        pthread_cond_signal( &SPU_Control_State[id].wake );
        pthread_mutex_unlock( &SPU_Control_State[id].mutex );
        return;
        
      case SPU_CONTROL_THREAD_TERMINATE:
        
        /* The query spu control thread was terminated.  Something is
           amiss! */
        
        pthread_mutex_unlock( &SPU_Control_State[id].mutex );
        ERROR(( "parallel_execute called while spu_halt is running - "
                "Should never happen in this restricted implementation." ));
        
      default:
        
        /* Spu control thread isn't ready to accept new tasks (it is
           about to wake-up to execute an already assigned task) */
        
        pthread_mutex_unlock( &SPU_Control_State[id].mutex );
        
      } /* switch */
    } /* if spu control thread might be available */
  } /* for */  
  /* Never get here */
}

/****************************************************************************
 *
 * Function: spu_control_thread(this) - internal use only
 *
 * The spu control thread side of the parallel_execute
 *
 * Arguments: this (void *) - pointer to this spu control thread's state
 *
 * Returns: (void *) NULL always.
 *
 * Dependencies: pthread_mutex_lock, pthread_mutex_unlock, pthread_cond_wait
 *
 * Globals read/altered: SPU_Control_State (rw)
 *
 * Notes:
 * - When spu_control_thread is not sleeping (SPU_CONTROL_THREAD_SLEEP
 *   state), it has a mutex lock over over SPU_Control_State[id].state
 *   to insure consistency in the loop flow control
 *
 ***************************************************************************/

void *
spu_control_thread( void *_this_control_thread ) {
  spu_control_state_t *this_control_thread = _this_control_thread;
  unsigned int entry;
  uint64_t env;

  /* Spu control thread state is SPU_CONTROL_THREAD_ACK and the spu
     control thread mutex is unlocked when entering.  Since
     pthread_cond_wait unlockes the spu control thread mutex when the
     spu control thread goes to sleep and the SPU_CONTROL_THREAD_ACK
     case writes the spu control thread state, the mutex needs to be
     locked first. */
     
  pthread_mutex_lock( &this_control_thread->mutex );

  /* Loop while the spu control thread is still executing */

  for(;;) {

    switch( this_control_thread->state ) {

    case SPU_CONTROL_THREAD_TERMINATE:

      /* Terminate the spu control thread. Unlock the mutex first (as
         pthread_cond_wait locks it when the spu control thread wakes
         up). */

      pthread_mutex_unlock( &this_control_thread->mutex );
      return NULL;

    case SPU_CONTROL_THREAD_EXECUTE:

      /* Note: the spu control thread mutex is locked while the spu
         control thread is executing a task. */

      /* Load the program into the context, execute the given task and
         set the completion flag if necessary. */

      if( spe_program_load( this_control_thread->context,
                            this_control_thread->pipeline ) )
        ERROR(( "spe_program_load failed!" ));

      entry = SPE_DEFAULT_ENTRY;
      env   = (((uint64_t)(this_control_thread->n_job))<<32) |
               ((uint64_t)(this_control_thread->job  ));
      if( spe_context_run( this_control_thread->context,
                           &entry,
                           0,
                           this_control_thread->args,
                           (void *)env,
                           NULL ) )
        ERROR(( "spe_context_run failed!" ));

      if( this_control_thread->flag ) *this_control_thread->flag = 1;

      /* Pass through into the next case */

    case SPU_CONTROL_THREAD_ACK:
    case SPU_CONTROL_THREAD_SLEEP:
    default:

      /* Go to sleep. Note: pthread_cond_wait unlocks the spu control
         thread mutex while the spu control thread is sleeping and
         locks it when the spu control thread wakes up */

      this_control_thread->state = SPU_CONTROL_THREAD_SLEEP;
      pthread_cond_wait( &this_control_thread->wake,
                         &this_control_thread->mutex );
      break;

    } /* switch */

  } /* for */

  return NULL; /* Never get here - Avoid compiler warning */
}

static void
spu_dispatch( spe_program_handle_t * pipeline,
              void * args,
              int size ) {
  int rank;

  if( spu.n_pipeline==0 ) ERROR(( "Boot the spu dispatcher first!" ));

  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may call spu_dispatch" ));

  if( (uint64_t)args & 15 ) ERROR(( "args must be 16-byte aligned" ));
  if( size & 15 )           ERROR(( "size must be a multiple of 16" ));

  if( Busy ) ERROR(( "Pipelines are busy!" ));

  for( rank=0; rank<spu.n_pipeline; rank++ ) {
    Done[rank] = 0;
    parallel_execute( pipeline,
                      ((char *)args) + rank*size,
                      size,
                      rank,
                      spu.n_pipeline,
                      &Done[rank] );
  }

  Busy = 1;
}
                 
static void
spu_wait( void ) {
  int rank;

  if( spu.n_pipeline==0 ) ERROR(( "Boot the spu dispatcher first!" ));

  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may call spu_wait" ));

  if( !Busy ) ERROR(( "Pipelines are not busy!" ));

  rank = 0;
  while( rank<spu.n_pipeline ) {
    if( Done[rank] ) rank++;
    else             nanodelay(5);
  }

  Busy = 0;
}

typedef void (*dispatcher_func_t)( pipeline_func_t, void *, int );

pipeline_dispatcher_t spu = {
  0,            /* n_pipeline */
  spu_boot,     /* boot */
  spu_halt,     /* halt */
  (dispatcher_func_t)spu_dispatch, /* dispatch */
  spu_wait      /* wait */
};

