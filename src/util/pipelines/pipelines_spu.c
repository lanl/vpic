#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)

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
// Modified extensively for VPIC-3P SPU dispatch by K. Bowers 4/20/2007

// TO DO:
// (?) Signal blocking in spu_control_threads
// (?) Timeouts in spu_halt, spu_boot (spin wait)

#include "pipelines.h"
#include <pthread.h>
#include <libspe2.h>

/* Must be same as SPU_COMPLETE in root_segment.c; Neil Peart rules!  */
#define SPU_COMPLETE ((uint32_t)2112)

static void *
spu_control_thread( void *_id );

extern spe_program_handle_t root_segment;

typedef struct spu_control_state {
  pthread_t handle;
  spe_context_ptr_t context;
  spe_program_handle_t * pipeline;
  void *args;
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
spu_boot( int num_pipe,
          int dispatch_to_host ) {
  size_t i;

  // Check if arguments are valid and dispatcher isn't already initialized

  if( spu.n_pipeline != 0 ) ERROR(( "Halt the spu dispatcher first!" ));
  if( num_pipe < 1 || num_pipe > MAX_PIPELINE )
    ERROR(( "Invalid number of pipelines requested" ));
  if( dispatch_to_host )
    ERROR(( "SPU pipelines cannot be dispatched to the host" ));

  // Initialize some global variables. Note: spu.n_pipeline = 0 here

  Host = pthread_self();
  Id = 0;

  MESSAGE(("In spu_boot with num_pipe = %d", num_pipe));

  // Initialize all the pipelines

  for( i=0; i<num_pipe; i++ ) {

	// Create SPE context
    if( ( SPU_Control_State[i].context=spe_context_create( 0, NULL ) )==NULL )
      ERROR(( "spe_context_create failed" ));

    // Set the spe pipeline to use the new root_segment generic
	// event loop
	SPU_Control_State[i].pipeline = &root_segment;

	// Create pthread to run SPE thread non-blocking
    if( pthread_create( &SPU_Control_State[i].handle,
                        NULL,
                        spu_control_thread,
                        &SPU_Control_State[i] ) )
      ERROR(( "pthread_create failed" ));

    // Nothing to do. After pthread returns we're done.
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

  // Make sure the host is the caller and there are spu control
  // threads to terminate

  if( !spu.n_pipeline ) ERROR(( "Boot the spu dispatcher first!" ));
  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may halt the spu dispatcher!" ));

  // Terminate the spu control threads

  uint32_t abort = 0;
  for( id=0; id<spu.n_pipeline; id++ ) {
    
    spe_in_mbox_write( SPU_Control_State[id].context, &abort,
                       1, SPE_MBOX_ANY_NONBLOCKING );

    if( pthread_join( SPU_Control_State[id].handle, NULL ) )
      ERROR(( "Unable to terminate spu control thread!" ));

  }

  // Free resources associated with dispatcher as all spu control
  // threads are now dead.  Note: This must be done in a separate loop
  // because non-terminated spu control threads calling
  // parallel_execute may try to access mutexes for destroyed spu
  // control threads and Id and with unknown results if the mutexes
  // were destroyed before all the spu control threads were
  // terminated.  This note does not apply for this restricted
  // dispatcher.

  for( id=0; id<spu.n_pipeline; id++ )
    if( spe_context_destroy( SPU_Control_State[id].context ) )
      ERROR(( "Unable to destroy spu control thread resources!" ));

  // Finish up

  spu.n_pipeline = Id = 0;
  Busy = 0;
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

static void *
spu_control_thread( void *_this_control_thread ) {
  spu_control_state_t *this_control_thread = _this_control_thread;
  unsigned int entry;

  MESSAGE(("Loading spe program"));
  if( spe_program_load( this_control_thread->context,
                        this_control_thread->pipeline ) ) {
    ERROR(( "spe_program_load failed!" ));
  } // if

  entry = SPE_DEFAULT_ENTRY;

  MESSAGE(("Running spe context"));
  if( spe_context_run( this_control_thread->context,
                       &entry,
                       0,
                       this_control_thread->args,
                       NULL,
                       NULL ) ) {
    ERROR(( "spe_context_run failed!" ));
  } // if

  pthread_exit(NULL);

  return NULL; // Never get here - Avoid compiler warning
}

static void
spu_dispatch( pipeline_func_t func,
              void * args,
              int sz_args ) {
  uint32_t rank;
  uint32_t data_hi, data_lo;

  if( spu.n_pipeline==0 ) ERROR(( "Boot the spu dispatcher first!" ));

  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may call spu_dispatch" ));

  if( (uint64_t)args & 15 ) ERROR(( "args must be 16-byte aligned" ));
  if( sz_args & 15 )        ERROR(( "size must be a multiple of 16" ));

  if( Busy ) ERROR(( "Pipelines are busy!" ));

  // Put address into a form that we can send over to the SPE threads
  uint32_t faddr = ((uint64_t)(func)) & 0xffffffff;

  for( rank=0; rank<spu.n_pipeline; rank++ ) {
    Done[rank] = 0;

    // Write pipeline to SPEs
    spe_in_mbox_write( SPU_Control_State[rank].context, &faddr,
                       1, SPE_MBOX_ANY_NONBLOCKING );

    // get low and hi bits to send context struct to SPE thread
    data_hi = ((uint64_t)(((char *)args) + rank*sz_args)) >> 32;
    data_lo = ((uint64_t)(((char *)args) + rank*sz_args)) & 0xffffffff;

    // Write context struct address to SPE thread
    spe_in_mbox_write( SPU_Control_State[rank].context, &data_lo,
                       1, SPE_MBOX_ANY_NONBLOCKING );
    spe_in_mbox_write( SPU_Control_State[rank].context, &data_hi,
                       1, SPE_MBOX_ANY_NONBLOCKING );

    // Write pipeline_rank and n_pipeline to SPE thread
    spe_in_mbox_write( SPU_Control_State[rank].context, &rank,
                       1, SPE_MBOX_ANY_NONBLOCKING );
	// This is klugey (Not sure why n_pipeline is not unsigned in rest of code)
	uint32_t n_pipeline = (uint32_t)(spu.n_pipeline);
    spe_in_mbox_write( SPU_Control_State[rank].context, &n_pipeline,
                       1, SPE_MBOX_ANY_NONBLOCKING );
  }

  Busy = 1;
}
                 
static void
spu_wait( void ) {
  int rank;
  uint32_t msg;

  if( spu.n_pipeline==0 ) ERROR(( "Boot the spu dispatcher first!" ));

  if( !pthread_equal( Host, pthread_self() ) )
    ERROR(( "Only the host may call spu_wait" ));

  if( !Busy ) ERROR(( "Pipelines are not busy!" ));

  for(rank=0; rank<spu.n_pipeline; rank++) {
    if( !Done[rank] ) {
      while( spe_out_mbox_status(SPU_Control_State[rank].context)<=0 );
      spe_out_mbox_read(SPU_Control_State[rank].context, &msg, 1);
      if( msg!=SPU_COMPLETE ) ERROR(("Bad SPE synchronization message"));
      Done[rank] = 1;
    } // if
  } // for

  Busy = 0;
}

pipeline_dispatcher_t spu = {
  0,            // n_pipeline
  spu_boot,     // boot
  spu_halt,     // halt
  spu_dispatch, // dispatch
  spu_wait      // wait
};

#endif
