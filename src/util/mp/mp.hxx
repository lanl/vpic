/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */
#ifndef mp_hxx
#define mp_hxx

#include "../util_base.h"
#include "mp_handle.h"

BEGIN_C_DECLS

void
mp_init_cxx( int argc,
             char ** argv );

void
mp_finalize_cxx( mp_handle h );

mp_handle
new_mp_cxx( void );

void
delete_mp_cxx( mp_handle *h );

int
mp_rank_cxx( mp_handle h );

int
mp_nproc_cxx( mp_handle h );

void * ALIGNED(16)
mp_recv_buffer_cxx( int tag,
                    mp_handle h );

void * ALIGNED(16)
mp_send_buffer_cxx( int tag,
                    mp_handle h );

void
mp_abort_cxx( int reason,
              mp_handle h );

void
mp_barrier_cxx( mp_handle h );

// Returns the time elapsed since the communicator was created.  Every
// rank gets the same value
double
mp_elapsed_cxx( mp_handle h );

// Stop watch. Different ranks may get different values.  First call
// to stop watch returns a measure of the overhead
double
mp_time00_cxx( mp_handle h );

double
mp_wtime_cxx( void );

void
mp_size_recv_buffer_cxx( int tag,
                         int size,
                         mp_handle h );

void
mp_size_send_buffer_cxx( int tag,
                         int size,
                         mp_handle h );

void
mp_begin_recv_cxx( int rbuf,
                   int size,
                   int sender,
                   int tag,
                   mp_handle h );

void
mp_begin_send_cxx( int sbuf,
                   int size,
                   int receiver,
                   int tag,
                   mp_handle h );

void
mp_end_recv_cxx( int rbuf,
                 mp_handle h );

void
mp_end_send_cxx( int sbuf,
                 mp_handle h );

void
mp_allsum_d_cxx( double *local,
                 double *global,
                 int n,
                 mp_handle h );

void
mp_allsum_i_cxx( int *local,
                 int *global,
                 int n,
                 mp_handle h );

void
mp_allgather_i_cxx( int *sbuf,
                    int *rbuf,
                    int n,
                    mp_handle h );

void
mp_allgather_i64_cxx( int64_t *sbuf,
                      int64_t *rbuf,
                      int n,
                      mp_handle h );

// We need blocking send/receive to implement turnstiles.

void
mp_send_i_cxx( int *buf,
               int n,
               int dst,
               mp_handle h );

void
mp_recv_i_cxx( int *buf,
               int n,
               int src,
               mp_handle h );

END_C_DECLS

#endif // mp_hxx
