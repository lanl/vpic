/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "MPWrapper.hxx"
#include "mp.hxx"

// FIXME: MP_HANDLE IS DUBIOUS WAY TO DO AN OPAQUE DATA TYPE

void
mp_init_cxx( int argc,
             char ** argv ) {
  MPWrapper::instance().mp_init( argc, argv );
} // mp_init

// Adding for clean termination
void
mp_finalize_cxx( mp_handle h ) {
  MPWrapper::instance().mp_finalize( h );
} // mp_finalize

mp_handle
new_mp_cxx( void ) {
  return MPWrapper::instance().new_mp();
} // new_mp

void
delete_mp_cxx( mp_handle *h ) {
  MPWrapper::instance().delete_mp( h );
} // delete_mp

int
mp_rank_cxx( mp_handle h ) {
  return MPWrapper::instance().mp_rank( h );
} // mp_rank

int
mp_nproc_cxx( mp_handle h ) {
  return MPWrapper::instance().mp_nproc( h );
} // mp_nproc

void * ALIGNED(16)
mp_recv_buffer_cxx( int tag,
                    mp_handle h ) {
  return MPWrapper::instance().mp_recv_buffer( tag, h );
} // mp_recv_buffer

void * ALIGNED(16)
mp_send_buffer_cxx( int tag,
                    mp_handle h ) {
  return MPWrapper::instance().mp_send_buffer( tag, h );
} // mp_send_buffer

void
mp_abort_cxx( int reason,
              mp_handle h ) {  
  MPWrapper::instance().mp_abort( reason, h );
} // mp_abort

void
mp_barrier_cxx( mp_handle h ) {
  MPWrapper::instance().mp_barrier( h );
} // mp_barrier

// Returns the time elapsed since the communicator was created.  Every
// rank gets the same value
double
mp_elapsed_cxx( mp_handle h ) {
  return MPWrapper::instance().mp_elapsed( h );
} // mp_elapsed

// Stop watch. Different ranks may get different values.  First call
// to stop watch returns a measure of the overhead.
double
mp_time00_cxx( mp_handle h ) {
  return MPWrapper::instance().mp_time00( h );
} // mp_time00

double
mp_wtime_cxx(void) {
  return MPWrapper::instance().mp_wtime();
} // mp_wtime

void
mp_size_recv_buffer_cxx( int tag,
                         int size,
                         mp_handle h ) {
  return MPWrapper::instance().mp_size_recv_buffer( tag, size, h );
} // mp_size_recv_buffer

void
mp_size_send_buffer_cxx( int tag,
                         int size,
                         mp_handle h ) {
  return MPWrapper::instance().mp_size_send_buffer( tag, size, h );
} // mp_size_send_buffer

void
mp_begin_recv_cxx( int rbuf,
                   int size,
                   int sender,
                   int tag,
                   mp_handle h ) {
  return MPWrapper::instance().mp_begin_recv( rbuf, size, sender, tag, h );
} // mp_begin_recv

void
mp_begin_send_cxx( int sbuf,
                   int size,
                   int receiver,
                   int tag,
                   mp_handle h ) {
  return MPWrapper::instance().mp_begin_send( sbuf, size,
                                              receiver, tag, h );
} // mp_begin_send

void
mp_end_recv_cxx( int rbuf,
                 mp_handle h ) {
  return MPWrapper::instance().mp_end_recv( rbuf, h );
} // mp_end_recv

void
mp_end_send_cxx( int sbuf,
                 mp_handle h ) {
  return MPWrapper::instance().mp_end_send( sbuf, h );
} // mp_end_send

void
mp_allsum_d_cxx( double *local,
                 double *global,
                 int n,
                 mp_handle h ) {
  return MPWrapper::instance().mp_allsum_d( local, global, n, h );
} // mp_allsum_d

void
mp_allsum_i_cxx( int *local,
                 int *global,
                 int n,
                 mp_handle h ) {
  return MPWrapper::instance().mp_allsum_i( local, global, n, h );
} // mp_allsum_i

void
mp_allgather_i_cxx( int *sbuf,
                    int *rbuf,
                    int n,
                    mp_handle h ) {
  return MPWrapper::instance().mp_allgather_i( sbuf, rbuf, n, h );
} // mp_allgather_i

void
mp_allgather_i64_cxx( int64_t *sbuf,
                      int64_t *rbuf,
                      int n,
                      mp_handle h ) {
  return MPWrapper::instance().mp_allgather_i64( sbuf, rbuf, n, h );
} // mp_allgather_i64

void
mp_gather_uc_cxx( unsigned char * sbuf,
                  unsigned char * rbuf,
                  int n,
                  mp_handle h ) {
  return MPWrapper::instance().mp_gather_uc( sbuf, rbuf, n, h );
} // mp_gather_uc_cxx

// We need blocking send/receive to implement turnstiles.

void
mp_send_i_cxx( int *buf,
               int n,
               int dst,
               mp_handle h ) {
  return MPWrapper::instance().mp_send_i( buf, n, dst, h );
} // mp_send_i

void
mp_recv_i_cxx( int *buf,
               int n,
               int src,
               mp_handle h ) {
  return MPWrapper::instance().mp_recv_i( buf, n, src, h );
} // mp_recv_i
