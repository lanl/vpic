#ifndef DMPPolicy_hxx
#define DMPPolicy_hxx

#include "mp_t.h"
#include "mp_dmp.h"

struct DMPPolicy {

  inline void
  mp_init( int argc,
           char ** argv ) {
    return mp_init_dmp( argc, argv );
  } // mp_init
  
  inline void
  mp_finalize( mp_handle h ) {
    return mp_finalize_dmp( h );
  } // mp_finalize
  
  inline mp_handle
  new_mp( void ) {
    return new_mp_dmp();
  } // new_mp
  
  inline void
  delete_mp( mp_handle *h ) {
    return delete_mp_dmp( h );
  } // delete_mp
  
  inline int
  mp_rank( mp_handle h ) {
    return mp_rank_dmp( h );
  } // mp_rank
  
  inline int
  mp_nproc( mp_handle h ) {
    return mp_nproc_dmp( h );
  } // mp_nproc
  
  inline void * ALIGNED(16)
  mp_recv_buffer( int rbuf,
                  mp_handle h ) {
    return mp_recv_buffer_dmp( rbuf, h );
  } // mp_recv_buffer
  
  inline void * ALIGNED(16)
  mp_send_buffer( int sbuf,
                  mp_handle h ) {
    return mp_send_buffer_dmp( sbuf, h );
  } // mp_send_buffer
  
  inline double
  mp_elapsed( mp_handle h ) {
    return mp_elapsed_dmp( h );
  } // mp_elapsed
  
  inline double
  mp_time00( mp_handle h ) {
    return mp_time00_dmp( h );
  } // mp_time00
  
  inline double
  mp_wtime() {
    return mp_wtime_dmp();
  } // mp_wtime
  
  inline void
  mp_abort( int reason,
            mp_handle h ) {
    return mp_abort_dmp( reason, h );
  } // mp_abort
  
  inline void
  mp_barrier( mp_handle h ) {
    return mp_barrier_dmp( h );
  } // mp_barrier
  
  inline void
  mp_size_recv_buffer( int rbuf,
                       int size,
                       mp_handle h ) {
    return mp_size_recv_buffer_dmp( rbuf, size, h );
  } // mp_size_recv_buffer
  
  inline void
  mp_size_send_buffer( int sbuf,
                       int size,
                       mp_handle h ) {
    return mp_size_send_buffer_dmp( sbuf, size, h );
  } // mp_size_send_buffer
  
  inline void
  mp_begin_recv( int rbuf,
                 int size,
                 int sender,
                 int tag,
                 mp_handle h ) {
    return mp_begin_recv_dmp( rbuf, size, sender, tag, h );
  } // mp_begin_recv
  
  inline void
  mp_begin_send( int sbuf,
                 int size,
                 int receiver,
                 int tag,
                 mp_handle h ) {
    return mp_begin_send_dmp( sbuf, size, receiver, tag, h );
  } // mp_begin_send
  
  inline void
  mp_end_recv( int rbuf,
               mp_handle h ) {
    return mp_end_recv_dmp( rbuf, h );
  } // mp_end_recv
  
  inline void
  mp_end_send( int sbuf,
               mp_handle h ) {
    return mp_end_send_dmp( sbuf, h );
  } // mp_end_send
  
  inline void
  mp_allsum_d( double *local,
               double *global,
               int n,
               mp_handle h ) {
    return mp_allsum_d_dmp( local, global, n, h );
  } // mp_allsum_d
  
  inline void
  mp_allsum_i( int *local,
               int *global,
               int n,
               mp_handle h ) {
    return mp_allsum_i_dmp( local, global, n, h );
  } // mp_allsum_i
  
  inline void
  mp_allgather_i( int *sbuf,
                  int *rbuf,
                  int n,
                  mp_handle h ) {
    return mp_allgather_i_dmp( sbuf, rbuf, n, h );
  } // mp_allgather_i
  
  inline void
  mp_allgather_i64( int64_t *sbuf,
                    int64_t *rbuf,
                    int n,
                    mp_handle h ) {
    return mp_allgather_i64_dmp( sbuf, rbuf, n, h );
  } // mp_allgather_i64
  
  inline void
  mp_send_i( int *buf,
             int n,
             int dst,
             mp_handle h ) {
    return mp_send_i_dmp( buf, n, dst, h );
  } // mp_send_i
  
  inline void
  mp_recv_i( int *buf,
             int n,
             int src,
             mp_handle h ) {
    return mp_recv_i_dmp( buf, n, src, h );
  } // mp_recv_i

}; // struct DMPPolicy

#endif // DMPPolicy_hxx
