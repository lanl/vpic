/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "mp.h"
#include "MPWrapper.hxx"

void boot_mp( int * pargc, char *** pargv ) {
  MPWrapper::instance().boot_mp( pargc, pargv );
}

void halt_mp( void ) { MPWrapper::instance().halt_mp(); }

void mp_abort( int reason ) { MPWrapper::instance().mp_abort( reason ); }

void mp_barrier( void ) { MPWrapper::instance().mp_barrier(); }

void mp_allsum_d( double *local, double *global, int n ) {
  return MPWrapper::instance().mp_allsum_d( local, global, n );
}

void mp_allsum_i( int *local, int *global, int n ) {
  return MPWrapper::instance().mp_allsum_i( local, global, n );
}

void mp_allgather_i( int *sbuf, int *rbuf, int n ) {
  return MPWrapper::instance().mp_allgather_i( sbuf, rbuf, n );
}

void mp_allgather_i64( int64_t *sbuf, int64_t *rbuf, int n ) {
  return MPWrapper::instance().mp_allgather_i64( sbuf, rbuf, n );
}

void mp_gather_uc( unsigned char * sbuf, unsigned char * rbuf, int n ) {
  return MPWrapper::instance().mp_gather_uc( sbuf, rbuf, n );
}

void mp_send_i( int *buf, int n, int dst ) {
  return MPWrapper::instance().mp_send_i( buf, n, dst );
}

void mp_recv_i( int *buf, int n, int src ) {
  return MPWrapper::instance().mp_recv_i( buf, n, src );
}

mp_t * new_mp( int n_port ) { return MPWrapper::instance().new_mp( n_port ); }

void delete_mp( mp_t * mp ) { MPWrapper::instance().delete_mp( mp ); }

void * ALIGNED(16) mp_recv_buffer( mp_t * mp, int tag ) {
  return MPWrapper::instance().mp_recv_buffer( mp, tag );
}

void * ALIGNED(16) mp_send_buffer( mp_t * mp, int tag ) {
  return MPWrapper::instance().mp_send_buffer( mp, tag );
}

void mp_size_recv_buffer( mp_t * mp, int tag, int size ) {
  MPWrapper::instance().mp_size_recv_buffer( mp, tag, size );
}

void mp_size_send_buffer( mp_t * mp, int tag, int size ) {
  MPWrapper::instance().mp_size_send_buffer( mp, tag, size );
}

void mp_begin_recv( mp_t * mp, int rbuf, int size, int sender, int tag ) {
  MPWrapper::instance().mp_begin_recv( mp, rbuf, size, sender, tag );
}

void mp_begin_send( mp_t * mp, int sbuf, int size, int receiver, int tag ) {
  MPWrapper::instance().mp_begin_send( mp, sbuf, size, receiver, tag );
}

void mp_end_recv( mp_t * mp, int rbuf ) {
  MPWrapper::instance().mp_end_recv( mp, rbuf );
}

void mp_end_send( mp_t * mp, int sbuf ) {
  MPWrapper::instance().mp_end_send( mp, sbuf );
}

