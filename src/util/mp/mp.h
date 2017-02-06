// FIXME: THIS API SHOULD NOT BE HARDWIRED TO USE ONLY THE WORLD
// COMM. FIXME: THIS API NEEDS A SERIOUS REVAMP (BUT AT LEAST IT IS
// LESS A HOUSE OF SHAME THAN PREVIOUSLY).

#ifndef mp_h
#define mp_h

#include "../util_base.h"

/* Opaque handle to the message passing buffers */

struct mp;
typedef struct mp mp_t;

/* Define a "turnstile".  At most up to n_turnstile processes can be
   in the turnstile at any given time.  Use this to implement
   critical sections and do other tricks liking limiting the number
   of simultaneous I/O operatorions on large jobs.  These macros use
   blocking send/receives to serialize writes.

   For example, to set up a turnstile that allows at most N
   simultaneous writes:

   BEGIN_TURNSTILE( N ) {
     ... do write ...
   } END_TURNSTILE

   BEGIN_TURNSTILE(1) (i.e., one turnstile) effectively serializes the
   code.  This construct is robust.  Turnstiles should not be nested.
   Code in turnstiles should not attempt to communicate with other
   processes.

   If everything were perfectly synchronous, then, when
   using a 10 turnstiles, processes 0:9 would enter the turnstile,
   followed by 10:19, followed by 20:29, ... */

#define BEGIN_TURNSTILE(n_turnstile) do {               \
   int _n_turnstile = (n_turnstile), _baton;            \
   if( world_rank>=_n_turnstile )                       \
      mp_recv_i( &_baton, 1, world_rank-_n_turnstile ); \
   do

#define END_TURNSTILE while(0);                         \
   if( world_rank+_n_turnstile < world_size )           \
     mp_send_i( &_baton, 1, world_rank+_n_turnstile );  \
 } while(0)

BEGIN_C_DECLS

void
boot_mp( int * pargc,
         char *** pargv );

void
halt_mp( void );

void
mp_abort( int reason );

/* Collective commucations */

void
mp_barrier( void );

void
mp_allsum_d( double * local,
             double * global,
             int n );

void
mp_allsum_i( int * local,
             int * global,
             int n );

void
mp_allgather_i( int * sbuf,
                int * rbuf,
                int n );

void
mp_allgather_i64( int64_t * sbuf,
                  int64_t * rbuf,
                  int n );

// FIXME: THIS API SHOULD TAKE THE ROOT NODE
void
mp_gather_uc( unsigned char * sbuf,
              unsigned char * rbuf,
              int n );

/* Turnstile communication primitives */
// FIXME: MESSAGE TAGGING ISSUES?

void
mp_send_i( int * buf,
           int n,
           int dst );

void
mp_recv_i( int * buf,
           int n,
           int src );

/* Buffered non-blocking point-to-point communications */

mp_t *
new_mp( int n_port );

void
delete_mp( mp_t * mp );

void * ALIGNED(128)
mp_recv_buffer( mp_t * mp,
                int port );

void * ALIGNED(128)
mp_send_buffer( mp_t * mp,
                int port );

void
mp_size_recv_buffer( mp_t * mp,
                     int port,
                     int size );

void
mp_size_send_buffer( mp_t * mp,
                     int port,
                     int size );

// FIXME: MP REALLY SHOULD HANDLE THE MESSAGE TAGGING
void
mp_begin_recv( mp_t * mp,
               int port,
               int sz,
               int src,
               int tag );

void
mp_begin_send( mp_t * mp,
               int port,
               int sz,
               int dst,
               int tag );

void
mp_end_recv( mp_t * mp,
             int rbuf );

void
mp_end_send( mp_t * mp,
             int sbuf );

END_C_DECLS

#endif /* mp_h */
