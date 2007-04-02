#ifndef _mp_h_
#define _mp_h_

#include <common.h>

/* Note: mp module assumes homogeneous cluster (in particular,
   bit-for-bit compatible data layouts between nodes */

/* Opaque data type so users do not have to see message passing internals */
typedef void *mp_handle;

BEGIN_C_DECLS

mp_handle new_mp(void);
void delete_mp( mp_handle *h );

int mp_rank( mp_handle h );
int mp_nproc( mp_handle h );

void * ALIGNED mp_recv_buffer( int rbuf, mp_handle h );
void * ALIGNED mp_send_buffer( int sbuf, mp_handle h );

double mp_elapsed( mp_handle h );
double mp_time00( mp_handle h );

void mp_abort( int reason, mp_handle h );
void mp_finalize( mp_handle h ); 
void mp_barrier( mp_handle h );

error_code mp_size_recv_buffer( int rbuf, int size, mp_handle h );
error_code mp_size_send_buffer( int sbuf, int size, mp_handle h );

error_code mp_begin_recv( int rbuf, int size,
                                 int sender, int tag, mp_handle h );
error_code mp_begin_send( int sbuf, int size,
                                 int receiver, int tag, mp_handle h );

error_code mp_end_recv( int rbuf, mp_handle h );
error_code mp_end_send( int sbuf, mp_handle h );

error_code mp_allsum_d( double *local, double *global, int n,
                               mp_handle h );
error_code mp_allsum_i( int *local, int *global, int n,
                               mp_handle h );
error_code mp_allgather_i( int *sbuf, int *rbuf, int n, mp_handle h );
error_code mp_allgather_i64( INT64_TYPE *sbuf, INT64_TYPE *rbuf, int n, mp_handle h );

error_code mp_send_i( int *buf, int n, int dst, mp_handle h );
error_code mp_recv_i( int *buf, int n, int src, mp_handle h );

END_C_DECLS

#endif
