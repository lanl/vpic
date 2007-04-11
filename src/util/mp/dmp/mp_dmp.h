#ifndef _mp_dmp_h_
#define _mp_dmp_h_

#include <util_base.h>
#include <mp_t.h>

/* Note: mp module assumes homogeneous cluster (in particular,
   bit-for-bit compatible data layouts between nodes */

BEGIN_C_DECLS

mp_handle new_mp_dmp(void);
void delete_mp_dmp( mp_handle *h );

int mp_rank_dmp( mp_handle h );
int mp_nproc_dmp( mp_handle h );

void * ALIGNED mp_recv_buffer_dmp( int rbuf, mp_handle h );
void * ALIGNED mp_send_buffer_dmp( int sbuf, mp_handle h );

double mp_elapsed_dmp( mp_handle h );
double mp_time00_dmp( mp_handle h );

void mp_abort_dmp( int reason, mp_handle h );
void mp_finalize_dmp( mp_handle h ); 
void mp_barrier_dmp( mp_handle h );

error_code mp_size_recv_buffer_dmp( int rbuf, int size, mp_handle h );
error_code mp_size_send_buffer_dmp( int sbuf, int size, mp_handle h );

error_code mp_begin_recv_dmp( int rbuf, int size,
                                 int sender, int tag, mp_handle h );
error_code mp_begin_send_dmp( int sbuf, int size,
                                 int receiver, int tag, mp_handle h );

error_code mp_end_recv_dmp( int rbuf, mp_handle h );
error_code mp_end_send_dmp( int sbuf, mp_handle h );

error_code mp_allsum_d_dmp( double *local, double *global, int n,
                               mp_handle h );
error_code mp_allsum_i_dmp( int *local, int *global, int n,
                               mp_handle h );
error_code mp_allgather_i_dmp( int *sbuf, int *rbuf, int n, mp_handle h );
error_code mp_allgather_i64_dmp( int64_t *sbuf, int64_t *rbuf,
	int n, mp_handle h );

error_code mp_send_i_dmp( int *buf, int n, int dst, mp_handle h );
error_code mp_recv_i_dmp( int *buf, int n, int src, mp_handle h );

END_C_DECLS

#endif // _mp_dmp_h_
