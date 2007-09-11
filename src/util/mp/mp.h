#ifndef _mp_h_
#define _mp_h_

#include <mp.hxx>

// Note: mp module assumes homogeneous cluster (in particular,
// bit-for-bit compatible data layouts between nodes

BEGIN_C_DECLS

static inline void mp_init(int argc, char ** argv) {
	mp_init_cxx(argc, argv);
} // mp_finalize

static inline void mp_finalize( mp_handle h ) {
	mp_finalize_cxx(h);
} // mp_finalize

static inline mp_handle new_mp(void) {
	return new_mp_cxx();
} // new_mp

static inline void delete_mp( mp_handle *h ) {
	delete_mp_cxx(h);
} // delete_mp

static inline int mp_rank( mp_handle h ) {
	return mp_rank_cxx(h);
} // mp_rank

static inline int mp_nproc( mp_handle h ) {
	return mp_nproc_cxx(h);
} // mp_nproc

static inline void * ALIGNED(16) mp_recv_buffer( int rbuf, mp_handle h ) {
	return mp_recv_buffer_cxx(rbuf, h);
} // mp_recv_buffer

static inline void * ALIGNED(16) mp_send_buffer( int sbuf, mp_handle h ) {
	return mp_send_buffer_cxx(sbuf, h);
} // mp_send_buffer

static inline void mp_abort( int reason, mp_handle h ) {
	mp_abort_cxx(reason, h);
} // mp_abort

static inline void mp_barrier( mp_handle h ) {
	mp_barrier_cxx(h);
} // mp_barrier

static inline double mp_elapsed( mp_handle h ) {
	return mp_elapsed_cxx(h);
} // mp_elapsed

static inline double mp_time00( mp_handle h ) {
	return mp_time00_cxx(h);
} // mp_time00

static inline double mp_wtime(void) {
	return mp_wtime_cxx();
} // mp_wtime

static inline error_code mp_size_recv_buffer( int rbuf, int size,
	mp_handle h ) {
	return mp_size_recv_buffer_cxx(rbuf, size, h);
} // mp_size_recv_buffer

static inline error_code mp_size_send_buffer( int sbuf, int size,
	mp_handle h ) {
	return mp_size_send_buffer_cxx(sbuf, size, h);
} // mp_size_send_buffer

static inline error_code mp_begin_recv( int rbuf, int size,
	int sender, int tag, mp_handle h ) {
	return mp_begin_recv_cxx(rbuf, size, sender, tag, h);
} // mp_begin_recv

static inline error_code mp_begin_send( int sbuf, int size, int receiver,
	int tag, mp_handle h ) {
	return mp_begin_send_cxx(sbuf, size, receiver, tag, h);
} // mp_begin_send

static inline error_code mp_end_recv( int rbuf, mp_handle h ) {
	return mp_end_recv_cxx(rbuf, h);
} // mp_end_recv

static inline error_code mp_end_send( int sbuf, mp_handle h ) {
	return mp_end_send_cxx(sbuf, h);
} // mp_end_send

static inline error_code mp_allsum_d( double *local, double *global, int n,
	mp_handle h ) {
	return mp_allsum_d_cxx(local, global, n, h);
} // mp_allsum_d

static inline error_code mp_allsum_i( int *local, int *global,
	int n, mp_handle h ) {
	return mp_allsum_i_cxx(local, global, n, h);
} // mp_allsum_i

static inline error_code mp_allgather_i( int *sbuf, int *rbuf,
	int n, mp_handle h ) {
	return mp_allgather_i_cxx(sbuf, rbuf, n, h);
} // mp_allgather_i

static inline error_code mp_allgather_i64( int64_t *sbuf,
	int64_t *rbuf, int n, mp_handle h ) {
	return mp_allgather_i64_cxx(sbuf, rbuf, n, h);
} // mp_allgather_i64

static inline error_code mp_send_i( int *buf, int n, int dst, mp_handle h ) {
	return mp_send_i_cxx(buf, n, dst, h);
} // mp_send_i

static inline error_code mp_recv_i( int *buf, int n, int src, mp_handle h ) {
	return mp_recv_i_cxx(buf, n, src, h);
} // mp_recv_i

END_C_DECLS

#endif
