#ifndef AAISPolicy_hxx
#define AAISPolicy_hxx

#include <mp_t.h>
#include <P2PConnection.hxx>

struct AAISPolicy {

	inline void mp_init(int argc, char ** argv) {
		P2PConnection & p2p = P2PConnection::instance();
		p2p.init(argc, argv);
	} // mp_init

	inline void mp_finalize() {
		P2PConnection & p2p = P2PConnection::instance();
		p2p.finalize();
	} // mp_finalize

	inline mp_handle new_mp(void) {
		// get p2p connection
		P2PConnection & p2p = P2PConnection::instance();

		mp_t * mp = static_cast<mp_t *>(malloc(sizeof(mp_t)));
		if(!mp) { return static_cast<mp_handle>(NULL); }

		// get rank and size from point-to-point connection
		mp->rank = p2p.rank();
		mp->nproc = p2p.size();

		// get wtime from point-to-point connection
		mp->elapsed_ref = p2p.wtime();
		mp->time00_ref = 0;
		mp->time00_toggle = 0;

		// NUM_BUF is a global defined in <mp_t.h>
		for(size_t i(0); i<NUM_BUF; i++) {
			mp->rbuf[i] = NULL;
			mp->sbuf[i] = NULL;
			mp->rbuf_size[i] = 0;
			mp->sbuf_size[i] = 0;
			// FIXME: Init rreq and sreq?
			mp->rreq_size[i] = 0;
			mp->sreq_size[i] = 0;
		} // for

		return static_cast<mp_handle>(mp);
	} // new_mp

	inline void delete_mp( mp_handle *h ) {
		mp_t * mp;

		if(!h) { return; }

		mp = static_cast<mp_t *>(*h);

		if(!mp) { return; }

		mp->rank = -1;
		mp->nproc = 0;
		mp->elapsed_ref = 0;
		mp->time00_ref = 0;
		mp->time00_toggle = 0;

		// NUM_BUF is a global defined in <mp_t.h>
		for(size_t i(0); i<NUM_BUF; i++) {
			// allocated with malloc
			if(mp->rbuf[i]) { free_aligned(mp->rbuf[i]); }
			if(mp->sbuf[i]) { free_aligned(mp->sbuf[i]); }

			mp->rbuf[i] = NULL;
			mp->sbuf[i] = NULL;

			mp->rbuf_size[i] = 0;
			mp->sbuf_size[i] = 0;
			// FIXME: Init rreq and sreq?
			mp->rreq_size[i] = 0;
			mp->sreq_size[i] = 0;
		} // for

		free(mp);
		*h = NULL;
	} // delete_mp

	inline int mp_rank( mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		return mp ? mp->rank : -1;	
	} // mp_rank

	inline int mp_nproc( mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		return mp ? mp->nproc : -1;	
	} // mp_nproc

	inline void * ALIGNED mp_recv_buffer( int tag, mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		if(!mp || tag<0 || tag>=NUM_BUF) { return NULL; }
		return mp->rbuf[tag];
	} // mp_recv_buffer

	inline void * ALIGNED mp_send_buffer( int tag, mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		if(!mp || tag<0 || tag>=NUM_BUF) { return NULL; }
		return mp->sbuf[tag];
	} // mp_send_buffer

	inline double mp_elapsed( mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
		double time = p2p.wtime() - mp->elapsed_ref;

		p2p.request(P2PTag::allreduce_max_double, 1);
		p2p.send(&time, 1);
		p2p.recv(&time, 1);

		return time;
	} // mp_elapsed

	inline double mp_time00( mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
		
		if(!mp) { return -1; }

		mp->time00_toggle ^= 1;

		if(mp->time00_toggle) {
			mp->time00_ref = p2p.wtime();
		} // if

		return(p2p.wtime() - mp->time00_ref);
	} // mp_time00

	inline void mp_abort( int reason, mp_handle h ) {
		P2PConnection & p2p = P2PConnection::instance();
		p2p.request(P2PTag::abort);
		p2p.send(&reason, 1);
	} // mp_abort

	inline void mp_barrier( mp_handle h ) {
		P2PConnection & p2p = P2PConnection::instance();
		p2p.request(P2PTag::barrier);
		p2p.sync();
	} // mp_barrier

	inline error_code mp_size_recv_buffer( int tag, int size, mp_handle h ) {
		return mp_size_recv_buffer_dmp(tag, size, h);
	} // mp_size_recv_buffer

	inline error_code mp_size_send_buffer( int tag, int size, mp_handle h ) {
		return mp_size_send_buffer_dmp(tag, size, h);
	} // mp_size_send_buffer

	inline error_code mp_begin_recv( int rbuf, int size, int sender,
		int tag, mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(rbuf<0 || rbuf>=NUM_BUF) { return ERROR_CODE("Bad recv_buf"); }
		if(size<=0) { return ERROR_CODE("Bad msg_size"); }
		if(sender<0 || sender>=mp->nproc) { return ERROR_CODE("Bad sender"); }
		if(mp->rbuf[rbuf]==NULL) { return ERROR_CODE("NULL recv_buf"); }
		if(mp->rbuf_size[rbuf]<size ) {
			return ERROR_CODE("recv_buf too small");
		} // if

		mp->rreq_size[rbuf] = size;

		// request an isend
		p2p.request(P2PTag::irecv);

		// send data non-blocking
		switch(p2p.irecv(static_cast<char *>(mp->rbuf[rbuf]), size,
			mp->rreq[rbuf])) {
			case MPI_SUCCESS:
				return NO_ERROR;
			case MPI_ERR_COMM:
				return ERROR_CODE("MPI_ERR_COMM");
			case MPI_ERR_COUNT:
				return ERROR_CODE("MPI_ERR_COUNT");
			case MPI_ERR_TYPE:
				return ERROR_CODE("MPI_ERR_TYPE");
			case MPI_ERR_TAG:
				return ERROR_CODE("MPI_ERR_TAG");
			case MPI_ERR_RANK:
				return ERROR_CODE("MPI_ERR_RANK");
			case MPI_ERR_OTHER:
				return ERROR_CODE("MPI_ERR_OTHER");
			default: break;
		} // switch

		return ERROR_CODE("Unknown MPI error");
	} // mp_begin_recv

	inline error_code mp_begin_send( int sbuf, int size, int receiver,
		int tag, mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(sbuf<0 || sbuf>=NUM_BUF) { return ERROR_CODE("Bad send_buf"); }
		if(size<=0) { return ERROR_CODE("Bad msg_size"); }
		if(receiver<0 || receiver>=mp->nproc) {
			return ERROR_CODE("Bad receiver");
		} // if
		if(mp->sbuf[sbuf]==NULL) { return ERROR_CODE("NULL send_buf"); }
		if(mp->sbuf_size[sbuf]<size ) {
			return ERROR_CODE("send_buf too small");
		} // if

		mp->sreq_size[sbuf] = size;

		// request an isend
		p2p.request(P2PTag::isend);

		// send data non-blocking
		switch(p2p.isend(static_cast<char *>(mp->sbuf[sbuf]), size,
			mp->sreq[sbuf])) {
			case MPI_SUCCESS:
				return NO_ERROR;
			case MPI_ERR_COMM:
				return ERROR_CODE("MPI_ERR_COMM");
			case MPI_ERR_COUNT:
				return ERROR_CODE("MPI_ERR_COUNT");
			case MPI_ERR_TYPE:
				return ERROR_CODE("MPI_ERR_TYPE");
			case MPI_ERR_TAG:
				return ERROR_CODE("MPI_ERR_TAG");
			case MPI_ERR_RANK:
				return ERROR_CODE("MPI_ERR_RANK");
			case MPI_ERR_OTHER:
				return ERROR_CODE("MPI_ERR_OTHER");
			default: break;
		} // switch

		return ERROR_CODE("Unknown MPI error");
	} // mp_begin_send

	inline error_code mp_end_recv( int rbuf, mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
		MPI_Status status;
		int size;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(rbuf<0 || rbuf>=NUM_BUF) { return ERROR_CODE("Bad recv_buf"); }

		switch(p2p.wait(mp->rreq[rbuf], status)) {
			case MPI_SUCCESS:
				break;
			case MPI_ERR_REQUEST:
				return ERROR_CODE("MPI_Wait - MPI_ERR_REQUEST");
			case MPI_ERR_ARG:
				return ERROR_CODE("MPI_Wait - MPI_ERR_ARG");
			default:
				return ERROR_CODE("MPI_Wait - Unknown MPI error");
		} // switch

		switch(p2p.get_count<char>(status, size)) {
			case MPI_SUCCESS:
				break;
			case MPI_ERR_ARG:
				return ERROR_CODE("MPI_Get_count - MPI_ERR_ARG");
			case MPI_ERR_TYPE:
				return ERROR_CODE("MPI_Get_count - MPI_ERR_TYPE");
			default:
				return ERROR_CODE("MPI_Wait - Unknown MPI error");
		} // switch

		if(mp->rreq_size[rbuf] != size) {
			return ERROR_CODE("Sizes do not match");
		} // if

		return NO_ERROR;
	} // mp_end_recv

	inline error_code mp_end_send( int sbuf, mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(sbuf<0 || sbuf>=NUM_BUF) { return ERROR_CODE("Bad send_buf"); }

		switch(p2p.wait(mp->sreq[sbuf], *MPI_STATUS_IGNORE)) {
			case MPI_SUCCESS:
				return NO_ERROR;
			case MPI_ERR_REQUEST:
				return ERROR_CODE("MPI_ERR_REQUEST");
			case MPI_ERR_ARG:
				return ERROR_CODE("MPI_ERR_ARG");
			default:
				break;
		} // switch

		return ERROR_CODE("Unknown MPI error");
	} // mp_end_send

	inline error_code mp_allsum_d( double *local, double *global,
		int n, mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(local==NULL) { return ERROR_CODE("Bad local"); }
		if(global==NULL) { return ERROR_CODE("Bad global"); }
		if(abs(local-global)<n) {
			return ERROR_CODE("Overlapping local and global");
		} // if

		p2p.request(P2PTag::allreduce_sum_double, n);
		p2p.send(local, n);
		p2p.recv(global, n);

		// need to fix error code propagation
		return NO_ERROR;
	} // mp_allsum_d

	inline error_code mp_allsum_i( int *local, int *global, int n,
		mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(local==NULL) { return ERROR_CODE("Bad local"); }
		if(global==NULL) { return ERROR_CODE("Bad global"); }
		if(abs(local-global)<n) {
			return ERROR_CODE("Overlapping local and global");
		} // if

		p2p.request(P2PTag::allreduce_sum_int, n);
		p2p.send(local, n);
		p2p.recv(global, n);

		// need to fix error code propagation
		return NO_ERROR;
	} // mp_allsum_i

	inline error_code mp_allgather_i( int *sbuf, int *rbuf, int n,
		mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(sbuf==NULL) { return ERROR_CODE("Bad send"); }
		if(rbuf==NULL) { return ERROR_CODE("Bad recv"); }
		if(n<1) { return ERROR_CODE("Bad n"); }

		p2p.request(P2PTag::allgather_int, n);
		p2p.send(sbuf, n);
		p2p.send(rbuf, n*p2p.size());

		// need to fix error code propagation
		return NO_ERROR;
	} // mp_allgather_i

	inline error_code mp_allgather_i64( int64_t *sbuf, int64_t *rbuf, int n,
		mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(sbuf==NULL) { return ERROR_CODE("Bad send"); }
		if(rbuf==NULL) { return ERROR_CODE("Bad recv"); }
		if(n<1) { return ERROR_CODE("Bad n"); }

		p2p.request(P2PTag::allgather_long_long, n);
		p2p.send(sbuf, n);
		p2p.send(rbuf, n*p2p.size());

		// need to fix error code propagation
		return NO_ERROR;
	} // mp_allgather_i64

	error_code mp_send_i( int *buf, int n, int dst, mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }

		p2p.request(P2PTag::send, n, dst);
		switch(p2p.send(buf, n)) {
			case MPI_SUCCESS:
				return NO_ERROR;
			case MPI_ERR_COMM:
				return ERROR_CODE("MPI_ERR_COMM");
			case MPI_ERR_COUNT:
				return ERROR_CODE("MPI_ERR_COUNT");
			case MPI_ERR_TYPE:
				return ERROR_CODE("MPI_ERR_TYPE");
			case MPI_ERR_OTHER:
				return ERROR_CODE("MPI_ERR_OTHER");
			default: break;
		} // switch

		return ERROR_CODE("Unknown MPI error");
	} // mp_send_i

	error_code mp_recv_i( int *buf, int n, int src, mp_handle h ) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }

		p2p.request(P2PTag::recv, n, src);
		switch(p2p.recv(buf, n)) {
			case MPI_SUCCESS:
				return NO_ERROR;
			case MPI_ERR_COMM:
				return ERROR_CODE("MPI_ERR_COMM");
			case MPI_ERR_COUNT:
				return ERROR_CODE("MPI_ERR_COUNT");
			case MPI_ERR_TYPE:
				return ERROR_CODE("MPI_ERR_TYPE");
			case MPI_ERR_OTHER:
				return ERROR_CODE("MPI_ERR_OTHER");
			default: break;
		} // switch

		return ERROR_CODE("Unknown MPI error");
	} // mp_recv_i

}; // struct AAISPolicy

#endif // AAISPolicy_hxx
