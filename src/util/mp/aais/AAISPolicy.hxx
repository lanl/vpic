#ifndef AAISPolicy_hxx
#define AAISPolicy_hxx

#include <util_base.h>
#include <mp_handle.h>
#include <ConnectionManager.hxx>
#include <P2PConnection.hxx>

static const double RESIZE_FACTOR = 1.1;

struct mp_t {
	int rank;
	int nproc;
	double elapsed_ref;
	double time00_ref;
	int time00_toggle;
	char * ALIGNED(16) rbuf[max_buffers];
	char * ALIGNED(16) sbuf[max_buffers];
	int rbuf_size[max_buffers];
	int sbuf_size[max_buffers];
	int rreq[max_buffers];
	int sreq[max_buffers];
	int rreq_size[max_buffers];
	int sreq_size[max_buffers];
}; // struct mp

struct AAISPolicy {

	inline void mp_init(int argc, char ** argv) {
		ConnectionManager::instance().init(argc, argv);
//		std::cerr << "WRAPPER: mp_init called" << std::endl;
	} // mp_init

	inline void mp_finalize(mp_handle h) {
		P2PConnection::instance().post(P2PTag::end);
		ConnectionManager::instance().finalize();
//		std::cerr << "WRAPPER: mp_finalize called" << std::endl;
	} // mp_finalize

	inline mp_handle new_mp(void) {
		// get p2p connection
		P2PConnection & p2p = P2PConnection::instance();

		mp_t * mp = static_cast<mp_t *>(malloc(sizeof(mp_t)));
		if(!mp) { return static_cast<mp_handle>(NULL); }

		// get rank and size from point-to-point connection
		mp->rank = p2p.global_id();
		mp->nproc = p2p.global_size();

//		std::cerr << "WRAPPER: new_mp called for rank " <<
//			mp->rank << std::endl;

		// get wtime from point-to-point connection
		//mp->elapsed_ref = p2p.wtime();
		mp->elapsed_ref = mp_wtime();
		mp->time00_ref = 0;
		mp->time00_toggle = 0;

		// max_buffers is a global defined in <mp_t.h>
		for(size_t i(0); i<max_buffers; i++) {
			mp->rbuf[i] = NULL;
			mp->sbuf[i] = NULL;
			mp->rbuf_size[i] = 0;
			mp->sbuf_size[i] = 0;
			mp->rreq[i] = 0;
			mp->sreq[i] = 0;
			mp->rreq_size[i] = 0;
			mp->sreq_size[i] = 0;
		} // for

		return static_cast<mp_handle>(mp);
	} // new_mp

	inline void delete_mp(mp_handle * h) {
//		std::cerr << "WRAPPER: delete_mp called" << std::endl;
		mp_t * mp;

		if(!h) { return; }

		mp = static_cast<mp_t *>(*h);

		if(!mp) { return; }

		mp->rank = -1;
		mp->nproc = 0;
		mp->elapsed_ref = 0;
		mp->time00_ref = 0;
		mp->time00_toggle = 0;

		// max_buffers is a global defined in <mp_t.h>
		for(size_t i(0); i<max_buffers; i++) {
			// allocated with malloc
			if(mp->rbuf[i]) { free_aligned(mp->rbuf[i]); }
			if(mp->sbuf[i]) { free_aligned(mp->sbuf[i]); }

			mp->rbuf[i] = NULL;
			mp->sbuf[i] = NULL;
			mp->rbuf_size[i] = 0;
			mp->sbuf_size[i] = 0;
			mp->rreq[i] = 0;
			mp->sreq[i] = 0;
			mp->rreq_size[i] = 0;
			mp->sreq_size[i] = 0;
		} // for

		free(mp);
		*h = NULL;
	} // delete_mp

	inline int mp_rank(mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
//		std::cerr << "WRAPPER: mp_rank called" << std::endl;
		return mp ? mp->rank : -1;	
	} // mp_rank

	inline int mp_nproc(mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
//		std::cerr << "WRAPPER: mp_nproc called" << std::endl;
		return mp ? mp->nproc : -1;	
	} // mp_nproc

	inline void * ALIGNED(16) mp_recv_buffer(int tag, mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
//		std::cerr << "WRAPPER: mp_recv_buffer called" << std::endl;
		if(!mp || tag<0 || uint64_t(tag)>=max_buffers) { return NULL; }
		return mp->rbuf[tag];
	} // mp_recv_buffer

	inline void * ALIGNED(16) mp_send_buffer(int tag, mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
//		std::cerr << "WRAPPER: mp_send_buffer called" << std::endl;
		if(!mp || tag<0 || uint64_t(tag)>=max_buffers) { return NULL; }
		return mp->sbuf[tag];
	} // mp_send_buffer

	inline double mp_elapsed(mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_elapsed called" << std::endl;
		//double time = p2p.wtime() - mp->elapsed_ref;
		double time = mp_wtime() - mp->elapsed_ref;

		MPRequest request(P2PTag::allreduce_max_double, P2PTag::data, 1, 0);
		p2p.post(request);
		p2p.send(&time, request.count, request.tag);
		p2p.recv(&time, request.count, request.tag, request.id);

		return time;
	} // mp_elapsed

	inline double mp_time00(mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		//P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_time00 called" << std::endl;
		
		if(!mp) { return -1; }

		mp->time00_toggle ^= 1;

		if(mp->time00_toggle) {
			//mp->time00_ref = p2p.wtime();
			mp->time00_ref = mp_wtime();
		} // if

		//return(p2p.wtime() - mp->time00_ref);
		return(mp_wtime() - mp->time00_ref);
	} // mp_time00

	inline double mp_wtime() {
		P2PConnection & p2p = P2PConnection::instance();
		double wtime;
		//std::cout << "requesting wtime " << std::endl;
		MPRequest request(P2PTag::wtime, P2PTag::data, 1, 0);
		p2p.post(request);
		p2p.recv(&wtime, request.count, request.tag, request.id);
		//std::cout << "received wtime " << std::endl;
		return wtime;
	} // mp_wtime

	inline void mp_abort(int reason, mp_handle h) {
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_abort called" << std::endl;
		MPRequest request(P2PTag::abort, P2PTag::data, 1, 0);
		p2p.post(request);
		p2p.send(&reason, 1, P2PTag::data);
		p2p.abort(reason);
	} // mp_abort

	inline void mp_barrier(mp_handle h) {
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_barrier called" << std::endl;
		p2p.post(P2PTag::barrier);
		p2p.barrier();
	} // mp_barrier

	inline error_code mp_size_recv_buffer(int tag, int size, mp_handle h) {
		mp_t *mp = (mp_t *)h;
		char * ALIGNED(16) rbuf;

		// Check input arguments
		if( mp==NULL            ) return ERROR_CODE("Bad handle");
		if( tag<0 || uint64_t(tag)>=max_buffers ) return ERROR_CODE("Bad tag");
		if( size<=0               ) return ERROR_CODE("Bad size");

		// If no buffer allocated for this tag
		if( mp->rbuf[tag]==NULL ) {
			mp->rbuf[tag] = (char * ALIGNED(16))
			malloc_aligned( size, 16 );
			if( mp->rbuf[tag]==NULL ) {
				return ERROR_CODE("malloc_aligned failed");
			}
			mp->rbuf_size[tag] = size;
			return NO_ERROR;
		}

		// Is there already a large enough buffer
		if( mp->rbuf_size[tag]>=size ) return NO_ERROR;

		// Try to reduce the number of realloc calls

		size = static_cast<int>(size*RESIZE_FACTOR);

		// Create the new recv buffer

		rbuf = (char * ALIGNED(16))malloc_aligned( size, 16 );
		if( rbuf==NULL ) return ERROR_CODE("malloc_aligned failed");

		// Preserve the old recv buffer data

		memcpy( rbuf, mp->rbuf[tag], mp->rbuf_size[tag] );

		// Free the old recv buffer

		free_aligned( mp->rbuf[tag] );
		mp->rbuf[tag]      = rbuf;
		mp->rbuf_size[tag] = size;

		return NO_ERROR;
	} // mp_size_recv_buffer

	inline error_code mp_size_send_buffer(int tag, int size, mp_handle h) {
		mp_t *mp = (mp_t *)h;
		char * ALIGNED(16) sbuf;

		// Check input arguments
		if( mp==NULL              ) return ERROR_CODE("Bad handle");
		if( tag<0 || uint64_t(tag)>=max_buffers ) return ERROR_CODE("Bad tag");
		if( size<=0               ) return ERROR_CODE("Bad size");

		// If no buffer allocated for this tag
		if( mp->sbuf[tag]==NULL ) {
			mp->sbuf[tag] = (char * ALIGNED(16))
			malloc_aligned( size, 16 );
			if( mp->sbuf[tag]==NULL ) {
				return ERROR_CODE("malloc_aligned failed");
			}
			mp->sbuf_size[tag] = size;
			return NO_ERROR;
		}

		// Is there already a large enough buffer
		if( mp->sbuf_size[tag]>=size ) return NO_ERROR;

		// Try to reduce the number of realloc calls

		size = static_cast<int>(size*RESIZE_FACTOR);

		// Create the new send buffer

		sbuf = (char * ALIGNED(16))malloc_aligned( size, 16 );
		if( sbuf==NULL ) return ERROR_CODE("malloc_aligned failed");

		// Preserve the old send buffer data

		memcpy( sbuf, mp->sbuf[tag], mp->sbuf_size[tag] );

		// Free the old recv buffer

		free_aligned( mp->sbuf[tag] );
		mp->sbuf[tag]      = sbuf;
		mp->sbuf_size[tag] = size;

		return NO_ERROR;
	} // mp_size_send_buffer

	inline error_code mp_begin_recv(int rbuf, int size, int sender,
		int tag, mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_begin_recv called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(rbuf<0 || uint64_t(rbuf)>=max_buffers)
			{ return ERROR_CODE("Bad recv_buf"); }
		if(size<=0) { return ERROR_CODE("Bad msg_size"); }
		if(sender<0 || sender>=mp->nproc) { return ERROR_CODE("Bad sender"); }
		if(mp->rbuf[rbuf]==NULL) { return ERROR_CODE("NULL recv_buf"); }
		if(mp->rbuf_size[rbuf]<size) {
			return ERROR_CODE("recv_buf too small");
		} // if

		mp->rreq_size[rbuf] = size;

		/*
		std::cerr << "WRAPPER: p2p.irecv rank " << p2p.global_id() <<
			" size " << size << " from " << sender <<
			" with tag " << tag << std::endl;
		*/

		MPRequest request(P2PTag::irecv, tag, size, rbuf, sender);
		p2p.post(request);

		
		/*
		!!!FIXME!!!
		int errcode =
			p2p.irecv(static_cast<char *>(mp->rbuf[rbuf]), size, tag, rbuf);
		*/
		p2p.irecv(static_cast<char *>(mp->rbuf[rbuf]), size, tag, rbuf);
		return NO_ERROR;

		// NEED TO DEAL WITH ERRORS!!!

		/*
		switch(p2p.irecv(static_cast<char *>(mp->rbuf[rbuf]), size,
			tag, rbuf)) {
			case MPI_SUCCESS:
//		std::cerr << "WRAPPER: p2p.irecv succeeded" << std::endl;
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
		*/
	} // mp_begin_recv

	inline error_code mp_begin_send(int sbuf, int size, int receiver,
		int tag, mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_begin_send called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(sbuf<0 || uint64_t(sbuf)>=max_buffers)
			{ return ERROR_CODE("Bad send_buf"); }
		if(size<=0) { return ERROR_CODE("Bad msg_size"); }
		if(receiver<0 || receiver>=mp->nproc) {
			return ERROR_CODE("Bad receiver");
		} // if
		if(mp->sbuf[sbuf]==NULL) { return ERROR_CODE("NULL send_buf"); }
		if(mp->sbuf_size[sbuf]<size) {
			return ERROR_CODE("send_buf too small");
		} // if

		mp->sreq_size[sbuf] = size;

		/*
		std::cerr << "WRAPPER: p2p.isend rank " << p2p.global_id() <<
			" size " << size << " to " << receiver <<
			" with tag " << tag << std::endl;
		*/

		MPRequest request(P2PTag::isend, tag, size, sbuf, receiver);
		p2p.post(request);

		/*
		!!!FIXME!!!
		int errcode = p2p.isend(static_cast<char *>(mp->sbuf[sbuf]),
			size, tag, sbuf);
		*/
		p2p.isend(static_cast<char *>(mp->sbuf[sbuf]), size, tag, sbuf);
		return NO_ERROR;

		/*
		switch(p2p.isend(static_cast<char *>(mp->sbuf[sbuf]), size,
			tag, sbuf)) {
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
		*/
	} // mp_begin_send

	inline error_code mp_end_recv(int rbuf, mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
		int size;
		//std::cerr << "WRAPPER: mp_end_recv called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(rbuf<0 || uint64_t(rbuf)>=max_buffers)
			{ return ERROR_CODE("Bad recv_buf"); }

		//std::cout << "WRAPPER: begin wait " << p2p.rank() << std::endl;

		MPRequest request(P2PTag::wait_recv, 0, 0, rbuf);
		p2p.post(request);

		int errcode = p2p.wait_recv(rbuf);

/*
		switch(p2p.wait_recv(rbuf)) {
			case MPI_SUCCESS:
		//std::cout << "WRAPPER: end wait" << std::endl;
				break;
			case MPI_ERR_REQUEST:
				return ERROR_CODE("MPI_Wait - MPI_ERR_REQUEST");
			case MPI_ERR_ARG:
				return ERROR_CODE("MPI_Wait - MPI_ERR_ARG");
			default:
				return ERROR_CODE("MPI_Wait - Unknown MPI error");
		} // switch
*/

		errcode = p2p.get_count<char>(rbuf, size);

/*
		switch(p2p.get_count<char>(rbuf, size)) {
			case MPI_SUCCESS:
		//std::cout << "WRAPPER: end get count" << std::endl;
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
*/

		return NO_ERROR;
	} // mp_end_recv

	inline error_code mp_end_send(int sbuf, mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
		//std::cerr << "WRAPPER: mp_end_send called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(sbuf<0 || uint64_t(sbuf)>=max_buffers)
			{ return ERROR_CODE("Bad send_buf"); }

		MPRequest request(P2PTag::wait_send, 0, 0, sbuf);
		p2p.post(request);

/*
	FIXME
		int errcode = p2p.wait_send(sbuf);
*/
		p2p.wait_send(sbuf);
		return NO_ERROR;

/*
		switch(p2p.wait_send(sbuf)) {
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
*/
	} // mp_end_send

	inline error_code mp_allsum_d(double *local, double *global,
		int n, mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_allsum_d called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(local==NULL) { return ERROR_CODE("Bad local"); }
		if(global==NULL) { return ERROR_CODE("Bad global"); }
		if(abs(local-global)<n) {
			return ERROR_CODE("Overlapping local and global");
		} // if

		MPRequest request(P2PTag::allreduce_sum_double, P2PTag::data, n, 0);
		p2p.post(request);
		p2p.send(local, request.count, request.tag);
		p2p.recv(global, request.count, request.tag, request.id);

		// need to fix error code propagation
		return NO_ERROR;
	} // mp_allsum_d

	inline error_code mp_allsum_i(int *local, int *global, int n,
		mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_allsum_i called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(local==NULL) { return ERROR_CODE("Bad local"); }
		if(global==NULL) { return ERROR_CODE("Bad global"); }
		if(abs(local-global)<n) {
			return ERROR_CODE("Overlapping local and global");
		} // if

		MPRequest request(P2PTag::allreduce_sum_int, P2PTag::data, n, 0);
		p2p.post(request);
		p2p.send(local, request.count, request.tag);
		p2p.recv(global, request.count, request.tag, request.id);

		// need to fix error code propagation
		return NO_ERROR;
	} // mp_allsum_i

	inline error_code mp_allgather_i(int *sbuf, int *rbuf, int n,
		mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_allgather_i called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(sbuf==NULL) { return ERROR_CODE("Bad send"); }
		if(rbuf==NULL) { return ERROR_CODE("Bad recv"); }
		if(n<1) { return ERROR_CODE("Bad n"); }

		MPRequest request(P2PTag::allgather_int, P2PTag::data, n, 0);
		p2p.post(request);
		p2p.send(sbuf, request.count, request.tag);
		p2p.recv(rbuf, request.count*p2p.global_size(),
			request.tag, request.id);

		// need to fix error code propagation
		return NO_ERROR;
	} // mp_allgather_i

	inline error_code mp_allgather_i64(int64_t *sbuf, int64_t *rbuf, int n,
		mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_allgather_i64 called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }
		if(sbuf==NULL) { return ERROR_CODE("Bad send"); }
		if(rbuf==NULL) { return ERROR_CODE("Bad recv"); }
		if(n<1) { return ERROR_CODE("Bad n"); }

		MPRequest request(P2PTag::allgather_int64, P2PTag::data, n, 0);
		p2p.post(request);
		p2p.send(sbuf, request.count, request.tag);
		p2p.recv(rbuf, request.count*p2p.global_size(),
			request.tag, request.id);

		// need to fix error code propagation
		return NO_ERROR;
	} // mp_allgather_i64

	error_code mp_send_i(int *buf, int n, int dst, mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_send_i called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }

		MPRequest request(P2PTag::send, P2PTag::data, n, 0, dst);
		p2p.post(request);

		/*
		FIXME
		int errcode = p2p.send(buf, request.count, request.tag);
		*/
		p2p.send(buf, request.count, request.tag);
		return NO_ERROR;

/*
		switch(p2p.send(buf, request.count, request.tag)) {
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
*/
	} // mp_send_i

	error_code mp_recv_i(int *buf, int n, int src, mp_handle h) {
		mp_t * mp = static_cast<mp_t *>(h);
		P2PConnection & p2p = P2PConnection::instance();
//		std::cerr << "WRAPPER: mp_recv_i called" << std::endl;

		if(mp==NULL) { return ERROR_CODE("Bad handle"); }

		MPRequest request(P2PTag::recv, P2PTag::data, n, 0, src);
		p2p.post(request);

		/*
		FIXME
		int errcode = p2p.recv(buf, request.count, request.tag, request.id);
		*/
		p2p.recv(buf, request.count, request.tag, request.id);
		return NO_ERROR;

/*
		switch(p2p.recv(buf, request.count, request.tag, request.id)) {
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
*/
	} // mp_recv_i

}; // struct AAISPolicy

#endif // AAISPolicy_hxx
