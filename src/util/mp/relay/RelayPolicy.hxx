#ifndef RelayPolicy_hxx
#define RelayPolicy_hxx

#include "../../util_base.h"
#include "../mp_handle.h"
#include "../../relay/ConnectionManager.hxx"
#include "../../relay/p2p/P2PConnection.hxx"

// FIXME TEMPORARY HACK
// FIXME: This is a hack to make the current processor rank available to the
// Error output messages during memory allocation.
extern int err_rank;

static const double RESIZE_FACTOR = 1.3125; // "silver ratio"

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

struct RelayPolicy {

    inline void mp_init(int argc, char ** argv) {
        ConnectionManager::instance().init(argc, argv);
    } // mp_init

    inline void mp_finalize(mp_handle h) {
        P2PConnection::instance().post(P2PTag::end);
        ConnectionManager::instance().finalize();
    } // mp_finalize

    inline mp_handle new_mp(void) {
        // get p2p connection
        P2PConnection & p2p = P2PConnection::instance();

        mp_t * mp;
        MALLOC( mp, 1 );
        CLEAR( mp, 1 );

        // get rank and size from point-to-point connection
        mp->rank = p2p.global_id();
        mp->nproc = p2p.global_size();

		// FIXME TEMPORARY HACK
		// FIXME: This is a hack to make the current processor rank
		// available to the Error output messages during memory allocation.
		err_rank = p2p.global_id();

        mp->elapsed_ref = mp_wtime();

        return static_cast<mp_handle>(mp);
    } // new_mp

    inline void delete_mp(mp_handle * h) {
        mp_t * mp;

        if(!h) return;

        mp = static_cast<mp_t *>(*h);
        if(!mp) return;

        mp->rank = -1;
        mp->nproc = 0;
        mp->elapsed_ref = 0;
        mp->time00_ref = 0;
        mp->time00_toggle = 0;

        // max_buffers is a global defined in <mp_t.h>
        for(size_t i(0); i<max_buffers; i++) {
            FREE_ALIGNED( mp->rbuf[i] );
            FREE_ALIGNED( mp->sbuf[i] );
            mp->rbuf_size[i] = 0;
            mp->sbuf_size[i] = 0;
            mp->rreq[i] = 0;
            mp->sreq[i] = 0;
            mp->rreq_size[i] = 0;
            mp->sreq_size[i] = 0;
        } // for

        FREE(mp);
        *h = NULL;
    } // delete_mp

    inline int mp_rank(mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        return mp ? mp->rank : -1;  
    } // mp_rank

    inline int mp_nproc(mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        return mp ? mp->nproc : -1; 
    } // mp_nproc

    inline void * ALIGNED(16) mp_recv_buffer(int tag, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        if(!mp || tag<0 || uint64_t(tag)>=max_buffers) { return NULL; }
        return mp->rbuf[tag];
    } // mp_recv_buffer

    inline void * ALIGNED(16) mp_send_buffer(int tag, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        if(!mp || tag<0 || uint64_t(tag)>=max_buffers) { return NULL; }
        return mp->sbuf[tag];
    } // mp_send_buffer

    inline double mp_elapsed(mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();
        double time = mp_wtime() - mp->elapsed_ref;

        MPRequest request(P2PTag::allreduce_max_double, P2PTag::data, 1, 0);
        p2p.post(request);
        p2p.send(&time, request.count, request.tag);
		double recv_time;
        p2p.recv(&recv_time, request.count, request.tag, request.id);

        return time;
    } // mp_elapsed

    inline double mp_time00(mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        
        if(!mp) { return -1; }

        mp->time00_toggle ^= 1;

        if(mp->time00_toggle) mp->time00_ref = mp_wtime(); // mp->time00_ref = p2p.wtime();

        return mp_wtime() - mp->time00_ref; // return p2p.wtime() - mp->time00_ref;
    } // mp_time00

    inline double mp_wtime() {
        P2PConnection & p2p = P2PConnection::instance();
        double wtime;
        MPRequest request(P2PTag::wtime, P2PTag::data, 1, 0);
        p2p.post(request);
        p2p.recv(&wtime, request.count, request.tag, request.id);
        return wtime;
    } // mp_wtime

    inline void mp_abort(int reason, mp_handle h) {
        P2PConnection & p2p = P2PConnection::instance();
        MPRequest request(P2PTag::abort, P2PTag::data, 1, 0);
        p2p.post(request);
        p2p.send(&reason, 1, P2PTag::data);
        p2p.abort(reason);
    } // mp_abort

    inline void mp_barrier(mp_handle h) {
        P2PConnection & p2p = P2PConnection::instance();
        p2p.post(P2PTag::barrier);
        p2p.barrier();
    } // mp_barrier

    inline void mp_size_recv_buffer(int tag, int size, mp_handle h) {
        mp_t *mp = (mp_t *)h;
        char * ALIGNED(16) rbuf;

        // Check input arguments
        if( mp==NULL            ) ERROR(( "Bad handle" ));
        if( tag<0 || uint64_t(tag)>=max_buffers ) ERROR(( "Bad tag" ));
        if( size<=0               ) ERROR(( "Bad size" ));

        // If no buffer allocated for this tag
        if( mp->rbuf[tag]==NULL ) {
          MALLOC_ALIGNED( mp->rbuf[tag], size, 16 );
          mp->rbuf_size[tag] = size;
          return;
        }

        // Is there already a large enough buffer
        if( mp->rbuf_size[tag]>=size ) return;

        // Try to reduce the number of realloc calls

        size = static_cast<int>(size*RESIZE_FACTOR);

        // Create the new recv buffer

        MALLOC_ALIGNED( rbuf, size, 16 );

        // Preserve the old recv buffer data
        // FIMXE: SILLY
        COPY( rbuf, mp->rbuf[tag], mp->rbuf_size[tag] );

        // Free the old recv buffer

        FREE_ALIGNED( mp->rbuf[tag] );
        mp->rbuf[tag]      = rbuf;
        mp->rbuf_size[tag] = size;
    } // mp_size_recv_buffer

    inline void mp_size_send_buffer(int tag, int size, mp_handle h) {
        mp_t *mp = (mp_t *)h;
        char * ALIGNED(16) sbuf;

        // Check input arguments
        if( mp==NULL              ) ERROR(( "Bad handle" ));
        if( tag<0 || uint64_t(tag)>=max_buffers ) ERROR(( "Bad tag" ));
        if( size<=0               ) ERROR(( "Bad size" ));

        // If no buffer allocated for this tag
        if( mp->sbuf[tag]==NULL ) {
          MALLOC_ALIGNED( mp->sbuf[tag], size, 16 );
          mp->sbuf_size[tag] = size;
          return;
        }

        // Is there already a large enough buffer
        if( mp->sbuf_size[tag]>=size ) return;

        // Try to reduce the number of realloc calls

        size = static_cast<int>(size*RESIZE_FACTOR);

        // Create the new send buffer

        MALLOC_ALIGNED( sbuf, size, 16 );

        // Preserve the old send buffer data
        // FIXME: SILLY
        COPY( sbuf, mp->sbuf[tag], mp->sbuf_size[tag] );

        // Free the old recv buffer

        FREE_ALIGNED( mp->sbuf[tag] );
        mp->sbuf[tag]      = sbuf;
        mp->sbuf_size[tag] = size;
    } // mp_size_send_buffer

    inline void mp_begin_recv(int rbuf, int size, int sender,
        int tag, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();

        if(mp==NULL) ERROR(( "Bad handle" ));
        if(rbuf<0 || uint64_t(rbuf)>=max_buffers) ERROR(( "Bad recv_buf" ));
        if(size<=0) ERROR(( "Bad msg_size" ));
        if(sender<0 || sender>=mp->nproc) ERROR(( "Bad sender" ));
        if(mp->rbuf[rbuf]==NULL) ERROR(( "NULL recv_buf" ));
        if(mp->rbuf_size[rbuf]<size) ERROR(( "recv_buf too small" ));

        mp->rreq_size[rbuf] = size;

        MPRequest request(P2PTag::irecv, tag, size, rbuf, sender);
        p2p.post(request);

        /*
        !!!FIXME!!!
        int errcode =
            p2p.irecv(static_cast<char *>(mp->rbuf[rbuf]), size, tag, rbuf);
        */
        p2p.irecv(static_cast<char *>(mp->rbuf[rbuf]), size, tag, rbuf);

        // NEED TO DEAL WITH ERRORS!!!

        /*
        switch(p2p.irecv(static_cast<char *>(mp->rbuf[rbuf]), size,
            tag, rbuf)) {
            case MPI_SUCCESS:
//      std::cerr << "WRAPPER: p2p.irecv succeeded" << std::endl;
                return;
            case MPI_ERR_COMM:
                ERROR(( "MPI_ERR_COMM" ));
            case MPI_ERR_COUNT:
                ERROR(( "MPI_ERR_COUNT" ));
            case MPI_ERR_TYPE:
                ERROR(( "MPI_ERR_TYPE" ));
            case MPI_ERR_TAG:
                ERROR(( "MPI_ERR_TAG" ));
            case MPI_ERR_RANK:
                ERROR(( "MPI_ERR_RANK" ));
            case MPI_ERR_OTHER:
                ERROR(( "MPI_ERR_OTHER" ));
            default: break;
        } // switch

        ERROR(( "Unknown MPI error" ));
        */
    } // mp_begin_recv

    inline void mp_begin_send(int sbuf, int size, int receiver,
        int tag, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();

        if(mp==NULL) ERROR(( "Bad handle" ));
        if(sbuf<0 || uint64_t(sbuf)>=max_buffers) ERROR(( "Bad send_buf" ));
        if(size<=0) ERROR(( "Bad msg_size" ));
        if(receiver<0 || receiver>=mp->nproc) ERROR(( "Bad receiver" ));
        if(mp->sbuf[sbuf]==NULL) ERROR(( "NULL send_buf" ));
        if(mp->sbuf_size[sbuf]<size) ERROR(( "send_buf too small" ));

        mp->sreq_size[sbuf] = size;

        MPRequest request(P2PTag::isend, tag, size, sbuf, receiver);
        p2p.post(request);

        /*
        !!!FIXME!!!
        int errcode = p2p.isend(static_cast<char *>(mp->sbuf[sbuf]),
            size, tag, sbuf);
        */
        p2p.isend(static_cast<char *>(mp->sbuf[sbuf]), size, tag, sbuf);

        /*
        switch(p2p.isend(static_cast<char *>(mp->sbuf[sbuf]), size,
            tag, sbuf)) {
            case MPI_SUCCESS:
                return;
            case MPI_ERR_COMM:
                ERROR(( "MPI_ERR_COMM" ));
            case MPI_ERR_COUNT:
                ERROR(( "MPI_ERR_COUNT" ));
            case MPI_ERR_TYPE:
                ERROR(( "MPI_ERR_TYPE" ));
            case MPI_ERR_TAG:
                ERROR(( "MPI_ERR_TAG" ));
            case MPI_ERR_RANK:
                ERROR(( "MPI_ERR_RANK" ));
            case MPI_ERR_OTHER:
                ERROR(( "MPI_ERR_OTHER" ));
            default: break;
        } // switch

        ERROR(( "Unknown MPI error" ));
        */
    } // mp_begin_send

    inline void mp_end_recv(int rbuf, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();
        int size;
        //std::cerr << "WRAPPER: mp_end_recv called" << std::endl;

        if(mp==NULL) ERROR(( "Bad handle" ));
        if(rbuf<0 || uint64_t(rbuf)>=max_buffers) ERROR(( "Bad recv_buf" ));

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
                ERROR(( "MPI_Wait - MPI_ERR_REQUEST" ));
            case MPI_ERR_ARG:
                ERROR(( "MPI_Wait - MPI_ERR_ARG" ));
            default:
                ERROR(( "MPI_Wait - Unknown MPI error" ));
        } // switch
*/

        errcode = p2p.get_count<char>(rbuf, size);

/*
        switch(p2p.get_count<char>(rbuf, size)) {
            case MPI_SUCCESS:
        //std::cout << "WRAPPER: end get count" << std::endl;
                break;
            case MPI_ERR_ARG:
                ERROR(( "MPI_Get_count - MPI_ERR_ARG" ));
            case MPI_ERR_TYPE:
                ERROR(( "MPI_Get_count - MPI_ERR_TYPE" ));
            default:
                ERROR(( "MPI_Wait - Unknown MPI error" ));
        } // switch

        if(mp->rreq_size[rbuf] != size) ERROR(( "Sizes do not match" ));
*/

    } // mp_end_recv

    inline void mp_end_send(int sbuf, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();
        //std::cerr << "WRAPPER: mp_end_send called" << std::endl;

        if(mp==NULL) ERROR(( "Bad handle" ));
        if(sbuf<0 || uint64_t(sbuf)>=max_buffers)
            ERROR(( "Bad send_buf" ));

        MPRequest request(P2PTag::wait_send, 0, 0, sbuf);
        p2p.post(request);

/*
    FIXME
        int errcode = p2p.wait_send(sbuf);
*/
        p2p.wait_send(sbuf);

/*
        switch(p2p.wait_send(sbuf)) {
            case MPI_SUCCESS:
                return;
            case MPI_ERR_REQUEST:
                ERROR(( "MPI_ERR_REQUEST" ));
            case MPI_ERR_ARG:
                ERROR(( "MPI_ERR_ARG" ));
            default:
                break;
        } // switch

        ERROR(( "Unknown MPI error" ));
*/
    } // mp_end_send

    inline void mp_allsum_d(double *local, double *global,
        int n, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();

        if(mp==NULL) ERROR(( "Bad handle" ));
        if(local==NULL) ERROR(( "Bad local" ));
        if(global==NULL) ERROR(( "Bad global" ));
        if(abs(local-global)<n) ERROR(( "Overlapping local and global" ));

        MPRequest request(P2PTag::allreduce_sum_double, P2PTag::data, n, 0);
        p2p.post(request);
        p2p.send(local, request.count, request.tag);
        p2p.recv(global, request.count, request.tag, request.id);

        // need to fix error code propagation
    } // mp_allsum_d

    inline void mp_allsum_i(int *local, int *global, int n,
        mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();

        if(mp==NULL) ERROR(( "Bad handle" ));
        if(local==NULL) ERROR(( "Bad local" ));
        if(global==NULL) ERROR(( "Bad global" ));
        if(abs(local-global)<n) ERROR(( "Overlapping local and global" ));

        MPRequest request(P2PTag::allreduce_sum_int, P2PTag::data, n, 0);
        p2p.post(request);
        p2p.send(local, request.count, request.tag);
        p2p.recv(global, request.count, request.tag, request.id);

        // need to fix error code propagation
    } // mp_allsum_i

	inline void mp_gather_uc(unsigned char * sbuf, unsigned char * rbuf,
		int n, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();

        if(mp==NULL) ERROR(( "Bad handle" ));
        if(sbuf==NULL) ERROR(( "Bad send" ));
        if(mp->rank == 0 && rbuf==NULL) ERROR(( "Bad recv" ));
        if(n<1) ERROR(( "Bad n" ));

        MPRequest request(P2PTag::gather_uc, P2PTag::data, n, 0);
        p2p.post(request);
        p2p.send(sbuf, request.count, request.tag);

		if(mp->rank == 0) {
        	p2p.recv(rbuf, request.count*p2p.global_size(),
				request.tag, request.id);
		} // if
	} // mp_gather_uc

    inline void mp_allgather_i(int *sbuf, int *rbuf, int n,
        mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();

        if(mp==NULL) ERROR(( "Bad handle" ));
        if(sbuf==NULL) ERROR(( "Bad send" ));
        if(rbuf==NULL) ERROR(( "Bad recv" ));
        if(n<1) ERROR(( "Bad n" ));

        MPRequest request(P2PTag::allgather_int, P2PTag::data, n, 0);
        p2p.post(request);
        p2p.send(sbuf, request.count, request.tag);
        p2p.recv(rbuf, request.count*p2p.global_size(),
            request.tag, request.id);

        // need to fix error code propagation
    } // mp_allgather_i

    inline void mp_allgather_i64(int64_t *sbuf, int64_t *rbuf, int n,
        mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();

        if(mp==NULL) ERROR(( "Bad handle" ));
        if(sbuf==NULL) ERROR(( "Bad send" ));
        if(rbuf==NULL) ERROR(( "Bad recv" ));
        if(n<1) ERROR(( "Bad n" ));

        MPRequest request(P2PTag::allgather_int64, P2PTag::data, n, 0);
        p2p.post(request);
        p2p.send(sbuf, request.count, request.tag);
        p2p.recv(rbuf, request.count*p2p.global_size(),
            request.tag, request.id);

        // need to fix error code propagation
    } // mp_allgather_i64

    void mp_send_i(int *buf, int n, int dst, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();

        if(mp==NULL) ERROR(( "Bad handle" ));

        MPRequest request(P2PTag::send, P2PTag::data, n, 0, dst);
        p2p.post(request);

        /*
        FIXME
        int errcode = p2p.send(buf, request.count, request.tag);
        */
        p2p.send(buf, request.count, request.tag);

/*
        switch(p2p.send(buf, request.count, request.tag)) {
            case MPI_SUCCESS:
                return;
            case MPI_ERR_COMM:
                ERROR(( "MPI_ERR_COMM" ));
            case MPI_ERR_COUNT:
                ERROR(( "MPI_ERR_COUNT" ));
            case MPI_ERR_TYPE:
                ERROR(( "MPI_ERR_TYPE" ));
            case MPI_ERR_OTHER:
                ERROR(( "MPI_ERR_OTHER" ));
            default: break;
        } // switch

        ERROR(( "Unknown MPI error" ));
*/
    } // mp_send_i

    void mp_recv_i(int *buf, int n, int src, mp_handle h) {
        mp_t * mp = static_cast<mp_t *>(h);
        P2PConnection & p2p = P2PConnection::instance();

        if(mp==NULL) ERROR(( "Bad handle" ));

        MPRequest request(P2PTag::recv, P2PTag::data, n, 0, src);
        p2p.post(request);

        /*
        FIXME
        int errcode = p2p.recv(buf, request.count, request.tag, request.id);
        */
        p2p.recv(buf, request.count, request.tag, request.id);

/*
        switch(p2p.recv(buf, request.count, request.tag, request.id)) {
            case MPI_SUCCESS:
                return;
            case MPI_ERR_COMM:
                ERROR(( "MPI_ERR_COMM" ));
            case MPI_ERR_COUNT:
                ERROR(( "MPI_ERR_COUNT" ));
            case MPI_ERR_TYPE:
                ERROR(( "MPI_ERR_TYPE" ));
            case MPI_ERR_OTHER:
                ERROR(( "MPI_ERR_OTHER" ));
            default: break;
        } // switch

        ERROR(( "Unknown MPI error" ));
*/
    } // mp_recv_i

}; // struct RelayPolicy

#endif // RelayPolicy_hxx
