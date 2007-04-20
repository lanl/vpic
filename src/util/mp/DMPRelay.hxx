/*
	Definition of DMPRelay class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef DMPRelay_hxx
#define DMPRelay_hxx

#include <P2PData.hxx>

const size_t default_buffer_size(4096);

/*!
	\class DMPRelay DMPRelay.h
	\brief  provides...
*/
class DMPRelay
	{
	public:

		//! Constructor
		DMPRelay()
			: cbuf_size_(default_buffer_size), cbuf_(new char[cbuf_size_]),
			ibuf_size_(default_buffer_size), ibuf_(new int[ibuf_size_]),
			lbuf_size_(default_buffer_size), lbuf_(new long long[lbuf_size_]),
			dbuf_size_(default_buffer_size), dbuf_(new double[dbuf_size_])
			{}

		//! Destructor
		~DMPRelay()
			{
			delete cbuf_;
			delete ibuf_;
			delete lbuf_;
			delete dbuf_;
			} // ~DMPRelay

		void start();

	private:

		template<typename T> void resize_buf(int count,
			int & current_count, T *& buffer);

		int cbuf_size_;
		char * cbuf_;

		int ibuf_size_;
		int * ibuf_;

		int lbuf_size_;
		long long * lbuf_;

		int dbuf_size_;
		double * dbuf_;

	}; // class DMPRelay

void DMPRelay::start()
	{
		P2PConnection & p2p = P2PConnection::instance();

		bool relay(true);
		P2PHeader header;
		MPI_Status status;
		int * tmp(NULL);

		while(relay) {
			switch(p2p.poll(header)) {
				case P2PTag::send:
				case P2PTag::isend:
					resize_buf(header.count, cbuf_size_, cbuf_);
					p2p.recv(cbuf_, header.count);
					MPI_Send(cbuf_, header.count, MPI_BYTE,
						header.rank, 0, p2p.dmp_comm());
					break;
				case P2PTag::recv:
				case P2PTag::irecv:
					resize_buf(header.count, cbuf_size_, cbuf_);
					MPI_Recv(cbuf_, header.count, MPI_BYTE,
						header.rank, 0, p2p.dmp_comm(), &status);
					p2p.send(cbuf_, header.count);
					break;
				case P2PTag::allreduce_max_double:
					resize_buf(header.count, dbuf_size_, dbuf_);
					p2p.recv(dbuf_, header.count);
					MPI_Allreduce(dbuf_, dbuf_,
						header.count, MPI_DOUBLE, MPI_MAX, p2p.dmp_comm());
					p2p.send(dbuf_, header.count);
					break;
				case P2PTag::allreduce_sum_double:
					resize_buf(header.count, dbuf_size_, dbuf_);
					p2p.recv(dbuf_, header.count);
					MPI_Allreduce(dbuf_, dbuf_,
						header.count, MPI_DOUBLE, MPI_SUM, p2p.dmp_comm());
					p2p.send(dbuf_, header.count);
					break;
				case P2PTag::allreduce_sum_int:
					resize_buf(header.count, ibuf_size_, ibuf_);
					p2p.recv(ibuf_, header.count);
					MPI_Allreduce(ibuf_, ibuf_,
						header.count, MPI_INT, MPI_SUM, p2p.dmp_comm());
					p2p.send(ibuf_, header.count);
					break;
				case P2PTag::allgather_int:
					resize_buf(header.count*p2p.size(), ibuf_size_, ibuf_);
					p2p.recv(ibuf_, header.count);
					MPI_Allgather(ibuf_, header.count, MPI_INT,
						ibuf_, header.count, MPI_INT, p2p.dmp_comm());
					p2p.send(ibuf_, header.count*p2p.size());
					break;
				case P2PTag::allgather_long_long:
					resize_buf(header.count*p2p.size(), lbuf_size_, lbuf_);
					p2p.recv(lbuf_, header.count);
					MPI_Allgather(lbuf_, header.count, MPI_LONG_LONG,
						lbuf_, header.count, MPI_LONG_LONG, p2p.dmp_comm());
					p2p.send(lbuf_, header.count*p2p.size());
					break;
				case P2PTag::wait:
					// nothing needs to be done on this side
					break;
				case P2PTag::barrier:
					MPI_Barrier(p2p.dmp_comm());
					p2p.sync();
					break;
				case P2PTag::abort:
					int reason;
					p2p.recv(&reason, 1);
					MPI_Abort(p2p.dmp_comm(), reason);
					MPI_Abort(p2p.peer_comm(), reason);
					break;
				case P2PTag::end:
					relay = false;
					break;
			} // switch
		} // while
	} // DMPRelay::start

template<typename T>
void DMPRelay::resize_buf(int count, int & current_count, T *& buffer)
	{
		if(count > current_count) {
			T * tmp = buffer;
			buffer = new T[count];
			memcpy(buffer, tmp, current_count*sizeof(T));
			current_count = count;
			delete tmp;
		} // if
	} // DMPRelay::resize_cbuf

#endif // DMPRelay_hxx
