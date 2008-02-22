/*
	Definition of MPRelay class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef MPRelay_hxx
#define MPRelay_hxx

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <P2PConnection.hxx>
#include <DMPConnection.hxx>
#include <MPData.hxx>
#include <FileIO.hxx>

// constant to determine how often to
// check pending messages
const size_t HANDLE_PENDING(20);

/*!
	\class MPRelay MPRelay.h
	\brief  provides...
*/
class MPRelay
	{
	public:

		//! Constructor
		MPRelay() {}

		//! Destructor
		~MPRelay() {}

		void start();

	private:

		MPBuffer<char> cbuf_send_[max_buffers];
		MPBuffer<char> cbuf_recv_[max_buffers];

		MPBuffer<int> ibuf_send_;
		MPBuffer<int> ibuf_recv_;

		MPBuffer<int64_t> lbuf_send_;
		MPBuffer<int64_t> lbuf_recv_;

		MPBuffer<double> dbuf_send_;
		MPBuffer<double> dbuf_recv_;

		MPRequest_T<MP_HOST> dmp_recv_request_[max_buffers];
		MPRequest_T<MP_HOST> dmp_send_request_[max_buffers];
		MPRequest_T<MP_HOST> p2p_send_request_[max_buffers];

		std::vector<int> pending_dmp_recv_;
		std::vector<int> pending_dmp_send_;
		std::vector<int> pending_p2p_send_;

		FileIO fileIO_;
		MPBuffer<char, filename_size> filename_;
		MPBuffer<char, io_buffer_size> io_buffer_;

	}; // class MPRelay

void MPRelay::start()
	{
		P2PConnection & p2p = P2PConnection::instance();
		DMPConnection & dmp = DMPConnection::instance();

		bool relay(true);
		int filesize;
		size_t pending_count(0);
		MPRequest_T<MP_HOST> request;

		while(relay) {
			switch(p2p.poll(request)) {

				case P2PTag::send:
/*
					std::cout << "send request" << std::endl;
					std::cout << request;
*/

					// resize buffer if necessary
					cbuf_send_[request.id].resize(request.count);

					// blocking receive from point-to-point peer
					p2p.recv(cbuf_send_[request.id].data(), request.count,
						request.tag, request.id);

					// blocking send to dmp peer
					dmp.send(cbuf_send_[request.id].data(), request.count,
						request.peer, request.tag);
					break;

				case P2PTag::recv:
/*
					std::cout << "recv request" << std::endl;
					std::cout << request;
*/

					// resize buffer if necessary
					cbuf_recv_[request.id].resize(request.count);

					// blocking receive from dmp peer
					dmp.recv(cbuf_recv_[request.id].data(), request.count,
						request.peer, request.tag, request.id);

					// blocking send to point-to-point peer
					p2p.send(cbuf_recv_[request.id].data(),
						request.count, request.tag);
					break;

					// this case is included for completeness and
					// in case the host process is ever required to
					// do real work.
				case P2PTag::isend:
/*
					std::cout << "isend request" << std::endl;
					std::cout << request;
*/

					// resize buffer if necessary
					cbuf_send_[request.id].resize(request.count);

					// save request
					dmp_send_request_[request.id] = request;

					// blocking receive from point-to-point peer
					p2p.recv(cbuf_send_[request.id].data(), request.count,
						request.tag, request.id);

					// non-blocking send to dmp peer
					dmp.isend(cbuf_send_[request.id].data(), request.count,
						request.peer, request.tag, request.id);

					// keep track of pending sends
					dmp_send_request_[request.id].state = pending;
					pending_dmp_send_.push_back(request.id);
					break;

				case P2PTag::irecv:
/*
					std::cout << "irecv request" << std::endl;
					std::cout << request;
*/

					// resize buffer if necessary
					cbuf_recv_[request.id].resize(request.count);

					// save request
					dmp_recv_request_[request.id] = request;

					// non-blocking receive to dmp peer
					dmp.irecv(cbuf_recv_[request.id].data(), request.count,
						request.peer, request.tag, request.id);

					// keep track of pending receives
					dmp_recv_request_[request.id].state = pending;
					pending_dmp_recv_.push_back(request.id);
					break;

				case P2PTag::wait_recv:
					if(dmp_recv_request_[request.id].state == pending) {
						// block on irecv if it hasn't finished
						dmp.wait_recv(dmp_recv_request_[request.id]);

						/*
						p2p.send(
							cbuf_recv_[dmp_recv_request_[request.id].id].data(),
							dmp_recv_request_[request.id].count,
							dmp_recv_request_[request.id].tag);
						*/
						if(p2p_send_request_[request.id].state == pending) {
							p2p.wait_send(request.id);
						} // if

						p2p_send_request_[request.id] =
							dmp_recv_request_[request.id];

						p2p.isend(
							cbuf_recv_[p2p_send_request_[request.id].id].data(),
							p2p_send_request_[request.id].count,
							p2p_send_request_[request.id].tag,
							p2p_send_request_[request.id].id);
						p2p_send_request_[request.id].state = pending;
						pending_p2p_send_.push_back(request.id);

						// remove from pending dmp receives
						std::vector<int>::iterator ita =
							std::find(pending_dmp_recv_.begin(),
							pending_dmp_recv_.end(), request.id);
						pending_dmp_recv_.erase(ita);
					}
					break;

				case P2PTag::wait_send:
					if(dmp_send_request_[request.id].state == pending) {
						// block on isend if it hasn't finished
						dmp.wait_send(dmp_send_request_[request.id]);

						// remove from pending sends
						std::vector<int>::iterator ita =
							std::find(pending_dmp_send_.begin(),
							pending_dmp_send_.end(), request.id);
						pending_dmp_send_.erase(ita);
					}
					break;

				case P2PTag::barrier:
					dmp.barrier();
					p2p.barrier();
					break;

				case P2PTag::allreduce_max_double:
					// resize buffers if necessary
					dbuf_send_.resize(request.count);
					dbuf_recv_.resize(request.count);

					// blocking receive from point-to-point peer
					p2p.recv(dbuf_send_.data(), request.count,
						request.tag, request.id);

					// call all reduce max with dmp peers
					dmp.allreduce_max(dbuf_send_.data(), dbuf_recv_.data(),
						request.count);

					// blocking send result back to point-to-point peer
					p2p.send(dbuf_recv_.data(), request.count, request.tag);
					break;

				case P2PTag::allreduce_sum_double:
					// resize buffers if necessary
					dbuf_send_.resize(request.count);
					dbuf_recv_.resize(request.count);

					// blocking receive from point-to-point peer
					p2p.recv(dbuf_send_.data(), request.count,
						request.tag, request.id);

					// call all reduce max with dmp peers
					dmp.allreduce_sum(dbuf_send_.data(), dbuf_recv_.data(),
						request.count);

					// blocking send result back to point-to-point peer
					p2p.send(dbuf_recv_.data(), request.count, request.tag);
					break;

				case P2PTag::allreduce_sum_int:
					// resize buffers if necessary
					ibuf_send_.resize(request.count);
					ibuf_recv_.resize(request.count);

					// blocking receive from point-to-point peer
					p2p.recv(ibuf_send_.data(), request.count,
						request.tag, request.id);

					// call all reduce max with dmp peers
					dmp.allreduce_sum(ibuf_send_.data(), ibuf_recv_.data(),
						request.count);

					// blocking send result back to point-to-point peer
					p2p.send(ibuf_recv_.data(), request.count, request.tag);
					break;

				case P2PTag::allgather_int:
					// resize buffers if necessary
					ibuf_send_.resize(request.count);
					ibuf_recv_.resize(request.count*dmp.global_size());

					// blocking receive from point-to-point peer
					p2p.recv(ibuf_send_.data(), request.count,
						request.tag, request.id);

					// call all gather with dmp peers
					dmp.allgather(ibuf_send_.data(), ibuf_recv_.data(),
						request.count);

					// blocking send result back to point-to-point peer
					p2p.send(ibuf_recv_.data(),
						request.count*dmp.global_size(), request.tag);
					break;

				case P2PTag::allgather_int64:
					// resize buffers if necessary
					lbuf_send_.resize(request.count);
					lbuf_recv_.resize(request.count*dmp.global_size());

					// blocking receive from point-to-point peer
					p2p.recv(lbuf_send_.data(), request.count,
						request.tag, request.id);

					// call all gather with dmp peers
					dmp.allgather(lbuf_send_.data(), lbuf_recv_.data(),
						request.count);

					// blocking send result back to point-to-point peer
					p2p.send(lbuf_recv_.data(),
						request.count*dmp.global_size(), request.tag);
					break;

				case P2PTag::abort:
					int reason;
					p2p.recv(&reason, 1, request.tag, request.id);
					dmp.abort(reason);
					p2p.abort(reason);
					break;

				case P2PTag::io_open_read:
					p2p.recv(filename_.data(), request.count,
						request.tag, request.id);
					fileIO_.open(filename_.data(), io_read);
					filesize = fileIO_.size();
					p2p.send(&filesize, 1, request.tag);
					break;

				case P2PTag::io_open_write:
					p2p.recv(filename_.data(), request.count,
						request.tag, request.id);
					fileIO_.open(filename_.data(), io_write);
					break;

				case P2PTag::io_open_write_append:
					p2p.recv(filename_.data(), request.count,
						request.tag, request.id);
					fileIO_.open(filename_.data(), io_write_append);
					break;

				case P2PTag::io_read:
					fileIO_.read(io_buffer_.data(), request.count);
					p2p.send(io_buffer_.data(), request.count, request.tag);
					break;

				case P2PTag::io_write:
					p2p.recv(io_buffer_.data(), request.count,
						request.tag, request.id);
					fileIO_.write(io_buffer_.data(), request.count);
					break;

				case P2PTag::io_close:
					fileIO_.close();
					break;

				case P2PTag::end:
					relay = false;
					break;

				case P2PTag::pending:
					if(pending_count++ == HANDLE_PENDING) {
						// reset counter
						pending_count = 0;

						// check for finished send communications
						for(std::vector<int>::iterator ita =
							pending_dmp_send_.begin();
							ita != pending_dmp_send_.end();) {

							// test for completion
							dmp.test_send(dmp_send_request_[*ita]);

							if(dmp_send_request_[*ita].state == complete) {
								ita = pending_dmp_send_.erase(ita);
							}
							else {
								++ita;
							} // if
						} // for

						// check for finished recv communications
						for(std::vector<int>::iterator ita =
							pending_dmp_recv_.begin();
							ita != pending_dmp_recv_.end();) {

							// test for completion
							dmp.test_recv(dmp_recv_request_[*ita]);

							if(dmp_recv_request_[*ita].state == complete) {
								// forward to accelerator
								p2p.send(
									cbuf_recv_[dmp_recv_request_[*ita].id].data(),
									dmp_recv_request_[*ita].count,
									dmp_recv_request_[*ita].tag);

								// remove from pending
								ita = pending_dmp_recv_.erase(ita);
							}
							else {
								++ita;
							} // if
						} // for

						for(std::vector<int>::iterator ita =
							pending_p2p_send_.begin();
							ita != pending_p2p_send_.end();) {

							// test for completion
							p2p.test_send(p2p_send_request_[*ita]);

							if(p2p_send_request_[*ita].state == complete) {
								ita = pending_p2p_send_.erase(ita);
							}
							else {
								++ita;
							} // if
						} // for
					} // if
					break;

				default:
					std::cerr << "ERROR: Unknown Tag" << std::endl;
					exit(1);
					break;

			} // switch
		} // while
	} // MPRelay::start

#endif // MPRelay_hxx
