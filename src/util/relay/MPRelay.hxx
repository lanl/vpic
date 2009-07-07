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

//#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>
#include "P2PConnection.hxx"
#include "DMPConnection.hxx"
#include "MPData.hxx"
#include "FileIO.hxx"
#include "FileUtils.hxx"

/*!
	\class MPRelay MPRelay.h
	\brief  provides...
*/
class MPRelay
	{
	public:

		//! Constructor
		MPRelay() : file_dsc_(0) {}

		//! Destructor
		~MPRelay() {}

		void start();

	private:

		MPBuffer<char> cbuf_send_[max_buffers];
		MPBuffer<char> cbuf_recv_[max_buffers];

		MPBuffer<unsigned char> ucbuf_send_;
		MPBuffer<unsigned char> ucbuf_recv_;

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

		int32_t file_dsc_;
		std::map<int, FileIO> file_;
		MPBuffer<char, filename_size> filename_;
		MPBuffer<char, io_buffer_size> io_buffer_;
		MPBuffer<char, utils_buffer_size> utils_buffer_;

	}; // class MPRelay

void MPRelay::start()
	{
		P2PConnection & p2p = P2PConnection::instance();
		DMPConnection & dmp = DMPConnection::instance();

		bool relay(true);
		int64_t filesize;
		int64_t foffset;
		int32_t fwhence;
		int32_t utils_return;
		int32_t fclose_return;
		double wtime;
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
					std::cerr << "isend request " << request;
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
					std::cerr << "irecv request " << request;
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
/*
					std::cout << "waiting recv requested" << std::endl;
*/
					if(dmp_recv_request_[request.id].state == pending) {
/*
					std::cerr << "waiting recv " << request;
*/

						// block on irecv if it hasn't finished
						dmp.wait_recv(dmp_recv_request_[request.id]);

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
/*
					std::cout << "waiting send requested" << std::endl;
*/
					if(dmp_send_request_[request.id].state == pending) {
/*
					std::cerr << "waiting send " << request;
*/

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

				case P2PTag::wtime:
					//std::cout << "wtime requested from dmp: count " <<
					//	request.count << " tag " << request.tag << std::endl;
					wtime = dmp.wtime();
					p2p.send(&wtime, request.count, request.tag);
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

				case P2PTag::gather_uc:
					// resize buffers if necessary
					ucbuf_send_.resize(request.count);
					ucbuf_recv_.resize(request.count*dmp.global_size());

					// blocking receive from point-to-point peer
					p2p.recv(ucbuf_send_.data(), request.count,
						request.tag, request.id);

					/*
					std::cerr << "#### request.count: " << request.count
						<< std::endl <<
						"global_size: " << dmp.global_size() << std::endl;;
					*/

					// call all gather with dmp peers
					dmp.gather(ucbuf_send_.data(), ucbuf_recv_.data(),
						request.count);

					if(dmp.global_id() == 0) {
						// blocking send result back to point-to-point peer
						p2p.send(ucbuf_recv_.data(),
							request.count*dmp.global_size(), request.tag);
					} // if
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

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" io_open_read " << filename_.data() <<
						" filedsc " << file_dsc_ << std::endl;
					*/

					// make sure that this id is not in use
					assert(file_.find(file_dsc_) == file_.end());

					file_[file_dsc_].open(filename_.data(), io_read);
					filesize = file_[file_dsc_].size();
					p2p.send(&file_dsc_, 1, request.tag);
					p2p.send(&filesize, 1, request.tag);

					// increment id
					++file_dsc_;
					break;

				case P2PTag::io_open_read_write:
					p2p.recv(filename_.data(), request.count,
						request.tag, request.id);

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" io_open_read_write " << filename_.data() << std::endl;
					*/

					// make sure that this id is not in use
					assert(file_.find(file_dsc_) == file_.end());

					file_[file_dsc_].open(filename_.data(), io_read_write);
					filesize = file_[file_dsc_].size();
					p2p.send(&file_dsc_, 1, request.tag);
					p2p.send(&filesize, 1, request.tag);

					// increment id
					++file_dsc_;
					break;

				case P2PTag::io_open_write:
					p2p.recv(filename_.data(), request.count,
						request.tag, request.id);

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" io_open_write " << filename_.data() <<
						" id " << file_dsc_ << std::endl;
					*/

					// make sure that this id is not in use
					assert(file_.find(file_dsc_) == file_.end());

					file_[file_dsc_].open(filename_.data(), io_write);
					p2p.send(&file_dsc_, 1, request.tag);

					// increment id
					++file_dsc_;
					break;

				case P2PTag::io_open_write_read:
					p2p.recv(filename_.data(), request.count,
						request.tag, request.id);

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" io_open_write_read " << filename_.data() << std::endl;
					*/

					// make sure that this id is not in use
					assert(file_.find(file_dsc_) == file_.end());

					file_[file_dsc_].open(filename_.data(), io_write_read);
					p2p.send(&file_dsc_, 1, request.tag);

					// increment id
					++file_dsc_;
					break;

				case P2PTag::io_open_append:
					p2p.recv(filename_.data(), request.count,
						request.tag, request.id);

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" io_open_append " << filename_.data() << std::endl;
					*/

					// make sure that this id is not in use
					assert(file_.find(file_dsc_) == file_.end());

					file_[file_dsc_].open(filename_.data(), io_append);
					p2p.send(&file_dsc_, 1, request.tag);

					// increment id
					++file_dsc_;
					break;

				case P2PTag::io_open_append_read:
					p2p.recv(filename_.data(), request.count,
						request.tag, request.id);

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" io_open_append_read " <<
						filename_.data() << std::endl;
					*/

					// make sure that this id is not in use
					assert(file_.find(file_dsc_) == file_.end());

					file_[file_dsc_].open(filename_.data(), io_append_read);
					p2p.send(&file_dsc_, 1, request.tag);

					// increment id
					++file_dsc_;
					break;

				case P2PTag::io_read:
					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" requesting io_read of " << request.count <<
						" bytes" << std::endl;
					*/

					// make sure that this id exists
					assert(file_.find(request.id) != file_.end());

					file_[request.id].read(io_buffer_.data(), request.count);
					p2p.send(io_buffer_.data(), request.count,
						request.tag);
					break;

				case P2PTag::io_write:
					p2p.recv(io_buffer_.data(), request.count,
						request.tag, request.id);

					// make sure that this id exists
					assert(file_.find(request.id) != file_.end());

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" writing " << request.count
						<< " to file id " << request.id << std::endl;
					*/

					file_[request.id].write(io_buffer_.data(), request.count);
					break;

				case P2PTag::io_seek:
					// make sure that this id exists
					assert(file_.find(request.id) != file_.end());

					p2p.recv(&foffset, 1, request.tag, request.id);
					p2p.recv(&fwhence, 1, request.tag, request.id);

					/*
					if(p2p.global_id() == 0) {
					std::cerr << "seek " << foffset << " " <<
						fwhence << std::endl;
					} // if
					*/

					file_[request.id].seek(foffset, fwhence);
					break;

				case P2PTag::io_tell:
					// make sure that this id exists
					assert(file_.find(request.id) != file_.end());

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" calling tell on " << filename_.data() << std::endl;
					*/

					foffset = file_[request.id].tell();

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" tell value " << foffset << std::endl;
					*/

					p2p.send(&foffset, 1, request.tag);
					break;

				case P2PTag::io_rewind:
					// make sure that this id exists
					assert(file_.find(request.id) != file_.end());

					file_[request.id].rewind();
					break;

				case P2PTag::io_size:
					// make sure that this id exists
					assert(file_.find(request.id) != file_.end());

					filesize = file_[request.id].size();

					p2p.send(&filesize, 1, request.tag);
					break;

				case P2PTag::io_close:

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" io_close id " << request.id << std::endl;
					*/

					// make sure that this id exists
					assert(file_.find(request.id) != file_.end());

					/*
					std::cerr << "rank: " << p2p.global_id() <<
						" writing " << request.count << " to file " <<
						std::endl;
					*/

					fclose_return = file_[request.id].close();
					p2p.send(&fclose_return, 1, request.tag);
					file_.erase(request.id);
					break;

				case P2PTag::utils_mkdir:
					p2p.recv(utils_buffer_.data(), request.count,
						request.tag, request.id);
					utils_return =
						FileUtils::makeDirectory(utils_buffer_.data());
					p2p.send(&utils_return, 1, request.tag);
					break;

				case P2PTag::utils_getcwd:
					utils_return = FileUtils::getCurrentWorkingDirectory(
							utils_buffer_.data(), request.count);

					p2p.send(io_buffer_.data(), request.count,
						request.tag);
					p2p.send(&utils_return, 1, request.tag);
					break;

				case P2PTag::end:
					relay = false;
					break;

				case P2PTag::pending:
					// if file map is empty, reset file_dsc_
					if(file_.size() == 0 && file_dsc_ > 0) {
						file_dsc_ = 0;
						/*
						std::cerr << "rank: " << p2p.global_id() <<
							" resetting map " << file_.size() << std::endl;
						*/
					}


					// check for finished send communications
					for(std::vector<int>::iterator ita =
						pending_dmp_send_.begin();
						ita != pending_dmp_send_.end();) {

						// test for completion
						dmp.test_send(dmp_send_request_[*ita]);

						if(dmp_send_request_[*ita].state == complete) {
/*
							std::cerr << "dmp send request complete" <<
								dmp_send_request_[*ita];
*/
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
/*
							std::cerr << "recv request complete" <<
								dmp_recv_request_[*ita];
*/
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
					break;

				default:
					std::cerr << "ERROR: Unknown Tag" << std::endl;
					exit(1);
					break;

			} // switch
		} // while
	} // MPRelay::start

#endif // MPRelay_hxx
