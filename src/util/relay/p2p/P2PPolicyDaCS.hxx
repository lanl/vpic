/*
	Definition of P2PPolicyDaCS class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef P2PPolicyDaCS_hxx
#define P2PPolicyDaCS_hxx

#include <iostream>
#include <assert.h>
#include <dacs.h>
#include <ConnectionManager.hxx>
#include <MPData.hxx>
#include <P2PTag.hxx>
#include <Type2DaCSSwapType.hxx>
#include <DaCSUtils.h>

/*!
	\struct P2PPolicyDaCS P2PPolicyDaCS.h
	\brief  provides...
*/
template<int ROLE> class P2PPolicyDaCS
	{
	public:

		// topology information
		inline int global_id()
			{ return ConnectionManager::instance().global_id(); }
		inline int global_size()
			{ return ConnectionManager::instance().global_size(); }

		// host side poll
		inline int poll(MPRequest_T<ROLE> & request);

		// accelerator side request
		inline int post(int p2ptag);
		inline int post(MPRequest_T<ROLE> & request);

		// send and recv methods
		template<typename T> int send(T * buffer, int count, int tag);
		template<typename T> int recv(T * buffer, int count,
			int tag, int id);
		template<typename T> int isend(T * buffer, int count,
			int tag, int id);
		template<typename T> int irecv(T * buffer, int count,
			int tag, int id);

		inline int wait_send(int id);
		inline int wait_recv(int id);
		inline int barrier();

		// return elements handled in status
		// dummy is a hack for broken ppu-g++
		template<typename T> int get_count(int id, int & count,
			T * dummy = NULL);

		double wtime();

		int abort(int reason);

	protected:

		P2PPolicyDaCS();
		~P2PPolicyDaCS();

	private:

		bool pending_;
		DACS_ERR_T errcode_;
		dacs_tag_t request_tag_;
		dacs_tag_t blocking_send_tag_;
		dacs_tag_t blocking_recv_tag_;

		dacs_tag_t send_tag_[max_buffers];
		dacs_tag_t recv_tag_[max_buffers];

	}; // class P2PPolicyDaCS

template<int ROLE>
P2PPolicyDaCS<ROLE>::P2PPolicyDaCS()
	:  pending_(false)
	{
		ConnectionManager & mgr = ConnectionManager::instance();

		// register request tag
		errcode_ = dacs_tag_reserve(&request_tag_,
			mgr.peer_de(), mgr.peer_pid());	
		process_dacs_errcode(errcode_);

		// register blocking send and recv tags
		errcode_ = dacs_tag_reserve(&blocking_send_tag_,
			mgr.peer_de(), mgr.peer_pid());	
		process_dacs_errcode(errcode_);

		errcode_ = dacs_tag_reserve(&blocking_recv_tag_,
			mgr.peer_de(), mgr.peer_pid());	
		process_dacs_errcode(errcode_);

		// register all of the send and recv tags
		for(size_t i(0); i<max_buffers; i++) {
			errcode_ = dacs_tag_reserve(&send_tag_[i],
				mgr.peer_de(), mgr.peer_pid());	
			process_dacs_errcode(errcode_);

			errcode_ = dacs_tag_reserve(&recv_tag_[i],
				mgr.peer_de(), mgr.peer_pid());	
			process_dacs_errcode(errcode_);
		} // for
	} // P2PPolicyDaCS<>::P2PPolicyDaCS

template<int ROLE>
P2PPolicyDaCS<ROLE>::~P2PPolicyDaCS()
	{
		ConnectionManager & mgr = ConnectionManager::instance();

		// release request tag
		errcode_ = dacs_tag_release(&request_tag_);
		process_dacs_errcode(errcode_);

		// release blocking send and recv tags
		errcode_ = dacs_tag_release(&blocking_send_tag_);
		process_dacs_errcode(errcode_);

		errcode_ = dacs_tag_release(&blocking_recv_tag_);
		process_dacs_errcode(errcode_);

		// release all of the send and recv tags
		for(size_t i(0); i<max_buffers; i++) {
			errcode_ = dacs_tag_release(&send_tag_[i]);
			process_dacs_errcode(errcode_);

			errcode_ = dacs_tag_release(&recv_tag_[i]);
			process_dacs_errcode(errcode_);
		} // for
	} // P2PPolicyDaCS<>::P2PPolicyDaCS

template<> inline
int P2PPolicyDaCS<MP_HOST>::poll(MPRequest_T<MP_HOST> & request)
	{
		ConnectionManager & mgr = ConnectionManager::instance();

		if(pending_) {
			// test for message completion
			errcode_ = dacs_test(request_tag_);

			switch(errcode_) {
				case DACS_ERR_TEST_READY:
					pending_ = false;
					return request.p2ptag;
				case DACS_ERR_TEST_BUSY:
					return P2PTag::pending;
				default:
					process_dacs_errcode(errcode_);
			} // switch
		}
		else {
			// initiated new recv operation for next request
			errcode_ = dacs_recv(&request, request_count(), mgr.peer_de(),
				mgr.peer_pid(), request_tag_, DACS_BYTE_SWAP_WORD);
			std::cerr << "host called dacs_recv(in relay) " << errcode_ << std::endl;
			process_dacs_errcode(errcode_);

			pending_ = true;
			return P2PTag::pending;
		} // if
	} // P2PPolicyDaCS<>::poll

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::post(int p2ptag)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		MPRequest_T<MP_ACCEL> request(p2ptag);

		errcode_ = dacs_send(mgr.peer_de(), mgr.peer_pid(), &request,
			request_count(), request_tag_, DACS_BYTE_SWAP_WORD);
		std::cerr << "accelerator called dacs_send(in relay) " << errcode_ << std::endl;
		process_dacs_errcode(errcode_);

		errcode_ = dacs_wait(request_tag_);
		process_dacs_errcode(errcode_);
	} // MPICommunicatorPolicy<>::request

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::post(MPRequest_T<MP_ACCEL> & request)
	{
		ConnectionManager & mgr = ConnectionManager::instance();

		errcode_ = dacs_send(mgr.peer_de(), mgr.peer_pid(), &request,
			request_count(), request_tag_, DACS_BYTE_SWAP_WORD);
		std::cerr << "accelerator called dacs_send(in relay) " << errcode_ << std::endl;
		process_dacs_errcode(errcode_);

		errcode_ = dacs_wait(request_tag_);
		process_dacs_errcode(errcode_);
	} // MPICommunicatorPolicy<>::request

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::send(T * buffer, int count, int tag)
	{
		std::cout << "requesting isend of " << count <<
			" elements" << std::endl;
		output_swap_type(Type2DaCSSwapType<T>::type());

		ConnectionManager & mgr = ConnectionManager::instance();
		errcode_ = dacs_send(mgr.peer_de(), mgr.peer_pid(), buffer,
			count*sizeof(T), blocking_send_tag_, Type2DaCSSwapType<T>::type());
		process_dacs_errcode(errcode_);

		errcode_ = dacs_wait(blocking_send_tag_);
		process_dacs_errcode(errcode_);

		return 0;
	} // MPICommunicatorPolicy<>::send

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::recv(T * buffer, int count, int tag, int id)
	{
		std::cout << "requesting irecv of " << count <<
			" elements" << std::endl;
		output_swap_type(Type2DaCSSwapType<T>::type());

		ConnectionManager & mgr = ConnectionManager::instance();
		errcode_ = dacs_recv(buffer, count*sizeof(T), mgr.peer_de(),
			mgr.peer_pid(), blocking_recv_tag_, Type2DaCSSwapType<T>::type());
		process_dacs_errcode(errcode_);

		errcode_ = dacs_wait(blocking_recv_tag_);
		process_dacs_errcode(errcode_);

		return 0;
	} // MPICommunicatorPolicy<>::recv

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::isend(T * buffer, int count, int tag, int id)
	{
		std::cout << "requesting isend of " << count <<
			" elements" << std::endl;
		output_swap_type(Type2DaCSSwapType<T>::type());

		ConnectionManager & mgr = ConnectionManager::instance();
		errcode_ = dacs_send(mgr.peer_de(), mgr.peer_pid(), buffer,
			count*sizeof(T), send_tag_[id], Type2DaCSSwapType<T>::type());
		process_dacs_errcode(errcode_);

		return 0;
	} // MPICommunicatorPolicy<>::isend

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::irecv(T * buffer, int count, int tag, int id)
	{
		std::cout << "requesting irecv of " << count <<
			" elements" << std::endl;
		output_swap_type(Type2DaCSSwapType<T>::type());
		ConnectionManager & mgr = ConnectionManager::instance();
		errcode_ = dacs_recv(buffer, count*sizeof(T), mgr.peer_de(),
			mgr.peer_pid(), recv_tag_[id], Type2DaCSSwapType<T>::type());
		process_dacs_errcode(errcode_);

		return 0;
	} // MPICommunicatorPolicy<>::irecv

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::wait_send(int id)
	{
		errcode_ = dacs_wait(send_tag_[id]);
		process_dacs_errcode(errcode_);

		return 0;
	} // P2PPolicyDaCS<>::wait

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::wait_recv(int id)
	{
		errcode_ = dacs_wait(recv_tag_[id]);
		process_dacs_errcode(errcode_);

		return 0;
	} // P2PPolicyDaCS<>::wait

template<int ROLE> inline
int P2PPolicyDaCS<ROLE>::barrier()
	{
		// hack until barriers are implemented in DaCS
		sleep(10);
		return 0;
	} // P2PPolicyDaCS<>::wait

template<>
template<typename T>
int P2PPolicyDaCS<MP_ACCEL>::get_count(int id, int & count, T * dummy)
	{
	} // P2PPolicyDaCS<>::get_count

template<int ROLE> inline
double P2PPolicyDaCS<ROLE>::wtime()
	{
	} // P2PPolicyDaCS<>::wtime

template<int ROLE>
int P2PPolicyDaCS<ROLE>::abort(int reason)
	{
	} // P2PPolicyDaCS<>::wait

#endif // P2PPolicyDaCS_hxx
