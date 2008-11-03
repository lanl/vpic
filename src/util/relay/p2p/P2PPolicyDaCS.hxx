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
#include "ConnectionManager.hxx"
#include "MPData.hxx"
#include "P2PTag.hxx"
#include "Type2DaCSSwapType.hxx"
#include "DaCSUtils.h"

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

		inline int role()
			{ return ROLE; }

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

		inline int test_send(MPRequest_T<MP_HOST> & request);
		inline int wait_send(int id);
		inline int wait_recv(int id);
		inline int barrier();

		// return elements handled in status
		// dummy is a hack for broken ppu-g++
		template<typename T> int get_count(int id, int & count,
			T * dummy = NULL);

		int abort(int reason);

	protected:

		P2PPolicyDaCS();
		~P2PPolicyDaCS();

	private:

		bool pending_;
		DACS_ERR_T errcode_;
		dacs_wid_t request_wid_;
		dacs_wid_t blocking_send_wid_;
		dacs_wid_t blocking_recv_wid_;

		dacs_wid_t send_wid_[max_buffers];
		dacs_wid_t recv_wid_[max_buffers];
		size_t recv_count_[max_buffers];

	}; // class P2PPolicyDaCS

template<int ROLE>
P2PPolicyDaCS<ROLE>::P2PPolicyDaCS()
	:  pending_(false)
	{
		//ConnectionManager & mgr = ConnectionManager::instance();

		// register request wid
		errcode_ = dacs_wid_reserve(&request_wid_);	
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// register blocking send and recv wids
		errcode_ = dacs_wid_reserve(&blocking_send_wid_);	
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wid_reserve(&blocking_recv_wid_);	
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// register all of the send and recv wids
		for(size_t i(0); i<max_buffers; i++) {
			errcode_ = dacs_wid_reserve(&send_wid_[i]);	
			process_dacs_errcode(errcode_, __FILE__, __LINE__);

			errcode_ = dacs_wid_reserve(&recv_wid_[i]);
			process_dacs_errcode(errcode_, __FILE__, __LINE__);

			recv_count_[i] = 0;
		} // for
	} // P2PPolicyDaCS<>::P2PPolicyDaCS

template<int ROLE>
P2PPolicyDaCS<ROLE>::~P2PPolicyDaCS()
	{
		//ConnectionManager & mgr = ConnectionManager::instance();

		// release request wid
		errcode_ = dacs_wid_release(&request_wid_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// release blocking send and recv wids
		errcode_ = dacs_wid_release(&blocking_send_wid_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wid_release(&blocking_recv_wid_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// release all of the send and recv wids
		for(size_t i(0); i<max_buffers; i++) {
			errcode_ = dacs_wid_release(&send_wid_[i]);
			process_dacs_errcode(errcode_, __FILE__, __LINE__);

			errcode_ = dacs_wid_release(&recv_wid_[i]);
			process_dacs_errcode(errcode_, __FILE__, __LINE__);
		} // for
	} // P2PPolicyDaCS<>::P2PPolicyDaCS

template<> inline
int P2PPolicyDaCS<MP_HOST>::poll(MPRequest_T<MP_HOST> & request)
	{
		ConnectionManager & mgr = ConnectionManager::instance();

		if(pending_) {
			// test for message completion
			errcode_ = dacs_test(request_wid_);

			switch(errcode_) {
				case DACS_WID_READY:
					pending_ = false;
					return request.p2ptag;
				case DACS_WID_BUSY:
					return P2PTag::pending;
				default:
					process_dacs_errcode(errcode_, __FILE__, __LINE__);
			} // switch
		}
		else {
			// initiated new recv operation for next request
			errcode_ = dacs_recv(&request, request_count(),
				mgr.peer_de(), mgr.peer_pid(), P2PTag::request,
				request_wid_, DACS_BYTE_SWAP_WORD);
#if 0
			std::cerr << "host called dacs_recv(in relay) " <<
				errcode_ << std::endl;
#endif
			process_dacs_errcode(errcode_, __FILE__, __LINE__);

			pending_ = true;
			return P2PTag::pending;
		} // if
	} // P2PPolicyDaCS<>::poll

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::post(int p2ptag)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		MPRequest_T<MP_ACCEL> request(p2ptag);

		errcode_ = dacs_send(&request, request_count(),
			mgr.peer_de(), mgr.peer_pid(), P2PTag::request,
			request_wid_, DACS_BYTE_SWAP_WORD);
#if 0
		std::cerr << "accelerator called dacs_send(in relay) " <<
			errcode_ << std::endl;
#endif
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wait(request_wid_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		return 0;
	} // MPICommunicatorPolicy<>::request

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::post(MPRequest_T<MP_ACCEL> & request)
	{
		ConnectionManager & mgr = ConnectionManager::instance();

		errcode_ = dacs_send(&request, request_count(),
			mgr.peer_de(), mgr.peer_pid(), P2PTag::request,
			request_wid_, DACS_BYTE_SWAP_WORD);
#if 0
		std::cerr << "accelerator called dacs_send(in relay) " <<
			errcode_ << std::endl;
#endif
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wait(request_wid_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		return 0;
	} // MPICommunicatorPolicy<>::request

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::send(T * buffer, int count, int tag)
	{
#if 0
		std::cout << "requesting isend of " << count <<
			" elements" << std::endl;
		output_swap_type(Type2DaCSSwapType<T>::type());
#endif

		ConnectionManager & mgr = ConnectionManager::instance();
		errcode_ = dacs_send(buffer, count*sizeof(T),
			mgr.peer_de(), mgr.peer_pid(), tag, blocking_send_wid_,
			Type2DaCSSwapType<T>::type());
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wait(blocking_send_wid_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		return 0;
	} // MPICommunicatorPolicy<>::send

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::recv(T * buffer, int count, int tag, int id)
	{
#if 0
		std::cout << "requesting irecv of " << count <<
			" elements" << std::endl;
		output_swap_type(Type2DaCSSwapType<T>::type());
#endif

		ConnectionManager & mgr = ConnectionManager::instance();
		errcode_ = dacs_recv(buffer, count*sizeof(T),
			mgr.peer_de(), mgr.peer_pid(), tag, blocking_recv_wid_,
			Type2DaCSSwapType<T>::type());
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wait(blocking_recv_wid_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		return 0;
	} // MPICommunicatorPolicy<>::recv

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::isend(T * buffer, int count, int tag, int id)
	{
#if 0
		std::cout << "requesting isend of " << count <<
			" elements" << std::endl;
		output_swap_type(Type2DaCSSwapType<T>::type());
#endif

		ConnectionManager & mgr = ConnectionManager::instance();
		errcode_ = dacs_send(buffer, count*sizeof(T),
			mgr.peer_de(), mgr.peer_pid(), tag, send_wid_[id],
			Type2DaCSSwapType<T>::type());
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		return 0;
	} // MPICommunicatorPolicy<>::isend

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::irecv(T * buffer, int count, int tag, int id)
	{
#if 0
		std::cout << "requesting irecv of " << count <<
			" elements" << std::endl;
		output_swap_type(Type2DaCSSwapType<T>::type());
#endif
		ConnectionManager & mgr = ConnectionManager::instance();
		errcode_ = dacs_recv(buffer, count*sizeof(T),
			mgr.peer_de(), mgr.peer_pid(), tag, recv_wid_[id],
			Type2DaCSSwapType<T>::type());
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		recv_count_[id] = count;

		return 0;
	} // MPICommunicatorPolicy<>::irecv

template<> inline
int P2PPolicyDaCS<MP_HOST>::test_send(MPRequest_T<MP_HOST> & request)
	{
		errcode_ = dacs_test(send_wid_[request.id]);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		if(errcode_ == DACS_WID_READY) {
			request.state = complete;	
		} // if

		return 0;
	} // P2PPolicyDaCS<>::test_send

template<int ROLE> inline
int P2PPolicyDaCS<ROLE>::wait_send(int id)
	{
		errcode_ = dacs_wait(send_wid_[id]);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		return 0;
	} // P2PPolicyDaCS<>::wait

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::wait_recv(int id)
	{
		errcode_ = dacs_wait(recv_wid_[id]);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		return 0;
	} // P2PPolicyDaCS<>::wait

template<int ROLE> inline
int P2PPolicyDaCS<ROLE>::barrier()
	{
		ConnectionManager & mgr = ConnectionManager::instance();

		//std::cout << "dacs_barrier_wait called" << std::endl;
		errcode_ = dacs_barrier_wait(mgr.group());
		process_dacs_errcode(errcode_, __FILE__, __LINE__);
		//std::cout << "dacs_barrier_wait passed" << std::endl;

		return 0;
	} // P2PPolicyDaCS<>::wait

template<>
template<typename T>
int P2PPolicyDaCS<MP_ACCEL>::get_count(int id, int & count, T * dummy)
	{
		// FIXME: this should probably just be removed from
		// the mp interface
		count = recv_count_[id];
		return 0;
	} // P2PPolicyDaCS<>::get_count

template<int ROLE>
int P2PPolicyDaCS<ROLE>::abort(int reason)
	{
		// FIXME
		exit(reason);
		return 0;
	} // P2PPolicyDaCS<>::wait

#endif // P2PPolicyDaCS_hxx
