/*
	Definition of P2PPolicyMPI class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef P2PPolicyMPI_hxx
#define P2PPolicyMPI_hxx

#include <mpi.h>
#include <ConnectionManager.hxx>
#include <P2PTag.hxx>
#include <MPData.hxx>
#include <Type2MPIType.hxx>

/*!
	\struct P2PPolicyMPI P2PPolicyMPI.h
	\brief  provides...
*/
template<int ROLE> class P2PPolicyMPI
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

		P2PPolicyMPI() {}
		~P2PPolicyMPI() {}

		MPI_Request send_request_[max_buffers];
		MPI_Request recv_request_[max_buffers];
		MPI_Status status_[max_buffers];

	private:

	}; // class P2PPolicyMPI

template<> inline
int P2PPolicyMPI<MP_HOST>::poll(MPRequest_T<MP_HOST> & request)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		MPI_Status status;
		int flag;

		MPI_Iprobe(mgr.peer_p2p_rank(), P2PTag::request, mgr.p2p_comm(),
			&flag, &status);

		if(flag) {
			MPI_Recv(&request, request_count(), MPI_INT, mgr.peer_p2p_rank(),
				P2PTag::request, mgr.p2p_comm(), &status);

			return request.p2ptag;
		}
		else {
			return P2PTag::pending;
		} // if
	} // P2PPolicyMPI<>::poll

template<> inline
int P2PPolicyMPI<MP_ACCEL>::post(int p2ptag)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		MPRequest_T<MP_ACCEL> request(p2ptag);
		return MPI_Send(&request, request_count(), MPI_INT,
			mgr.peer_p2p_rank(), P2PTag::request, mgr.p2p_comm());
	} // MPICommunicatorPolicy<>::request

template<> inline
int P2PPolicyMPI<MP_ACCEL>::post(MPRequest_T<MP_ACCEL> & request)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Send(&request, request_count(), MPI_INT,
			mgr.peer_p2p_rank(), P2PTag::request, mgr.p2p_comm());
	} // MPICommunicatorPolicy<>::request

template<int ROLE>
template<typename T>
int P2PPolicyMPI<ROLE>::send(T * buffer, int count, int tag)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Send(buffer, count, Type2MPIType<T>::type(),
			mgr.peer_p2p_rank(), tag, mgr.p2p_comm());
	} // MPICommunicatorPolicy<>::isend

template<int ROLE>
template<typename T>
int P2PPolicyMPI<ROLE>::recv(T * buffer, int count, int tag, int id)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Recv(buffer, count, Type2MPIType<T>::type(),
			mgr.peer_p2p_rank(), tag, mgr.p2p_comm(), &status_[id]);
	} // MPICommunicatorPolicy<>::isend

template<int ROLE>
template<typename T>
int P2PPolicyMPI<ROLE>::isend(T * buffer, int count, int tag, int id)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Isend(buffer, count, Type2MPIType<T>::type(),
			mgr.peer_p2p_rank(), tag, mgr.p2p_comm(), &send_request_[id]);
	} // MPICommunicatorPolicy<>::isend

template<int ROLE>
template<typename T>
int P2PPolicyMPI<ROLE>::irecv(T * buffer, int count, int tag, int id)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Irecv(buffer, count, Type2MPIType<T>::type(),
			mgr.peer_p2p_rank(), tag, mgr.p2p_comm(), &recv_request_[id]);
	} // MPICommunicatorPolicy<>::irecv

template<> inline
int P2PPolicyMPI<MP_HOST>::test_send(MPRequest_T<MP_HOST> & request)
	{
		int flag;
		MPI_Status status; // FIXME
		int err = MPI_Test(&send_request_[request.id], &flag, &status);

		if(flag) {
			request.state = complete;
		} // if

		return err;
	} // P2PPolicyMPI<>::test_send

template<int ROLE> inline
int P2PPolicyMPI<ROLE>::wait_send(int id)
	{
		return MPI_Wait(&send_request_[id], &status_[id]);
	} // P2PPolicyMPI<>::wait_send

template<> inline
int P2PPolicyMPI<MP_ACCEL>::wait_recv(int id)
	{
		return MPI_Wait(&recv_request_[id], &status_[id]);
	} // P2PPolicyMPI<>::wait_recv

template<int ROLE> inline
int P2PPolicyMPI<ROLE>::barrier()
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Barrier(mgr.p2p_comm());
	} // P2PPolicyMPI<>::barrier

template<>
template<typename T>
int P2PPolicyMPI<MP_ACCEL>::get_count(int id, int & count, T * dummy)
	{
		return MPI_Get_count(&status_[id], Type2MPIType<T>::type(),
			&count);
	} // P2PPolicyMPI<>::get_count

template<int ROLE>
int P2PPolicyMPI<ROLE>::abort(int reason)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Abort(mgr.p2p_comm(), reason);
	} // P2PPolicyMPI<>::abort

#endif // P2PPolicyMPI_hxx
