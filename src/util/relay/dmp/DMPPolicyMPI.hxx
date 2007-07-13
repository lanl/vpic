/*
	Definition of DMPPolicyMPI_T class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$
	vim: set ts=3 :
*/

#ifndef DMPPolicyMPI_hxx
#define DMPPolicyMPI_hxx

#include <mpi.h>
#include <MPData.hxx>
#include <ConnectionManager.hxx>
#include <Type2MPIType.hxx>

/*!
	\class DMPPolicyMPI DMPPolicyMPI.h
	\brief  provides...
*/
class DMPPolicyMPI
	{
	public:

		// topology information
		inline int global_id()
			{ return ConnectionManager::instance().global_id(); }
		inline int global_size()
			{ return ConnectionManager::instance().global_size(); }

		template<typename T> int send(T * buffer, int count,
			int dst, int tag);
		template<typename T> int recv(T * buffer, int count,
			int src, int tag, int id);
		template<typename T> int isend(T * buffer, int count,
			int dst, int tag, int id);
		template<typename T> int irecv(T * buffer, int count,
			int src, int tag, int id);

		int test_recv(MPRequest_T<MP_HOST> & request);
		int test_send(MPRequest_T<MP_HOST> & request);
		int wait_recv(MPRequest_T<MP_HOST> & request);
		int wait_send(MPRequest_T<MP_HOST> & request);
		int barrier();

		template<typename T> int allreduce_max(T * src, T * dest,
			int count);
		template<typename T> int allreduce_sum(T * src, T * dest,
			int count);
		template<typename T> int allgather(T * src, T * dest, int count);

		int abort(int reason);

	protected:

		DMPPolicyMPI() {}
		~DMPPolicyMPI() {}

	private:

		MPI_Request send_request_[max_buffers];
		MPI_Request recv_request_[max_buffers];
		MPI_Status recv_status_[max_buffers];

	}; // class DMPPolicyMPI

template<typename T>
int DMPPolicyMPI::send(T * buffer, int count, int dst, int tag)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Send(buffer, count, Type2MPIType<T>::type(),
			dst, tag, mgr.dmp_comm());
	} // DMPPolicyMPI::send

template<typename T>
int DMPPolicyMPI::recv(T * buffer, int count, int src,
	int tag, int id)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Recv(buffer, count, Type2MPIType<T>::type(),
			src, tag, mgr.dmp_comm(), &recv_status_[id]);
	} // DMPPolicyMPI::irecv

template<typename T>
int DMPPolicyMPI::isend(T * buffer, int count, int dst,
	int tag, int id)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Isend(buffer, count, Type2MPIType<T>::type(),
			dst, tag, mgr.dmp_comm(), &send_request_[id]);
	} // DMPPolicyMPI::send

template<typename T>
int DMPPolicyMPI::irecv(T * buffer, int count, int src,
	int tag, int id)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Irecv(buffer, count, Type2MPIType<T>::type(),
			src, tag, mgr.dmp_comm(), &recv_request_[id]);
	} // DMPPolicyMPI::irecv

int DMPPolicyMPI::test_recv(MPRequest_T<MP_HOST> & request)
	{
		int flag;
		MPI_Status status; // FIXME
		int err = MPI_Test(&recv_request_[request.id], &flag, &status);

		if(flag) {
			request.state = complete;
		} // if

		return err;
	} // DMPPolicyMPI::test

int DMPPolicyMPI::test_send(MPRequest_T<MP_HOST> & request)
	{
		int flag;
		MPI_Status status; // FIXME
		int err = MPI_Test(&send_request_[request.id], &flag, &status);

		if(flag) {
			request.state = complete;
		} // if

		return err;
	} // DMPPolicyMPI::test

int DMPPolicyMPI::wait_recv(MPRequest_T<MP_HOST> & request)
	{
		MPI_Status status; // FIXME
		return MPI_Wait(&recv_request_[request.id], &status);
	} // DMPPolicyMPI::test

int DMPPolicyMPI::wait_send(MPRequest_T<MP_HOST> & request)
	{
		MPI_Status status; // FIXME
		return MPI_Wait(&send_request_[request.id], &status);
	} // DMPPolicyMPI::test

int DMPPolicyMPI::barrier()
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Barrier(mgr.dmp_comm());
	} // DMPPolicyMPI::test

template<typename T>
int  DMPPolicyMPI::allreduce_max(T * src, T * dest, int count)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Allreduce(src, dest, count, Type2MPIType<T>::type(),
			MPI_MAX, mgr.dmp_comm());
	} // DMPPolicyMPI::allreduce_max

template<typename T>
int  DMPPolicyMPI::allreduce_sum(T * src, T * dest, int count)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Allreduce(src, dest, count, Type2MPIType<T>::type(),
			MPI_SUM, mgr.dmp_comm());
	} // DMPPolicyMPI::allreduce_max

template<typename T>
int  DMPPolicyMPI::allgather(T * src, T * dest, int count)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Allgather(src, count, Type2MPIType<T>::type(),
			dest, count, Type2MPIType<T>::type(), mgr.dmp_comm());
	} // DMPPolicyMPI::allreduce_max

int DMPPolicyMPI::abort(int reason)
	{
		ConnectionManager & mgr = ConnectionManager::instance();
		return MPI_Abort(mgr.dmp_comm(), reason);
	} // DMPPolicyMPI::test

		int abort(int reason);

#endif // DMPPolicyMPI_hxx
