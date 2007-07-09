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

#include <dacs.h>
#include <P2PTag.hxx>
#include <MPData.hxx>

/*!
	\struct P2PPolicyDaCS P2PPolicyDaCS.h
	\brief  provides...
*/
template<int ROLE> class P2PPolicyDaCS
	{
	public:

		void init(int argc, char ** argv);
		void finalize();

		// topology information
		inline int global_id()
			{ return MPIConnectionManager::instance().global_id(); }
		inline int global_size()
			{ return MPIConnectionManager::instance().global_size(); }

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

		P2PPolicyDaCS() {}
		~P2PPolicyDaCS() {}

	private:

	}; // class P2PPolicyDaCS

template<int ROLE> inline
void P2PPolicyDaCS<ROLE>::init(int argc, char ** argv)
	{
	} // P2PPolicyDaCS<>::init

template<int ROLE> inline
void P2PPolicyDaCS<ROLE>::finalize()
	{
	} // P2PPolicyDaCS<>::init

template<> inline
int P2PPolicyDaCS<MP_HOST>::poll(MPRequest_T<MP_HOST> & request)
	{
	} // P2PPolicyDaCS<>::poll

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::post(int p2ptag)
	{
	} // MPICommunicatorPolicy<>::request

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::post(MPRequest_T<MP_ACCEL> & request)
	{
	} // MPICommunicatorPolicy<>::request

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::send(T * buffer, int count, int tag)
	{
	} // MPICommunicatorPolicy<>::isend

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::recv(T * buffer, int count, int tag, int id)
	{
	} // MPICommunicatorPolicy<>::isend

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::isend(T * buffer, int count, int tag, int id)
	{
	} // MPICommunicatorPolicy<>::isend

template<int ROLE>
template<typename T>
int P2PPolicyDaCS<ROLE>::irecv(T * buffer, int count, int tag, int id)
	{
	} // MPICommunicatorPolicy<>::isend

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::wait_send(int id)
	{
	} // P2PPolicyDaCS<>::wait

template<> inline
int P2PPolicyDaCS<MP_ACCEL>::wait_recv(int id)
	{
	} // P2PPolicyDaCS<>::wait

template<int ROLE> inline
int P2PPolicyDaCS<ROLE>::barrier()
	{
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
