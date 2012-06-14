/*
	Definition of CMPolicyDaCS class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef CMPolicyDaCS_hxx
#define CMPolicyDaCS_hxx

#include <iostream>
#include <assert.h>
#include <dacs.h>
#include <DaCSUtils.h>
#include <MPData.hxx>

/*!
	\struct CMPolicyDaCS CMPolicyDaCS.h
	\brief  provides...
*/
class CMPolicyDaCS
	{
	public:

		inline void init(int * pargc, char *** pargv);
		inline void finalize();

		inline int global_id() { return id_; }
		inline int global_size() { return size_; }

		inline de_id_t peer_de() { return peer_de_; }
		inline dacs_process_id_t peer_pid() { return peer_pid_; }
		inline dacs_group_t group() { return group_; }

	protected:

		CMPolicyDaCS()
			: id_(-1), size_(-1), peer_de_(DACS_DE_PARENT),
			peer_pid_(DACS_PID_PARENT) {}
		CMPolicyDaCS(const CMPolicyDaCS & mgr) {}
		~CMPolicyDaCS() {}

	private:

		int id_;
		int size_;

		DACS_ERR_T errcode_;
		de_id_t peer_de_;
		dacs_process_id_t peer_pid_;
		dacs_group_t group_;

	}; // class CMPolicyDaCS

void CMPolicyDaCS::init(int * pargc, char *** pargv)
	{
		// initialize runtime environment
		//errcode_ = dacs_runtime_init(NULL, NULL);	
		errcode_ = dacs_init(DACS_INIT_SINGLE_THREADED);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// this process must accept the HE offer of group membership
		// execution will block until the HE has made the offer and
		// all members have accepted
		errcode_ = dacs_group_accept(peer_de(), peer_pid(), &group_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// this gets the MPI rank and size information from
		// the HE process
		dacs_wid_t wid;
		errcode_ = dacs_wid_reserve(&wid);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		int data[2] __attribute__ ((aligned (16)));

		/*
		void * mem(NULL);
		posix_memalign((void **)(&mem), 16, 2*sizeof(int));
		int * data __attribute__ ((aligned (16))) = (int *)mem;
		*/

		errcode_ = dacs_recv(data, 2*sizeof(int), peer_de_,
			peer_pid_, 0, wid, DACS_BYTE_SWAP_WORD);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wait(wid);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wid_release(&wid);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		id_ = data[0];
		size_ = data[1];

		//free(mem);
	} // CMPolicyDaCS::init

void CMPolicyDaCS::finalize()
	{
		// this process must request permission to leave the
		// group before the group can be destroyed on the
		// HE process
		errcode_ = dacs_group_leave(&group_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		//errcode_ = dacs_runtime_exit();
		errcode_ = dacs_exit();
		process_dacs_errcode(errcode_, __FILE__, __LINE__);
	} // CMPolicyDaCS::finalize

#endif // CMPolicyDaCS_hxx
