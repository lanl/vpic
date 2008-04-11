/*
	Definition of CMPolicyMPIDaCS class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef CMPolicyMPIDaCS_hxx
#define CMPolicyMPIDaCS_hxx

#include <iostream>
#include <assert.h>
#include <mpi.h>
#include <dacs.h>
#include <DaCSUtils.h>
#include <MPData.hxx>

/*!
	\struct CMPolicyMPIDaCS CMPolicyMPIDaCS.h
	\brief  provides...
*/
class CMPolicyMPIDaCS
	{
	public:

		void init(int argc, char ** argv);
		void finalize();

		inline int global_id() { return rank_; }
		inline int global_size() { return size_; }

		inline de_id_t peer_de() { return peer_de_; }
		inline dacs_process_id_t peer_pid() { return peer_pid_; }
		inline dacs_group_t group() { return group_; }

		inline MPI_Comm dmp_comm() { return MPI_COMM_WORLD; }

	protected:

		CMPolicyMPIDaCS()
			: rank_(-1), size_(-1) {}
		CMPolicyMPIDaCS(const CMPolicyMPIDaCS & mgr) {}
		~CMPolicyMPIDaCS() {}

	private:

		int rank_;
		int size_;

		DACS_ERR_T errcode_;
		de_id_t peer_de_;
		dacs_process_id_t peer_pid_;
		dacs_group_t group_;

	}; // class CMPolicyMPIDaCS

void CMPolicyMPIDaCS::init(int argc, char ** argv)
	{
		// check command-line arguments
		if(argc < 2) {
			std::cerr << "Usage: " << argv[0] <<
				" PPE executable (with full path)" << std::endl;
			exit(1);
		} // if

		// MPI initialization
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
		MPI_Comm_size(MPI_COMM_WORLD, &size_);

		// DaCS initialization
		//errcode_ = dacs_runtime_init(NULL, NULL);	
		errcode_ = dacs_init(DACS_INIT_SINGLE_THREADED);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		std::cerr << "host process initialized" << std::endl;

		// create a group for this process and its child
		errcode_ = dacs_group_init(&group_, 0);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// add this process to group
		errcode_ = dacs_group_add_member(DACS_DE_SELF,
			DACS_PID_SELF, group_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);
		
		// check for available children and try to reserve one
		uint32_t children;
		errcode_ = dacs_get_num_avail_children(DACS_DE_CBE, &children);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);
		assert(children>0);
		std::cerr << "children " << children << std::endl;
		children = 1;

		errcode_ = dacs_reserve_children(DACS_DE_CBE, &children, &peer_de_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);


        // get args to pass on to child
        const char ** args = argc > 2 ?
            const_cast<const char **>(&argv[2]) : NULL;

#if 1
		errcode_ = dacs_de_start(peer_de_, argv[1], args, NULL,
			DACS_PROC_LOCAL_FILE, &peer_pid_);	
#else
		errcode_ = dacs_de_start(peer_de_, argv[1], args, NULL,
			DACS_PROC_REMOTE_FILE, &peer_pid_);	
#endif
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// add the new child process to the group
		errcode_ = dacs_group_add_member(peer_de_, peer_pid_, group_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// need to call this to finalize the group
		// after this, no new members can be added
		// this call block until all children have
		// accepted their invitations
		errcode_ = dacs_group_close(group_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		dacs_wid_t wid;
		errcode_ = dacs_wid_reserve(&wid);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		int data[2] __attribute__ ((aligned (16)));

		/*
		void * mem(NULL);
		posix_memalign((void **)(&mem), 16, 2*sizeof(int));
		int * data __attribute__ ((aligned (16))) = (int *)mem;
		*/

		data[0] = rank_;
		data[1] = size_;
		errcode_ = dacs_send(data, 2*sizeof(int),
			peer_de_, peer_pid_, 0, wid, DACS_BYTE_SWAP_WORD);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wait(wid);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		errcode_ = dacs_wid_release(&wid);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		//free(mem);
	} // CMPolicyMPIDaCS::init

void CMPolicyMPIDaCS::finalize()
	{
		// destroy the group, this will implicitly remove
		// this process from the group
		errcode_ = dacs_group_destroy(&group_);
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// DaCS finalization
		//errcode_ = dacs_runtime_exit();
		errcode_ = dacs_exit();
		process_dacs_errcode(errcode_, __FILE__, __LINE__);

		// MPI finalization
		MPI_Finalize();
	} // CMPolicyMPIDaCS::finalize

#endif // CMPolicyMPIDaCS_hxx
