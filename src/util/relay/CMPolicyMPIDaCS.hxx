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
		de_id_t de_;
		de_id_t peer_de_;
		dacs_process_id_t peer_pid_;

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
		errcode_ = dacs_runtime_init(NULL, NULL, &de_);	
		process_dacs_errcode(errcode_);

		std::cerr << "host process initialized" << std::endl;

///* THIS WILL CHANGE WITH THE NEXT DACS VERSION
		errcode_ = dacs_topology_get_child(de_, 0, &peer_de_);
		process_dacs_errcode(errcode_);

		DACS_DE_TYPE_T type;
		errcode_ = dacs_topology_get_type(peer_de_, &type);
		process_dacs_errcode(errcode_);

		if(type == DACS_DE_CELL_BLADE) {
			std::cerr << "encountered blade" << std::endl;
			de_id_t blade_de = peer_de_;
			errcode_ = dacs_topology_get_child(blade_de, 0, &peer_de_);
			process_dacs_errcode(errcode_);

			errcode_ = dacs_topology_get_type(peer_de_, &type);
			process_dacs_errcode(errcode_);
			assert(type == DACS_DE_CBE);
		} // if
//*/ THIS WILL CHANGE WITH THE NEXT DACS VERSION

		errcode_ = dacs_topology_reserve(peer_de_);
		process_dacs_errcode(errcode_);

		errcode_ = dacs_de_start(peer_de_, argv[1], NULL, NULL,
			DACS_PROC_SYSTEMX_FILE, &peer_pid_);	
		std::cerr << "host called dacs_de_start " << argv[1] << " " << errcode_ << std::endl;

		dacs_tag_t tag;
		errcode_ = dacs_tag_reserve(&tag, peer_de_, peer_pid_);
		std::cerr << "host called dacs_tag_reserve " << errcode_ << std::endl;
		process_dacs_errcode(errcode_);

		int data[2] __attribute__ ((aligned (16)));

		/*
		void * mem(NULL);
		posix_memalign((void **)(&mem), 16, 2*sizeof(int));
		int * data __attribute__ ((aligned (16))) = (int *)mem;
		*/

		data[0] = rank_;
		data[1] = size_;
		errcode_ = dacs_send(peer_de_, peer_pid_, data,
			2*sizeof(int), tag, DACS_BYTE_SWAP_WORD);

		std::cerr << "host called dacs_send " << errcode_ << std::endl;
		sleep(2);
		process_dacs_errcode(errcode_);

		errcode_ = dacs_wait(tag);
		std::cerr << "host called dacs_wait " << errcode_ << std::endl;
		sleep(2);
		process_dacs_errcode(errcode_);

		errcode_ = dacs_tag_release(&tag);
		std::cerr << "host called dacs_tag_release" << std::endl;
		process_dacs_errcode(errcode_);

		//free(mem);
	} // CMPolicyMPIDaCS::init

void CMPolicyMPIDaCS::finalize()
	{
		// DaCS finalization
		errcode_ = dacs_runtime_exit();
		process_dacs_errcode(errcode_);

		// MPI finalization
		MPI_Finalize();
	} // CMPolicyMPIDaCS::finalize

#endif // CMPolicyMPIDaCS_hxx
