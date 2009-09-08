/*
	Definition of CMPolicyMultipleContextMPI class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef CMPolicyMultipleContextMPI_hxx
#define CMPolicyMultipleContextMPI_hxx

#include <mpi.h>
#include <assert.h>
#include <iostream>
#include "MPData.hxx"

/*!
	\struct CMPolicyMultipleContextMPI_T CMPolicyMultipleContextMPI_T.h
	\brief  provides...
*/
template<int ROLE> class CMPolicyMultipleContextMPI_T
	{
	public:

		static CMPolicyMultipleContextMPI_T & instance() {
			static CMPolicyMultipleContextMPI_T mpi;
			return mpi;
		} // instance

		void init(int * pargc, char *** pargv);
		void finalize();

		// return local id and size
		inline int local_id() { return p2p_rank_; }
		inline int local_size() { return p2p_size_; }

		// return global id and size
		inline int global_id() { return dmp_rank_; }
		inline int global_size() { return dmp_size_; }

		inline int peer_p2p_rank() { return peer_p2p_rank_; }
		inline int peer_world_rank() { return peer_world_rank_; }

		inline MPI_Comm p2p_comm() { return p2p_comm_; }
		inline MPI_Comm dmp_comm() { return dmp_comm_; }

	protected:

		CMPolicyMultipleContextMPI_T()
			: world_rank_(-1), world_size_(-1), peer_world_rank_(-1),
			p2p_rank_(-1), p2p_size_(-1), peer_p2p_rank_(-1),
			dmp_rank_(-1), dmp_size_(-1), p2p_comm_(NULL) {}
		CMPolicyMultipleContextMPI_T(const CMPolicyMultipleContextMPI_T & mpi)
			{}
		~CMPolicyMultipleContextMPI_T() {}

	private:

		int world_rank_;
		int world_size_;
		int peer_world_rank_;
		int p2p_rank_;
		int p2p_size_;
		int peer_p2p_rank_;
		int dmp_rank_;
		int dmp_size_;

		MPI_Comm p2p_comm_;
		MPI_Comm dmp_comm_;

		// Borrowed from Paul Henning
		inline void initialize_communicators();

		inline void propagate_dmp_rank();

	}; // class CMPolicyMultipleContextMPI_T

template<int ROLE>
void CMPolicyMultipleContextMPI_T<ROLE>::init(int * pargc, char *** pargv)
	{
		// initialize mpi
		MPI_Init(pargc, pargv);

		// get world rank
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);

		// get world size
		MPI_Comm_size(MPI_COMM_WORLD, &world_size_);

		// setup communicators
		initialize_communicators();
	} // CMPolicyMultipleContextMPI_T<>::initialize

template<int ROLE>
void CMPolicyMultipleContextMPI_T<ROLE>::finalize()
	{
		MPI_Finalize();
	} // CMPolicyMultipleContextMPI_T<>::finalize

// Borrowed from Paul Henning
template<int ROLE> inline
void CMPolicyMultipleContextMPI_T<ROLE>::initialize_communicators()
	{
		// no static assertions yet
		assert(ROLE == 0 || ROLE == 1);

		int * role_type[2] = {NULL, NULL};	
		int role_size[2] = {0, 0};

		int * rank_role = new int[world_size_];

		// gather roles for all ranks
		int role = ROLE;
		MPI_Allgather(&role, 1, MPI_INT, rank_role, 1,
			MPI_INT, MPI_COMM_WORLD);

		for(int i(0); i<world_size_; i++) {
			assert(rank_role[i] > -1 && rank_role[i] < 2);
			role_size[rank_role[i]] += 1;
		} // for

		// check to make sure that there are enough
		// of each role type
		if(role_size[0] != role_size[1]) {
			std::cerr << "Role sizes must agree!" << std::endl;
			exit(1);
		} // if

		// allocate some space for sorting roles
		role_type[0] = new int[role_size[0]];
		role_type[1] = new int[role_size[1]];

		int current[2] = {0, 0};

		// index roles by rank
		for(int i(0); i<world_size_; i++) {
			register const int role = rank_role[i];
			role_type[role][current[role]++] = i;
		} // for

		// find my offset in my role type
		int offset(0);
		for(; offset<role_size[ROLE]; offset++) {
			if(role_type[ROLE][offset] == world_rank_) {
				break;
			} // if
		} // for

		assert(offset<role_size[ROLE]);

		// get my peer rank
		peer_world_rank_ = role_type[ROLE^1][offset];

		// get world group
		MPI_Group world_group;
		MPI_Comm_group(MPI_COMM_WORLD, &world_group);

		// create dmp communicator
		MPI_Group dmp_group;
		MPI_Group_incl(world_group, role_size[MP_HOST],
			role_type[MP_HOST], &dmp_group);
		MPI_Comm_create(MPI_COMM_WORLD, dmp_group, &dmp_comm_);

		// create a private communicator
		int p2p_ranks[2];
		p2p_ranks[ROLE%2] = world_rank_;
		p2p_ranks[(ROLE+1)%2] = peer_world_rank_;

		MPI_Group p2p_group;
		MPI_Group_incl(world_group, 2, p2p_ranks, &p2p_group);
		MPI_Comm_create(MPI_COMM_WORLD, p2p_group, &p2p_comm_);

		// get p2p rank and size
		MPI_Comm_rank(p2p_comm_, &p2p_rank_);
		peer_p2p_rank_ = p2p_rank_^1;
		MPI_Comm_size(p2p_comm_, &p2p_size_);
		
		// set dmp rank and size
		propagate_dmp_rank();
		
		// cleanup
		MPI_Group_free(&world_group);
		MPI_Group_free(&dmp_group);
		MPI_Group_free(&p2p_group);

		delete rank_role;
		delete role_type[0];
		delete role_type[1];
	} // CMPolicyMultipleContextMPI_T<>::initialize_communicators

template<> inline
void CMPolicyMultipleContextMPI_T<MP_HOST>::propagate_dmp_rank()
	{
		MPI_Comm_rank(dmp_comm_, &dmp_rank_);
		MPI_Comm_size(dmp_comm_, &dmp_size_);

		// send to accelerator
		int info[2];
		info[0] = dmp_rank_;
		info[1] = dmp_size_;
		MPI_Send(info, 2, MPI_INT, peer_p2p_rank_, 0, p2p_comm_);
	} // CMPolicyMultipleContextMPI_T<>::propagate_dmp_rank

template<> inline
void CMPolicyMultipleContextMPI_T<MP_ACCEL>::propagate_dmp_rank()
	{
		// receive from host
		int info[2];
		MPI_Status status;
		MPI_Recv(info, 2, MPI_INT, peer_p2p_rank_, 0,
			p2p_comm_, &status);
		dmp_rank_ = info[0];
		dmp_size_ = info[1];
	} // CMPolicyMultipleContextMPI_T<>::propagate_dmp_rank

#if defined HOST_BUILD
	typedef CMPolicyMultipleContextMPI_T<MP_HOST> CMPolicyMultipleContextMPI;
#else
	typedef CMPolicyMultipleContextMPI_T<MP_ACCEL> CMPolicyMultipleContextMPI;
#endif // build type

#endif // CMPolicyMultipleContextMPI_hxx
