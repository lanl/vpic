/*
	Definition of MPICommunicationPolicy class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef MPICommunicationPolicy_hxx
#define MPICommunicationPolicy_hxx

#include <mpi.h>
#include <iostream>
#include <string>
#include <assert.h>
#include <P2PData.hxx>
#include <Type2MPIType.hxx>

#define HOST 0
#define ACCEL 1

/*!
	\struct MPICommunicationPolicy MPICommunicationPolicy.h
	\brief  provides...
*/
template<int ROLE> class MPICommunicationPolicy
	{
	public:

		void init(int argc, char ** argv);
		void finalize() { MPI_Finalize(); }

		// return dmp rank and size
		inline int rank() { return dmp_rank_; }
		inline int size() { return dmp_size_; }

		int peer_rank() { return p2p_peer_rank_; }

		int poll(P2PHeader & header);
		int request(int p2ptag, int tag = P2PTag::data,
			int count = 0, int rank = 0);

		template<typename T> int send(T * buffer, int count,
			int tag = P2PTag::data);
		template<typename T> int recv(T * buffer, int count,
			int tag = P2PTag::data);

		template<typename T> int isend(T * buffer, int count,
			int tag, MPI_Request & request);
		template<typename T> int irecv(T * buffer, int count,
			int tag, MPI_Request & request);

		int wait(MPI_Request & request, MPI_Status & status);
		template<typename T> int get_count(MPI_Status & status,
			int & count, T * dummy = NULL);
		int sync();

		double wtime() { return MPI_Wtime(); }

		MPI_Comm & peer_comm() { return p2p_comm_; }
		MPI_Comm & dmp_comm() { return dmp_comm_; }

	protected:

		MPICommunicationPolicy()
			: world_rank_(-1), world_size_(-1), peer_world_rank_(-1),
			p2p_rank_(-1), p2p_size_(-1), p2p_peer_rank_(-1),
			dmp_rank_(-1), dmp_size_(-1), p2p_comm_(NULL), dmp_comm_(NULL)
			{
			} // MPICommunicationPolicy::MPICommunicationPolicy

		~MPICommunicationPolicy()
			{}

		int world_rank_;
		int world_size_;
		int peer_world_rank_;
		int p2p_rank_;
		int p2p_size_;
		int p2p_peer_rank_;
		int dmp_rank_;
		int dmp_size_;

		MPI_Comm p2p_comm_;
		MPI_Comm dmp_comm_;

	private:

		// Borrowed from Paul Henning
		void initialize_communicators();

		void propagate_dmp_rank();

	}; // class MPICommunicationPolicy

template<int ROLE>
void MPICommunicationPolicy<ROLE>::init(int argc, char ** argv)
	{
		// initialize mpi
		MPI_Init(&argc, &argv);

		// get world rank
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);

		// get world size
		MPI_Comm_size(MPI_COMM_WORLD, &world_size_);

		// setup communicators
		initialize_communicators();
		
		/*
		// output some info
		std::string role_str = ROLE == HOST ? "host" : "accelerator";
		std::cout << role_str << ": " << 
			" world rank: " << world_rank_ << 
			" world size: " << world_size_ << 
			" peer world rank: " << peer_world_rank_ << 
			" p2p rank: " << p2p_rank_ << 
			" p2p size: " << p2p_size_ << 
			" peer p2p rank: " << p2p_peer_rank_ <<
			" dmp rank: " << dmp_rank_ <<
			" dmp size: " << dmp_size_ << std::endl;
		*/
	} // MPICommunicationPolicy<>::initialize

template<>
void MPICommunicationPolicy<HOST>::finalize()
	{
		MPI_Finalize();
	} // MPICommunicationPolicy<>::finalize

template<>
void MPICommunicationPolicy<ACCEL>::finalize()
	{
		P2PHeader header;
		MPI_Send(&header, 3, MPI_INT, p2p_peer_rank_, P2PTag::end, p2p_comm_);
		MPI_Finalize();
	} // MPICommunicationPolicy<>::finalize

template<>
int MPICommunicationPolicy<HOST>::poll(P2PHeader & header)
	{
		MPI_Status status;

		MPI_Probe(p2p_peer_rank_, MPI_ANY_TAG, p2p_comm_, &status);

		const int tag = status.MPI_TAG;

		MPI_Recv(&header, 3, MPI_INT, p2p_peer_rank_, tag,
			p2p_comm_, &status);

		return tag;
	} // MPICommunicationPolicy<>::poll

template<>
int MPICommunicationPolicy<ACCEL>::request(int p2ptag, int tag,
	int count, int rank)
	{
		P2PHeader header(tag, count, rank);
		return MPI_Send(&header, 3, MPI_INT, p2p_peer_rank_,
			p2ptag, p2p_comm_);
	} // MPICommunicatorPolicy<>::request

template<int ROLE>
template<typename T>
int MPICommunicationPolicy<ROLE>::send(T * buffer, int count, int tag)
	{
		return MPI_Send(buffer, count, Type2MPIType<T>::type(),
			p2p_peer_rank_, tag, p2p_comm_);
	} // MPICommunicatorPolicy<>::send

template<int ROLE>
template<typename T>
int MPICommunicationPolicy<ROLE>::recv(T * buffer, int count, int tag)
	{
		MPI_Status status;
		return MPI_Recv(buffer, count, Type2MPIType<T>::type(),
			p2p_peer_rank_, tag, p2p_comm_, &status);
	} // MPICommunicatorPolicy<>::send

template<>
template<typename T>
int MPICommunicationPolicy<ACCEL>::isend(T * buffer,
	int count, int tag, MPI_Request & request)
	{
		/*
		return MPI_Send(buffer, count, Type2MPIType<T>::type(),
			p2p_peer_rank_, tag, p2p_comm_);
		*/
		///*
		return MPI_Isend(buffer, count, Type2MPIType<T>::type(),
			p2p_peer_rank_, tag, p2p_comm_, &request);
		//*/
	} // MPICommunicatorPolicy<>::send

template<>
template<typename T>
int MPICommunicationPolicy<ACCEL>::irecv(T * buffer,
	int count, int tag, MPI_Request & request)
	{
		/*
		MPI_Status status;
		return MPI_Recv(buffer, count, Type2MPIType<T>::type(),
			p2p_peer_rank_, tag, p2p_comm_, &status);
		*/
		///*
		return MPI_Irecv(buffer, count, Type2MPIType<T>::type(),
			p2p_peer_rank_, tag, p2p_comm_, &request);
		//*/
	} // MPICommunicatorPolicy<>::send

template<int ROLE>
int MPICommunicationPolicy<ROLE>::wait(MPI_Request & request,
	MPI_Status & status)
	{
		return MPI_Wait(&request, &status);
		//return MPI_SUCCESS;
	} // MPICommunicationPolicy<>::wait

template<>
template<typename T>
int MPICommunicationPolicy<ACCEL>::get_count(MPI_Status & status,
	int & count, T * dummy)
	{
		return MPI_Get_count(&status, Type2MPIType<T>::type(), &count);
	} // MPICommunicationPolicy<>::wait

template<>
int MPICommunicationPolicy<HOST>::sync()
	{
		P2PHeader header;
		return MPI_Send(&header, 2, MPI_INT, p2p_peer_rank_,
			P2PTag::data, p2p_comm_);
	} // MPICommunicationPolicy<>::sync

template<>
int MPICommunicationPolicy<ACCEL>::sync()
	{
		MPI_Status status;
		P2PHeader header;
		return MPI_Recv(&header, 2, MPI_INT, p2p_peer_rank_,
			P2PTag::data, p2p_comm_, &status);
	} // MPICommunicationPolicy<>::sync

// Borrowed from Paul Henning
template<int ROLE>
void MPICommunicationPolicy<ROLE>::initialize_communicators()
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
		MPI_Group_incl(world_group, role_size[HOST],
			role_type[HOST], &dmp_group);
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
		p2p_peer_rank_ = p2p_rank_^1;
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
	} // MPICommunicationPolicy<>::initialize_communicators

template<>
void MPICommunicationPolicy<HOST>::propagate_dmp_rank()
	{
		MPI_Comm_rank(dmp_comm_, &dmp_rank_);
		MPI_Comm_size(dmp_comm_, &dmp_size_);

		// send to accelerator
		int info[2];
		info[0] = dmp_rank_;
		info[1] = dmp_size_;
		MPI_Send(info, 2, MPI_INT, p2p_peer_rank_, 0, p2p_comm_);
	} // MPICommunicationPolicy<>::propagate_dmp_rank

template<>
void MPICommunicationPolicy<ACCEL>::propagate_dmp_rank()
	{
		// receive from host
		int info[2];
		MPI_Status status;
		MPI_Recv(info, 2, MPI_INT, p2p_peer_rank_, 0,
			p2p_comm_, &status);
		dmp_rank_ = info[0];
		dmp_size_ = info[1];
	} // MPICommunicationPolicy<>::propagate_dmp_rank

#endif // MPICommunicationPolicy_hxx
