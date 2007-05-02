/*------------------------------------------------------------------------------
------------------------------------------------------------------------------*/

#include <iostream>
#include <P2PConnection.hxx>

struct Buffer {
	int rank;
	int size;
}; // struct buffer

int main(int argc, char *argv[]) {

	// get an instance of the manager
	P2PConnection & p2p = P2PConnection::instance();		

	// initialize everything
	p2p.init(argc, argv);

	int rank, size;
	Buffer * msg;

	rank = p2p.rank();
	size = p2p.size();

	MPI_Request * sreq = new MPI_Request[size];
	MPI_Request * rreq = new MPI_Request[size];
	msg = new Buffer[size];

	// send to other ranks
	for(int i(0); i<size; i++) {
		if(i != rank) {
			// send a message to every other accelerator
			// through the host relay
			msg[i].rank = rank;
			msg[i].size = size;
			p2p.request(P2PTag::isend, P2PTag::data, 8, i);
			p2p.isend(reinterpret_cast<char *>(&msg[i]), 8,
				P2PTag::data, sreq[i]);
		} // if
	} // for

	// receive from other ranks
	for(int i(0); i<size; i++) {
		if(i != rank) {
			// receive a message from every other accelerator
			// through the host relay
			p2p.request(P2PTag::irecv, P2PTag::data, 8, i);
			p2p.irecv(reinterpret_cast<char *>(&msg[i]), 8,
				P2PTag::data, rreq[i]);
		} // if
	} // for

	// wait on messages
	MPI_Status status;
	for(int i(0); i<size; i++) {
		if(i != rank) {
			p2p.wait(sreq[i], status);
			p2p.wait(rreq[i], status);
		} // if
	} // for

	// do a barrier
	std::cout << "rank " << rank << " reached barrier" << std::endl;
	p2p.request(P2PTag::barrier);
	p2p.sync();

	// do all reduce max	
	double val = static_cast<double>(rank);
	p2p.request(P2PTag::allreduce_max_double, P2PTag::data, 1);
	p2p.send(&val, 1);
	p2p.recv(&val, 1);
	std::cout << "all reduce max " << val << std::endl;

	// do all reduce sum	
	val = static_cast<double>(rank);
	p2p.request(P2PTag::allreduce_sum_double, P2PTag::data, 1);
	p2p.send(&val, 1);
	p2p.recv(&val, 1);
	std::cout << "all reduce sum " << val << std::endl;

	// do all gather
	//int * gather = new int[size*10];
	int * gather = (int *)malloc(size*10*sizeof(int));
	for(size_t i(0); i<10; i++) { gather[i] = 10*rank + i; }
	p2p.request(P2PTag::allgather_int, P2PTag::data, 10);
	p2p.send(gather, 10);
	p2p.recv(gather, size*10);

	if(rank == 0) {
		for(int r(0); r<size; r++) {
			std::cout << "gathered";

			for(size_t i(0); i<10; i++) {
				std::cout << " " << gather[r*10 + i]; 
			} // for
			std::cout << std::endl;
		} // for
	} // if

	// print results
	for(int i(0); i<size; i++) {
		if(i != rank) {
			std::cout << "rank " << rank << " got message " << msg[i].rank <<
				" " << msg[i].size << " from " << i << std::endl;
		} // if
	} // for

	/*
	int reason = 0;
	p2p.request(P2PTag::abort);
	p2p.send(&reason, 1);
	*/

	// finalize
	p2p.finalize();

	// cleanup
	delete[] sreq;
	delete[] rreq;
	delete[] msg;

	return 0;
} // main
