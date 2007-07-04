/*------------------------------------------------------------------------------
------------------------------------------------------------------------------*/

#include <sstream>
#include <string>
#include <P2PConnection.hxx>

int main(int argc, char *argv[]) {

	MPBuffer<float> send_buffer[max_buffers];
	MPBuffer<float> recv_buffer[max_buffers];
	MPRequest request;
	double dbuf[5];
	int ibuf_send[5];
	int * ibuf_recv(NULL);
	int64_t lbuf_send[5];
	int64_t * lbuf_recv(NULL);

	// get an instance of the p2p connection
	P2PConnection & p2p = P2PConnection::instance();

	// initialize everything
	p2p.init(argc, argv);


	int rank = p2p.global_id();
	int size = p2p.global_size();

	// initialize integer arrays
	ibuf_recv = new int[5*size];
	lbuf_recv = new int64_t[5*size];

	int count = 2*sizeof(float);

	for(size_t i(0); i<size; i++) {
		if(i != rank) {
			request.set(P2PTag::irecv, 0, count, i, i);
			send_buffer[i].data()[0] = static_cast<float>(rank);
			send_buffer[i].data()[1] = static_cast<float>(size);
			p2p.post(request);
			p2p.irecv(reinterpret_cast<char *>(recv_buffer[i].data()),
				request.count, request.tag, request.id);
		} // if
	} // for
	
	for(size_t i(0); i<size; i++) {
		if(i != rank) {
			request.set(P2PTag::send, 0, count, i, i);
			send_buffer[i].data()[0] = static_cast<float>(rank);
			send_buffer[i].data()[1] = static_cast<float>(size);

			// request that host do a blocking send
			// but do a non-blocking point-to-point send
			// so that we can keep working
			p2p.post(request);
			p2p.isend(reinterpret_cast<char *>(send_buffer[i].data()),
				request.count, request.tag, request.id);
		} // if
	} // for
	
	for(size_t i(0); i<size; i++) {
		if(i != rank) {
			p2p.wait_send(i);
		} // if
	} // for

	for(size_t i(0); i<size; i++) {
		if(i != rank) {
			request.set(P2PTag::wait_recv, 0, count, i, i);
			p2p.post(request);
			p2p.wait_recv(i);
		} // if
	} // for

	p2p.post(P2PTag::barrier);
	p2p.barrier();

	for(size_t i(0); i<size; i++) {
		if(i != rank) {
			std::cout << recv_buffer[i].data()[0] <<
				" " << recv_buffer[i].data()[1] << std::endl;
		} // if
	} // for

	// do all reduce max on doubles
	for(size_t i(0); i<5; i++) {
		dbuf[i] = static_cast<double>(rank*10 + i);
	} // for
	request.set(P2PTag::allreduce_max_double, 0, 5, rank);
	p2p.post(request);
	p2p.send(dbuf, request.count, request.tag);
	p2p.recv(dbuf, request.count, request.tag, request.id);

	p2p.post(P2PTag::barrier);
	p2p.barrier();

	{
		std::ostringstream ostr;
		ostr << dbuf[0];
		for(size_t i(1); i<5; i++) {
			ostr << " " << dbuf[i];
		} // for
		std::cout << "allreduce max double: " << ostr.str() << std::endl;
	} // scope

	// do all reduce sum on doubles
	request.set(P2PTag::allreduce_sum_double, 0, 5, rank);
	p2p.post(request);
	p2p.send(dbuf, request.count, request.tag);
	p2p.recv(dbuf, request.count, request.tag, request.id);

	p2p.post(P2PTag::barrier);
	p2p.barrier();

	{
		std::ostringstream ostr;
		ostr << dbuf[0];
		for(size_t i(1); i<5; i++) {
			ostr << " " << dbuf[i];
		} // for
		std::cout << "allreduce sum double: " << ostr.str() << std::endl;
	} // scope

	// do all reduce sum on integers
	for(size_t i(0); i<5; i++) {
		ibuf_send[i] = rank*10 + i;
	} // for
	request.set(P2PTag::allreduce_sum_int, 0, 5, rank);
	p2p.post(request);
	p2p.send(ibuf_send, request.count, request.tag);
	p2p.recv(ibuf_recv, request.count, request.tag, request.id);

	p2p.post(P2PTag::barrier);
	p2p.barrier();

	{
		std::ostringstream ostr;
		ostr << ibuf_recv[0];
		for(size_t i(1); i<5; i++) {
			ostr << " " << ibuf_recv[i];
		} // for
		std::cout << "allreduce max integer: " << ostr.str() << std::endl;
	} // scope

	// do all gather on integers
	for(size_t i(0); i<5; i++) {
		ibuf_send[i] = rank*10 + i;
	} // for
	request.set(P2PTag::allgather_int, 0, 5, rank);
	p2p.post(request);
	p2p.send(ibuf_send, request.count, request.tag);
	p2p.recv(ibuf_recv, request.count*size, request.tag, request.id);

	p2p.post(P2PTag::barrier);
	p2p.barrier();

	{
		std::ostringstream ostr;
		ostr << ibuf_recv[0];
		for(size_t i(1); i<5*size; i++) {
			ostr << " " << ibuf_recv[i];
		} // for
		std::cout << "allgather integer: " << ostr.str() << std::endl;
	} // scope

	// do all gather on int64_t
	for(size_t i(0); i<5; i++) {
		lbuf_send[i] = rank*10 + i;
	} // for
	request.set(P2PTag::allgather_int64, 0, 5, rank);
	p2p.post(request);
	p2p.send(lbuf_send, request.count, request.tag);
	p2p.recv(lbuf_recv, request.count*size, request.tag, request.id);

	p2p.post(P2PTag::barrier);
	p2p.barrier();

	{
		std::ostringstream ostr;
		ostr << lbuf_recv[0];
		for(size_t i(1); i<5*size; i++) {
			ostr << " " << lbuf_recv[i];
		} // for
		std::cout << "allgather int64: " << ostr.str() << std::endl;
	} // scope

/*
	// test abort
	int reason = 0;
	p2p.post(P2PTag::abort);
	p2p.send(&reason, 1, 0);
*/

	// ask relay to stop
	p2p.post(P2PTag::end);

	// finalize communication
	p2p.finalize();

	delete[] ibuf_recv;
	delete[] lbuf_recv;

	return 0;

} // main
