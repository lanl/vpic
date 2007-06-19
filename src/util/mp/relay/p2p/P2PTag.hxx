/*
	Definition of P2PTag data structure

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef P2PTag_hxx
#define P2PTag_hxx

struct P2PTag {

	enum {
		end,
		data,
		send,
		recv,
		isend,
		irecv,
		wait_send,
		wait_recv,
		abort,
		barrier,
		allreduce_max_double,
		allreduce_sum_double,
		allreduce_sum_int,
		allgather_int,
		allgather_int64,
		io_open_read,
		io_open_write,
		io_open_write_append,
		io_write,
		io_read,
		io_close
	}; // enum

}; // class P2PTag

#endif // P2PTag_hxx
