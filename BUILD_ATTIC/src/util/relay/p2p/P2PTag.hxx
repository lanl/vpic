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
		request,
		pending,
		data,
		send,
		recv,
		isend,
		irecv,
		wait_send,
		wait_recv,
		abort,
		barrier,
		wtime,
		allreduce_max_double,
		allreduce_sum_double,
		allreduce_sum_int,
		gather_uc,
		allgather_int,
		allgather_int64,
		io_open_read,
		io_open_read_write,
		io_open_write,
		io_open_write_read,
		io_open_append,
		io_open_append_read,
		io_write,
		io_read,
		io_seek,
		io_tell,
		io_rewind,
		io_size,
		io_close,
		utils_mkdir,
		utils_getcwd
	}; // enum

}; // class P2PTag

#endif // P2PTag_hxx
