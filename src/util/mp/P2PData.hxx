/*
	Definition of P2PData data structures

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef P2PData_hxx
#define P2PData_hxx

struct P2PTag {
	enum {
		end,
		data,
		send,
		recv,
		isend,
		irecv,
		wait,
		abort,
		barrier,
		allreduce_max_double,
		allreduce_sum_double,
		allreduce_sum_int,
		allgather_int,
		allgather_int64
	}; // enum
}; // class P2PTag

struct P2PHeader {
	P2PHeader(int tag_ = P2PTag::data, int count_ = 0, int rank_ = 0)
		: tag(tag_), count(count_), rank(rank_) {}

	int tag;
	int count;
	int rank;
}; // struct P2PHeader

#endif // P2PData_hxx
