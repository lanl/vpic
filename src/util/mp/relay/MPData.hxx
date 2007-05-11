/*
	Definition of MPData class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$
	vim: set ts=3 :
*/

#ifndef MPData_h
#define MPData_h

// max number of buffers that can be handled
static const int max_buffers(27);

// default buffer size
static const int default_buffer_size(4096);

// state of a relay request
enum MPState {
	pending,
	complete,
	error
}; // MPState

#define MP_HOST 0
#define MP_ACCEL 1

// message passing request
template<int ROLE> struct MPRequest_T {};

template<> struct MPRequest_T<MP_HOST> {

	MPRequest_T(int tag_ = 0, int count_ = 0, int id_ = 0, int peer_ = 0,
		MPState state_ = pending)
		: tag(tag_), count(count_), id(id_),
		peer(peer_), state(state_) {}

	void set(int tag_ = 0, int count_ = 0, int id_ = 0,
		int peer_ = 0, MPState state_ = pending) {
		tag = tag_;
		count = count_;
		id = id_;
		peer = peer_;
		state = state_;
	} // set

	int tag;
	int count;
	int id;
	int peer;
	MPState state;

}; // struct MPRequest_T

template<> struct MPRequest_T<MP_ACCEL> {

	MPRequest_T(int tag_ = 0, int count_ = 0, int id_ = 0, int peer_ = 0)
		: tag(tag_), count(count_), id(id_), peer(peer_) {}

	void set(int tag_ = 0, int count_ = 0, int id_ = 0, int peer_ = 0) {
		tag = tag_;
		count = count_;
		id = id_;
		peer = peer_;
	} // set

	int tag;
	int count;
	int id;
	int peer;

}; // struct MPRequest_T

#if defined HOST_BUILD
	typedef MPRequest_T<MP_HOST> MPRequest;
#else
	typedef MPRequest_T<MP_ACCEL> MPRequest;
#endif // BUILD

// this is a kluge that avoids defining an MPI_TYPE for
// the MPRequest_T data structure
int request_count() { return sizeof(MPRequest_T<MP_ACCEL>)/sizeof(int); }

// message passing buffer
template<typename T> class MPBuffer
	{
	public:

		MPBuffer()
			: count_(default_buffer_size),
			data_(new T[default_buffer_size]) {}
		~MPBuffer()
			{ count_ = 0; delete[] data_; }

		int count() { return count_; }
		T * data() { return data_; }

		void resize(int count) {
			if(count > count_) {
				T * tmp = data_;
				data_ = new T[count];
				memcpy(data_, tmp, count_*sizeof(T));
				count_ = count;
				delete[] tmp;
			} // if
		} // resize

	private:

		int count_;
		T * data_;

}; // class MPBuffer

#endif // MPData_h
