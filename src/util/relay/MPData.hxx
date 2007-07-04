/*
	Definition of MPData class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$
	vim: set ts=3 :
*/

#ifndef MPData_hxx
#define MPData_hxx

static const uint64_t max_buffers(27);

static const uint64_t default_buffer_size(4096);

static const uint64_t filename_size(256);

static const uint64_t io_line_size(256);

// this should change to a -D compile option
static const uint64_t io_buffer_size(4096);
//static const uint64_t io_buffer_size(2048);
//static const uint64_t io_buffer_size(1024);
//static const uint64_t io_buffer_size(64);
//static const uint64_t io_buffer_size(32);
//static const uint64_t io_buffer_size(16);
//static const uint64_t io_buffer_size(8);

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

	MPRequest_T(int p2ptag_ = 0, int tag_ = 0, int count_ = 0,
		int id_ = 0, int peer_ = 0, MPState state_ = pending)
		: p2ptag(p2ptag_), tag(tag_), count(count_), id(id_),
		peer(peer_), state(state_) {}

	void set(int p2ptag_ = 0, int tag_ = 0, int count_ = 0,
		int id_ = 0, int peer_ = 0, MPState state_ = pending) {
		p2ptag = p2ptag_;
		tag = tag_;
		count = count_;
		id = id_;
		peer = peer_;
		state = state_;
	} // set

	int p2ptag;
	int tag;
	int count;
	int id;
	int peer;
	MPState state;

}; // struct MPRequest_T

template<> struct MPRequest_T<MP_ACCEL> {

	MPRequest_T(int p2ptag_ = 0, int tag_ = 0, int count_ = 0,
		int id_ = 0, int peer_ = 0)
		: p2ptag(p2ptag_), tag(tag_), count(count_), id(id_), peer(peer_) {}

	void set(int p2ptag_ = 0, int tag_ = 0, int count_ = 0,
		int id_ = 0, int peer_ = 0) {
		p2ptag = p2ptag_;
		tag = tag_;
		count = count_;
		id = id_;
		peer = peer_;
	} // set

	int p2ptag;
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
inline int request_count() {
	return sizeof(MPRequest_T<MP_ACCEL>)/sizeof(int);
} // request_count

// message passing buffer
template<typename T, int BS = default_buffer_size> class MPBuffer
	{
	public:

		MPBuffer()
			: size_(BS),
			data_(new T[BS]) {}
		~MPBuffer()
			{ size_ = 0; delete[] data_; }

		int size() { return size_; }
		T * data() { return data_; }

		void resize(int size) {
			if(size > size_) {
				T * tmp = data_;
				data_ = new T[size];
				memcpy(data_, tmp, size_*sizeof(T));
				size_ = size;
				delete[] tmp;
			} // if
		} // resize

	private:

		int size_;
		T * data_;

}; // class MPBuffer

#endif // MPData_hxx
