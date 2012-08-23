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

//static const uint64_t default_buffer_size(4096);
static const uint64_t default_buffer_size(51200);

static const uint64_t filename_size(256);

static const uint64_t io_line_size(1024);

// this should change to a -D compile option
static const uint64_t io_buffer_size(51200);
static const uint64_t utils_buffer_size(4096);
//static const uint64_t io_buffer_size(2048);
//static const uint64_t io_buffer_size(1024);
//static const uint64_t io_buffer_size(64);
//static const uint64_t io_buffer_size(32);
//static const uint64_t io_buffer_size(16);
//static const uint64_t io_buffer_size(8);

// state of a relay request
enum MPState {
	unitialized,
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
		int id_ = 0, int peer_ = 0, MPState state_ = unitialized)
		: p2ptag(p2ptag_), tag(tag_), count(count_), id(id_),
		peer(peer_), state(state_) {}

	void set(int p2ptag_ = 0, int tag_ = 0, int count_ = 0,
		int id_ = 0, int peer_ = 0, MPState state_ = unitialized) {
		p2ptag = p2ptag_;
		tag = tag_;
		count = count_;
		id = id_;
		peer = peer_;
		state = state_;
	} // set

	/*
	int p2ptag __attribute__ ((aligned (16)));
	int tag __attribute__ ((aligned (16)));
	int count __attribute__ ((aligned (16)));
	int id __attribute__ ((aligned (16)));
	int peer __attribute__ ((aligned (16)));
	MPState state __attribute__ ((aligned (16)));
	*/

	int p2ptag;
	int tag;
	int count;
	int id;
	int peer;
	MPState state;

} __attribute__ ((aligned (16))); // struct MPRequest_T

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

	/*
	int p2ptag __attribute__ ((aligned (16)));
	int tag __attribute__ ((aligned (16)));
	int count __attribute__ ((aligned (16)));
	int id __attribute__ ((aligned (16)));
	int peer __attribute__ ((aligned (16)));
	*/
	int p2ptag;
	int tag;
	int count;
	int id;
	int peer;

} __attribute__ ((aligned (16))); // struct MPRequest_T

#if defined HOST_BUILD
	typedef MPRequest_T<MP_HOST> MPRequest;
#else
	typedef MPRequest_T<MP_ACCEL> MPRequest;
#endif // BUILD

// print request information to stdout
template<int ROLE>
std::ostream & operator << (std::ostream & stream,
	const MPRequest_T<ROLE> & request)
	{
		stream <<
			"p2ptag: " << request.p2ptag <<
			" tag: " << request.tag <<
			" count: " << request.count <<
			" id: " << request.id <<
			" peer: " << request.peer << std::endl;

		return stream;
	} // operator <<

// this is a kluge that avoids defining an MPI_TYPE for
// the MPRequest_T data structure
inline int request_count() {
#if defined USE_DACS_P2P
	// DaCS needs the size in bytes
	return sizeof(MPRequest_T<MP_ACCEL>);
#else
	// MPI needs the size in elements
	return sizeof(MPRequest_T<MP_ACCEL>)/sizeof(int);
#endif // USE_DACS_P2P
} // request_count

#include <unistd.h>
#include <sys/syscall.h>

// message passing buffer
template<typename T, int BS = default_buffer_size, int ROLE = 0> class MPBuffer
	{
	public:

		MPBuffer() : id_(0), size_(BS), data_(new T[BS]) {
			/*
			dataOrig_ = data_;
			int pid = syscall(SYS_gettid);
			printf("Role %d(%d - %p - %d): called constructor with data: %p\n",
				ROLE, pid, (void *)this, id_, (void *)data_);
			*/
		}

		~MPBuffer() {
			/*
			int pid = syscall(SYS_gettid);
			if(data_ != dataOrig_) {
			printf("ERROR!!! Role %d(%d - %p - %d): called destructor with data: %p\n",
				ROLE, pid, (void *)this, id_, (void *)data_);
			}
			else {
			printf("Role %d(%d - %p - %d): called destructor with data: %p\n",
				ROLE, pid, (void *)this, id_, (void *)data_);
			} // else
			*/
			delete[] data_;
		}

		void setId(uint32_t id) { id_ = id; }

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

		uint32_t id_;
		int size_;
		T * data_;
		T * dataOrig_;

}; // class MPBuffer

#endif // MPData_hxx
