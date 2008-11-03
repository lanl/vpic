#ifndef P2PUtilsPolicy_hxx
#define P2PUtilsPolicy_hxx

#include <P2PConnection.hxx>

class P2PUtilsPolicy
	{
	public:

		P2PUtilsPolicy() {}
		~P2PUtilsPolicy() {}

		static int makeDirectory(const char * dirname);

	private:

	}; // class P2PUtilsPolicy

inline int P2PUtilsPolicy::makeDirectory(const char * dirname)
	{
		P2PConnection & p2p = P2PConnection::instance();

		size_t msg_size = strlen(dirname)+1;
		int retval;
		MPRequest request(P2PTag::utils_mkdir, P2PTag::data, msg_size);

		p2p.post(request);
		p2p.send(const_cast<char *>(dirname), request.count, request.tag);
		p2p.recv(&retval, 1, request.tag, request.id);

		return retval;
	} // P2PUtilsPolicy::makeDirectory

#endif // P2PUtilsPolicy_hxx
