#ifndef P2PUtilsPolicy_h
#define P2PUtilsPolicy_h

#include "P2PConnection.h"

class P2PUtilsPolicy
{
  public:
    P2PUtilsPolicy() {}
    ~P2PUtilsPolicy() {}

    static int makeDirectory( const char *dirname );
    static int getCurrentWorkingDirectory( char *dirname, size_t size );

  private:
}; // class P2PUtilsPolicy

inline int P2PUtilsPolicy::makeDirectory( const char *dirname )
{
    P2PConnection &p2p = P2PConnection::instance();

    size_t msg_size = strlen( dirname ) + 1;
    int retval;
    MPRequest request( P2PTag::utils_mkdir, P2PTag::data, msg_size );

    p2p.post( request );
    p2p.send( const_cast<char *>( dirname ), request.count, request.tag );
    p2p.recv( &retval, 1, request.tag, request.id );

    return retval;
} // P2PUtilsPolicy::makeDirectory

inline int P2PUtilsPolicy::getCurrentWorkingDirectory( char *dirname,
                                                       size_t size )
{
    P2PConnection &p2p = P2PConnection::instance();

    size_t msg_size = size;
    int retval;
    MPRequest request( P2PTag::utils_mkdir, P2PTag::data, msg_size );

    p2p.post( request );
    p2p.recv( dirname, size, request.count, request.tag );
    p2p.recv( &retval, 1, request.tag, request.id );

    return retval;
} // P2PUtilsPolicy::getCurrentWorkingDirectory

#endif // P2PUtilsPolicy_h
