/*
        Definition of P2PIOPolicy class

        Author: Benjamin Karl Bergen

        $Revision$
        $LastChangedBy$
        $LastChangedDate$
        vim: set ts=3 :
*/

#ifndef P2PIOPolicy_h
#define P2PIOPolicy_h

#include <cstdarg>
#include <string>

#include "FileIOData.h"
#include "MPData.h"
#include "P2PConnection.h"
#include "swap.h"

/*!
        \class P2PIOPolicy P2PIOPolicy.h
        \brief  provides...
*/
template <bool swapped>
class P2PIOPolicy
{
  public:
    //! Constructor
    P2PIOPolicy()
        : id_( -1 )
        , mode_( io_closed )
        , current_( 0 )
        , write_buffer_offset_( 0 )
        , read_buffer_offset_( 0 )
        , file_size_( 0 )
    {
        pending_[0] = false;
        pending_[1] = false;
        buffer_fill_[0] = 0;
        buffer_fill_[1] = 0;

        io_line_.setId( 1 );
    }

    //! Destructor
    ~P2PIOPolicy() {}

    FileIOStatus open( const char *filename, FileIOMode mode );
    int32_t close();

    bool isOpen() { return ( id_ >= 0 ); }

    int64_t size();

    void print( const char *format, va_list &args );

    template <typename T>
    size_t read( T *data, size_t elements );
    template <typename T>
    size_t write( const T *data, size_t elements );

    int64_t seek( uint64_t offset, int32_t whence );
    int64_t tell();
    void rewind();
    void flush();

  private:
    template <typename T>
    inline void swap_bytes( T *data, size_t elements );

    void send_write_block( uint32_t buffer );
    int64_t wait_write_block( uint32_t buffer );

    void request_read_block( uint32_t buffer );
    int64_t wait_read_block( uint32_t buffer );

    MPBuffer<char, io_buffer_size, 1> io_buffer_[2];
    MPBuffer<char, io_line_size, 1> io_line_;

    int32_t id_;
    FileIOMode mode_;
    uint32_t current_;
    uint64_t write_buffer_offset_;
    uint64_t read_buffer_offset_;
    uint64_t buffer_fill_[2];
    bool pending_[2];
    int request_id[2];
    MPRequest request_[2];

    int64_t file_size_;
    ldiv_t read_blocks_;

}; // class P2PIOPolicy

template <bool swapped>
FileIOStatus P2PIOPolicy<swapped>::open( const char *filename, FileIOMode mode )
{
    assert( id_ < 0 );
    P2PConnection &p2p = P2PConnection::instance();

    /*
    std::cerr << "PPE rank: " << p2p.global_id() <<
            " opening " << filename << std::endl;
    */

    // this sends the string terminator as well
    size_t msg_size = strlen( filename ) + 1;
    MPRequest request;

    // re-initialize some values
    write_buffer_offset_ = 0;
    read_buffer_offset_ = 0;
    pending_[0] = false;
    pending_[1] = false;
    buffer_fill_[0] = 0;
    buffer_fill_[1] = 0;
    current_ = 0;

    // file io mode
    switch ( mode )
    {

    case io_read:
        request.set( P2PTag::io_open_read, P2PTag::data, msg_size );
        p2p.post( request );
        break;

    case io_read_write:
        request.set( P2PTag::io_open_read_write, P2PTag::data, msg_size );
        p2p.post( request );
        break;

    case io_write:
        request.set( P2PTag::io_open_write, P2PTag::data, msg_size );
        p2p.post( request );
        break;

    case io_write_read:
        request.set( P2PTag::io_open_write_read, P2PTag::data, msg_size );
        p2p.post( request );
        break;

    case io_append:
        request.set( P2PTag::io_open_append, P2PTag::data, msg_size );
        p2p.post( request );
        break;

    case io_append_read:
        request.set( P2PTag::io_open_append_read, P2PTag::data, msg_size );
        p2p.post( request );
        break;

    default:
        return fail;

    } // switch

    // save this for flush logic
    mode_ = mode;

    // send the filename to peer
    p2p.send( const_cast<char *>( filename ), request.count, request.tag );

    // get file descriptor id
    p2p.recv( &id_, 1, request.tag, request.id );
    assert( id_ >= 0 );

    if ( mode_ == io_read || mode_ == io_read_write )
    {
        p2p.recv( &file_size_, 1, request.tag, id_ );
        read_blocks_ = ldiv( file_size_, io_buffer_size );

        /*
        std::cerr << "PPE rank: " << p2p.global_id() <<
                " file size: " << file_size_ << " read blocks: " <<
                read_blocks_.quot << std::endl;
        */

        // request block
        request_read_block( current_ );

        // request next block
        // This will only do a real request if
        // read_blocks_ > 1
        request_read_block( current_ ^ 1 );

        // wait on the first block
        wait_read_block( current_ );
    } // if

    return ok;
} // P2PIOPolicy<>::open

template <bool swapped>
int32_t P2PIOPolicy<swapped>::close()
{

    P2PConnection &p2p = P2PConnection::instance();
    /*
    std::cerr << "PPE rank: " << p2p.global_id() <<
            " closing file" << std::endl;
    */

    // force write if current block hasn't been written
    if ( write_buffer_offset_ > 0 )
    {
        flush();
    } // if

    /*
    std::cerr << "PPE rank: " << p2p.global_id() <<
            " after flush" << std::endl;
    */

    MPRequest request( P2PTag::io_close, P2PTag::data, 1, id_ );
    p2p.post( request );

    int32_t status( 0 );
    p2p.recv( &status, 1, request.tag, request.id );

    id_ = -1;

    return status;
} // P2PIOPolicy<>::close

template <bool swapped>
int64_t P2PIOPolicy<swapped>::size()
{
    assert( id_ >= 0 );

    P2PConnection &p2p = P2PConnection::instance();

    MPRequest request( P2PTag::io_size, P2PTag::data, 1, id_ );
    P2PConnection::instance().post( request );

    int64_t tmp;
    p2p.recv( &tmp, 1, request.tag, request.id );

    return tmp;
} // P2PIOPolicy<>::size

template <bool swapped>
void P2PIOPolicy<swapped>::print( const char *format, va_list &args )
{
    assert( id_ >= 0 );

    // sprintf to local buffer
    vsprintf( io_line_.data(), format, args );

    /*
    P2PConnection & p2p = P2PConnection::instance();
    std::cerr << "PPE rank: " << p2p.global_id() <<
            " printing " << io_line_.data() << std::endl;
    */

    // use write function to do actual work
    P2PIOPolicy::write( io_line_.data(), strlen( io_line_.data() ) );
} // P2PIOPolicy<>::print

template <bool swapped>
template <typename T>
size_t P2PIOPolicy<swapped>::read( T *data, size_t elements )
{
    assert( id_ >= 0 );

    // everything is done in bytes
    uint64_t bytes = elements * sizeof( T );
    char *bdata = reinterpret_cast<char *>( data );
    uint64_t bdata_offset( 0 );
    uint64_t read_bytes( 0 );

    do
    {
        const int64_t over_run =
            ( read_buffer_offset_ + bytes ) - buffer_fill_[current_];

        if ( over_run > 0 )
        {
            const uint64_t under_run =
                buffer_fill_[current_] - read_buffer_offset_;

            // copy remainder of current buffer to data
            memcpy( bdata + bdata_offset,
                    io_buffer_[current_].data() + read_buffer_offset_,
                    under_run );
            bdata_offset += under_run;
            bytes -= under_run;

            // re-fill current buffer
            request_read_block( current_ );
            current_ ^= 1;
            if ( pending_[current_] )
            {
                read_bytes += wait_read_block( current_ );
            } // if
            read_buffer_offset_ = 0;
        }
        else
        {
            memcpy( bdata + bdata_offset,
                    io_buffer_[current_].data() + read_buffer_offset_, bytes );
            read_buffer_offset_ += bytes;
            bytes = 0;
        } // if
    } while ( bytes > 0 );

    // this will only do something if
    // this class was instantiated as P2PIOPolicy<true>
    swap_bytes( data, elements );

    return static_cast<size_t>( read_bytes );
} // P2PIOPolicy<>::read

template <bool swapped>
template <typename T>
size_t P2PIOPolicy<swapped>::write( const T *data, size_t elements )
{
    assert( id_ >= 0 );
    assert( mode_ != io_read );

    // book-keeping is done in bytes
    uint64_t bytes( elements * sizeof( T ) );
    const char *bdata = reinterpret_cast<const char *>( data );
    uint64_t bdata_offset( 0 );
    uint64_t write_bytes( 0 );

    /*
    P2PConnection & p2p = P2PConnection::instance();
    std::cerr << "PPE rank: " << p2p.global_id() <<
            " bytes " << bytes << std::endl;
    */
    do
    {
        const int64_t over_run =
            ( write_buffer_offset_ + bytes ) - io_buffer_[current_].size();

        /*
        std::cerr << "PPE rank: " << p2p.global_id() <<
                " over_run " << over_run << std::endl;
        */
        if ( over_run > 0 )
        {
            const uint64_t under_run =
                io_buffer_[current_].size() - write_buffer_offset_;

            // because of the possiblity of byte swapping
            // we need to make sure that only even multiples
            // of the type are copied at once
            const uint64_t copy_bytes =
                ( under_run / sizeof( T ) ) * sizeof( T );

            /*
            printf("PPE rank: %d dst %p src %p bytes %ld\n",
                    p2p.global_id(),
                    io_buffer_[current_].data() + write_buffer_offset_,
                    bdata + bdata_offset, copy_bytes);
            */

            memcpy( io_buffer_[current_].data() + write_buffer_offset_,
                    bdata + bdata_offset, copy_bytes );

            // need to force type here to get correct swapping
            swap_bytes<T>( reinterpret_cast<T *>( io_buffer_[current_].data() +
                                                  write_buffer_offset_ ),
                           copy_bytes / sizeof( T ) );

            bdata_offset += copy_bytes;
            bytes -= copy_bytes;

            request_[current_].set( P2PTag::io_write, P2PTag::data,
                                    write_buffer_offset_ + copy_bytes, id_ );
            send_write_block( current_ );
            current_ ^= 1;
            if ( pending_[current_] )
            {
                write_bytes += wait_write_block( current_ );
            } // if
            write_buffer_offset_ = 0;
        }
        else
        {
            /*
            printf("PPE rank: %d dst %p src %p bytes %ld\n",
                    p2p.global_id(),
                    io_buffer_[current_].data() + write_buffer_offset_,
                    bdata + bdata_offset, bytes);
            */

            memcpy( io_buffer_[current_].data() + write_buffer_offset_,
                    bdata + bdata_offset, bytes );

            // need to force type here to get correct swapping
            swap_bytes( reinterpret_cast<T *>( io_buffer_[current_].data() +
                                               write_buffer_offset_ ),
                        bytes / sizeof( T ) );

            write_buffer_offset_ += bytes;
            bytes = 0;
        } // if
    } while ( bytes > 0 );

    return static_cast<size_t>( write_bytes );
} // P2PIOPolicy<>::write

template <bool swapped>
int64_t P2PIOPolicy<swapped>::seek( uint64_t offset, int32_t whence )
{
    assert( id_ >= 0 );

    MPRequest request( P2PTag::io_seek, P2PTag::data, 0, id_ );
    P2PConnection &p2p = P2PConnection::instance();

    p2p.post( request );
    p2p.send( &offset, 1, P2PTag::data );
    p2p.send( &whence, 1, P2PTag::data );

    // FIXME: need real return
    return 0;
} // P2PIOPolicy<>::seek

template <bool swapped>
int64_t P2PIOPolicy<swapped>::tell()
{
    assert( id_ >= 0 );

    MPRequest request( P2PTag::io_tell, P2PTag::data, 0, id_ );
    P2PConnection &p2p = P2PConnection::instance();

    p2p.post( request );

    int64_t offset;
    p2p.recv( &offset, 1, request.tag, request.id );

    return offset;
} // P2PIOPolicy<>::tell

template <bool swapped>
void P2PIOPolicy<swapped>::rewind()
{
    assert( id_ >= 0 );

    MPRequest request( P2PTag::io_rewind, P2PTag::data, 0, id_ );
    P2PConnection &p2p = P2PConnection::instance();

    p2p.post( request );
} // P2PIOPolicy<>::rewind

template <bool swapped>
void P2PIOPolicy<swapped>::request_read_block( uint32_t buffer )
{
    P2PConnection &p2p = P2PConnection::instance();

    if ( read_blocks_.quot > 0 )
    {
        request_[buffer].set( P2PTag::io_read, P2PTag::data,
                              io_buffer_[buffer].size(), id_ );

        /*
        std::cerr << "PPE rank: " << p2p.global_id() <<
                " requesting " << request_[buffer].count <<
                " bytes " << std::endl;
        */

        p2p.post( request_[buffer] );
        p2p.irecv( io_buffer_[buffer].data(), request_[buffer].count,
                   request_[buffer].tag, request_[buffer].id );
        pending_[buffer] = true;
        buffer_fill_[buffer] = request_[buffer].count;

        read_blocks_.quot--;
    }
    else if ( read_blocks_.rem > 0 )
    {
        request_[buffer].set( P2PTag::io_read, P2PTag::data, read_blocks_.rem,
                              id_ );

        /*
        std::cerr << "PPE rank: " << p2p.global_id() <<
                " requesting " << request_[buffer].count <<
                " bytes " << std::endl;
        */

        p2p.post( request_[buffer] );
        p2p.irecv( io_buffer_[buffer].data(), request_[buffer].count,
                   request_[buffer].tag, request_[buffer].id );
        pending_[buffer] = true;
        buffer_fill_[buffer] = request_[buffer].count;

        read_blocks_.rem = 0;
    } // if
} // P2PIOPolicy<>::request_read_block

template <bool swapped>
int64_t P2PIOPolicy<swapped>::wait_read_block( uint32_t buffer )
{
    P2PConnection &p2p = P2PConnection::instance();
    p2p.wait_recv( request_[buffer].id );
    pending_[buffer] = false;
    return request_[buffer].count;
} // P2PIOPolicy<>::wait_read_block

template <bool swapped>
void P2PIOPolicy<swapped>::send_write_block( uint32_t buffer )
{
    P2PConnection &p2p = P2PConnection::instance();

    p2p.post( request_[buffer] );
    /*
    std::cerr << "PPE rank: " << p2p.global_id() <<
            " sending " << request_[buffer].count << std::endl;
    */
    p2p.isend( io_buffer_[buffer].data(), request_[buffer].count,
               request_[buffer].tag, request_[buffer].id );
    pending_[buffer] = true;
} // P2PIOPolicy<>::send_write_block

template <bool swapped>
int64_t P2PIOPolicy<swapped>::wait_write_block( uint32_t buffer )
{
    P2PConnection &p2p = P2PConnection::instance();
    p2p.wait_send( request_[buffer].id );
    pending_[buffer] = false;
    return request_[buffer].count;
} // P2PIOPolicy<>::wait_write_block

template <bool swapped>
void P2PIOPolicy<swapped>::flush()
{
    /*
    P2PConnection & p2p = P2PConnection::instance();
    std::cerr << "PPE rank: " << p2p.global_id() <<
            " write_buffer_offset_: " << write_buffer_offset_ << std::endl;
    */

    // check to see if we need to flush the current buffer
    if ( write_buffer_offset_ > 0 )
    {
        request_[current_].set( P2PTag::io_write, P2PTag::data,
                                write_buffer_offset_, id_ );
        send_write_block( current_ );

        current_ ^= 1;
    } // if

    // wait on the other buffer if it is in flight
    if ( pending_[current_] )
    {
        wait_write_block( current_ );
    } // if

    // wait on the one we just sent
    current_ ^= 1;
    wait_write_block( current_ );
} // P2PIOPolicy<>::flush

template <>
template <typename T>
inline void P2PIOPolicy<true>::swap_bytes( T *data, size_t elements )
{
    for ( size_t i( 0 ); i < elements; i++ )
    {
        utils::swap( data[i] );
    } // for
} // P2PIOPolicy<>::swap_bytes

template <>
template <typename T>
inline void P2PIOPolicy<false>::swap_bytes( T *data, size_t elements )
{
} // P2PIOPolicy<>::swap_bytes

#endif // P2PIOPolicy_h
