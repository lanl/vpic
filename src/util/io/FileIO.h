/*
        Definition of FileIO class

        Author: Benjamin Karl Bergen

        $Revision$
        $LastChangedBy$
        $LastChangedDate$
        vim: set ts=3 :
*/

#ifndef FileIO_h
#define FileIO_h

#include "FileIOData.h"
#include <stdarg.h>

/*!
        \class FileIO FileIO.h
        \brief  provides...
*/
template <class ReadWritePolicy>
struct FileIO_T : public ReadWritePolicy
{
    //! Constructor
    FileIO_T() {}

    //! Destructor
    ~FileIO_T() {}

    FileIOStatus open( const char* filename, FileIOMode mode )
    {
        return ReadWritePolicy::open( filename, mode );
    }
    int32_t close() { return ReadWritePolicy::close(); }

    bool isOpen() { return ReadWritePolicy::isOpen(); }

    int64_t size() { return ReadWritePolicy::size(); }

    void print( const char* format, ... )
    {
        va_list args;
        va_start( args, format );
        ReadWritePolicy::print( format, args );
    }

    template <typename T>
    size_t read( T* data, size_t elements )
    {
        return ReadWritePolicy::read( data, elements );
    }
    template <typename T>
    size_t write( const T* data, size_t elements )
    {
        return ReadWritePolicy::write( data, elements );
    }

    int64_t seek( uint64_t offset, int32_t whence )
    {
        return ReadWritePolicy::seek( offset, whence );
    }
    int64_t tell() { return ReadWritePolicy::tell(); }
    void rewind() { ReadWritePolicy::rewind(); }

}; // class FileIO_T

#if defined USE_MPRELAY

#if defined HOST_BUILD
#include "StandardIOPolicy.h"

typedef FileIO_T<StandardIOPolicy> FileIO;
#else
#include "P2PIOPolicy.h"

// typedef FileIO_T<P2PIOPolicy<true> > FileIOSwapped;
// typedef FileIO_T<P2PIOPolicy<true> > FileIO;
typedef FileIO_T<P2PIOPolicy<true>> FileIO;
typedef FileIO_T<P2PIOPolicy<false>> FileIOUnswapped;
#endif // BUILD

#else
#include "StandardIOPolicy.h"
typedef FileIO_T<StandardIOPolicy> FileIO;
typedef FileIO_T<StandardIOPolicy> FileIOUnswapped;
#endif // MP Implementation

#endif // FileIO_h
