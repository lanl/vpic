#ifndef CheckPtIO_h
#define CheckPtIO_h

#include "../io/FileIO.h"
#include "checkpt_private.h"

struct CheckPtIO
{

    static checkpt_t* checkpt_open_rdonly( const char* name )
    {
        if ( !name )
            ERROR( ( "NULL name" ) );

        FileIO* fileIO = new FileIO;

        if ( fileIO->open( name, io_read ) != ok )
        {
            ERROR( ( "Unable to open \"%s\" for checkpt read", name ) );
        } // if

        return reinterpret_cast<checkpt_t*>( fileIO );
    } // checkpt_open_rdonly

    static checkpt_t* checkpt_open_wronly( const char* name )
    {
        if ( !name )
            ERROR( ( "NULL name" ) );

        FileIO* fileIO = new FileIO;

        if ( fileIO->open( name, io_write ) != ok )
        {
            ERROR( ( "Unable to open \"%s\" for checkpt read", name ) );
        } // if

        return reinterpret_cast<checkpt_t*>( fileIO );
    } // checkpt_open_wronly

    static void checkpt_close( checkpt_t* checkpt )
    {
        FileIO* fileIO = reinterpret_cast<FileIO*>( checkpt );

        int32_t err = fileIO->close();

        if ( err != 0 )
        {
            ERROR( ( "Error closing file (%d)", err ) );
        } // if

        delete fileIO;
    } // checkpt_close

    static void checkpt_read( checkpt_t* checkpt, void* data, size_t sz )
    {
        if ( !sz )
            return;
        if ( !checkpt || !data )
            ERROR( ( "Invalid checkpt_read request" ) );

        FileIO* fileIO = reinterpret_cast<FileIO*>( checkpt );

        // FIXME: add return values
        fileIO->read( reinterpret_cast<char*>( data ), sz );
    } // checkpt_read

    static void checkpt_write( checkpt_t* checkpt, const void* data, size_t sz )
    {
        if ( !sz )
            return;
        if ( !checkpt || !data )
            ERROR( ( "Invalid checkpt_read request" ) );

        FileIO* fileIO = reinterpret_cast<FileIO*>( checkpt );

        // FIXME: add return values
        fileIO->write( reinterpret_cast<const char*>( data ), sz );
    } // checkpt_write

}; // struct CheckPtIO

#endif // CheckPtIO_h
