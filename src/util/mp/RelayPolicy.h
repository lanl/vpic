#ifndef RelayPolicy_h
#define RelayPolicy_h

#include "ConnectionManager.h"
#include "P2PConnection.h"
#include "checkpt.h"

/* Define this comm and mp opaque handles */
/* FIXME: PARENT, COLOR AND KEY ARE FOR FUTURE EXPANSION */

struct collective
{
    collective_t *parent;
    int color;
    int key;
    /* FIXME: RELAY HAS NO WAY OF DOING AN MPI_COMM_SPLIT RIGHT NOW */
};

struct mp
{
    int n_port;
    char *ALIGNED( 128 ) * rbuf;
    char *ALIGNED( 128 ) * sbuf;
    int *rbuf_sz;
    int *sbuf_sz;
    int *rreq_sz;
    int *sreq_sz;
    int *rreq;
    int *sreq;
};

/* Create the world collective */

static collective_t __world = {NULL, 0, 0};
collective_t *_world = &__world;
int _world_rank = 0;
int _world_size = 1;

/* collective checkpointer */
/* FIXME: SINCE RIGHT NOW, THERE IS ONLY THE WORLD COLLECTIVE AND NO WAY
   TO CREATE CHILDREN COLLECTIVES (NOT EVEN IN PRINCIPLE WITH THE CURRENT
   STATE OF RELAY), THIS IS BASICALLY A PLACEHOLDER. */

void checkpt_collective( const collective_t *comm )
{
    CHECKPT_VAL( int, world_rank );
    CHECKPT_VAL( int, world_size );
}

collective_t *restore_collective( void )
{
    int rank, size;
    RESTORE_VAL( int, rank );
    RESTORE_VAL( int, size );
    if ( size != world_size )
        ERROR( (
            "The number of processes that made this checkpt (%i) is different "
            "from the number of processes currently (%i)",
            size, world_size ) );
    if ( rank != world_rank )
        ERROR(
            ( "This process (%i) is reading a checkpoint previously written by "
              "a different process (%i)",
              rank, world_rank ) );
    return world;
}

/* mp checkpointer */

void checkpt_mp( mp_t *mp )
{
    int port;
    CHECKPT( mp, 1 );
    CHECKPT( mp->rbuf, mp->n_port );
    CHECKPT( mp->sbuf, mp->n_port );
    CHECKPT( mp->rbuf_sz, mp->n_port );
    CHECKPT( mp->sbuf_sz, mp->n_port );
    CHECKPT( mp->rreq_sz, mp->n_port );
    CHECKPT( mp->sreq_sz, mp->n_port );
    CHECKPT( mp->rreq, mp->n_port );
    CHECKPT( mp->sreq, mp->n_port );
    for ( port = 0; port < mp->n_port; port++ )
    {
        CHECKPT_ALIGNED( mp->rbuf[port], mp->rbuf_sz[port], 128 );
        CHECKPT_ALIGNED( mp->sbuf[port], mp->sbuf_sz[port], 128 );
    }
}

mp_t *restore_mp( void )
{
    mp_t *mp;
    int port;
    RESTORE( mp );
    RESTORE( mp->rbuf );
    RESTORE( mp->sbuf );
    RESTORE( mp->rbuf_sz );
    RESTORE( mp->sbuf_sz );
    RESTORE( mp->rreq_sz );
    RESTORE( mp->sreq_sz );
    RESTORE( mp->rreq );
    RESTORE( mp->sreq );
    for ( port = 0; port < mp->n_port; port++ )
    {
        RESTORE_ALIGNED( mp->rbuf[port] );
        RESTORE_ALIGNED( mp->sbuf[port] );
    }
    return mp;
}

struct RelayPolicy
{

    // FIXME: ERROR CODE HANDLING

    // FIXME-KJB: The whole sizing process in here is kinda silly and should
    // be removed in the long haul.

#define RESIZE_FACTOR 1.3125

    inline void boot_mp( int *pargc, char ***pargv )
    {
        ConnectionManager::instance().init( pargc, pargv );
        P2PConnection &p2p = P2PConnection::instance();
        __world.parent = NULL, __world.color = 0, __world.key = 0;
        _world_rank = p2p.global_id();
        _world_size = p2p.global_size();
        REGISTER_OBJECT( &__world, checkpt_collective, restore_collective,
                         NULL );
    }

    inline void halt_mp( void )
    {
        UNREGISTER_OBJECT( &__world );
        __world.parent = NULL, __world.color = 0, __world.key = 0;
        _world_size = 1;
        _world_rank = 0;
        P2PConnection::instance().post( P2PTag::end );
        ConnectionManager::instance().finalize();
    }

    inline void mp_abort( int reason )
    {
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::abort, P2PTag::data, 1, 0 );
        p2p.post( request );
        p2p.send( &reason, 1, P2PTag::data );
        p2p.abort( reason );
    }

    inline void mp_barrier( void )
    {
        P2PConnection &p2p = P2PConnection::instance();
        p2p.post( P2PTag::barrier );
        p2p.barrier();
    }

    inline void mp_allsum_d( double *local, double *global, int n )
    {
        if ( !local || !global || n < 1 || abs( local - global ) < n )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::allreduce_sum_double, P2PTag::data, n, 0 );
        p2p.post( request );
        p2p.send( local, request.count, request.tag );
        p2p.recv( global, request.count, request.tag, request.id );
    }

    inline void mp_allsum_i( int *local, int *global, int n )
    {
        if ( !local || !global || n < 1 || abs( local - global ) < n )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::allreduce_sum_int, P2PTag::data, n, 0 );
        p2p.post( request );
        p2p.send( local, request.count, request.tag );
        p2p.recv( global, request.count, request.tag, request.id );
    }

    inline void mp_allgather_i( int *sbuf, int *rbuf, int n )
    {
        if ( !sbuf || !rbuf || n < 1 )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::allgather_int, P2PTag::data, n, 0 );
        p2p.post( request );
        p2p.send( sbuf, request.count, request.tag );
        p2p.recv( rbuf, request.count * world_size, request.tag, request.id );
    }

    inline void mp_allgather_i64( int64_t *sbuf, int64_t *rbuf, int n )
    {
        if ( !sbuf || !rbuf || n < 1 )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::allgather_int64, P2PTag::data, n, 0 );
        p2p.post( request );
        p2p.send( sbuf, request.count, request.tag );
        p2p.recv( rbuf, request.count * world_size, request.tag, request.id );
    }

    inline void mp_gather_uc( unsigned char *sbuf, unsigned char *rbuf, int n )
    {
        if ( !sbuf || ( !rbuf && world_rank == 0 ) || n < 1 )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::gather_uc, P2PTag::data, n, 0 );
        p2p.post( request );
        p2p.send( sbuf, request.count, request.tag );
        if ( world_rank == 0 )
            p2p.recv( rbuf, request.count * world_size, request.tag,
                      request.id );
    }

    inline void mp_send_i( int *buf, int n, int dst )
    {
        if ( !buf || n < 1 || dst < 0 || dst >= world_size )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::send, P2PTag::data, n, 0, dst );
        p2p.post( request );
        p2p.send( buf, request.count, request.tag );
    }

    inline void mp_recv_i( int *buf, int n, int src )
    {
        if ( !buf || n < 1 || src < 0 || src >= world_size )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::recv, P2PTag::data, n, 0, src );
        p2p.post( request );
        p2p.recv( buf, request.count, request.tag, request.id );
    }

    /* ---- BEGIN EXACT CUT-AND-PASTE JOB FROM DMPPOLICY ---- */
    /* FIXME-KJB: AT THIS POINT, MUCH OF MP IN DMP AND RELAY COULD BE EXTRACTED
       INTO A UNIFIED IMPLEMENTATION (AND, AT THE SAME TIME, THE API FIXED) */

    inline mp_t *new_mp( int n_port )
    {
        mp_t *mp;
        if ( n_port < 1 )
            ERROR( ( "Bad args" ) );
        MALLOC( mp, 1 );
        mp->n_port = n_port;
        MALLOC( mp->rbuf, n_port );
        MALLOC( mp->sbuf, n_port );
        MALLOC( mp->rbuf_sz, n_port );
        MALLOC( mp->sbuf_sz, n_port );
        MALLOC( mp->rreq_sz, n_port );
        MALLOC( mp->sreq_sz, n_port );
        MALLOC( mp->rreq, n_port );
        MALLOC( mp->sreq, n_port );
        CLEAR( mp->rbuf, n_port );
        CLEAR( mp->sbuf, n_port );
        CLEAR( mp->rbuf_sz, n_port );
        CLEAR( mp->sbuf_sz, n_port );
        CLEAR( mp->rreq_sz, n_port );
        CLEAR( mp->sreq_sz, n_port );
        CLEAR( mp->rreq, n_port );
        CLEAR( mp->sreq, n_port );
        REGISTER_OBJECT( mp, checkpt_mp, restore_mp, NULL );
        return mp;
    }

    inline void delete_mp( mp_t *mp )
    {
        int port;
        if ( !mp )
            return;
        UNREGISTER_OBJECT( mp );
        for ( port = 0; port < mp->n_port; port++ )
        {
            FREE_ALIGNED( mp->rbuf[port] );
            FREE_ALIGNED( mp->sbuf[port] );
        }
        FREE( mp->rreq );
        FREE( mp->sreq );
        FREE( mp->rreq_sz );
        FREE( mp->sreq_sz );
        FREE( mp->rbuf_sz );
        FREE( mp->sbuf_sz );
        FREE( mp->rbuf );
        FREE( mp->sbuf );
        FREE( mp );
    }

    inline void *ALIGNED( 128 ) mp_recv_buffer( mp_t *mp, int port )
    {
        if ( !mp || port < 0 || port >= mp->n_port )
            ERROR( ( "Bad args" ) );
        return mp->rbuf[port];
    }

    inline void *ALIGNED( 128 ) mp_send_buffer( mp_t *mp, int port )
    {
        if ( !mp || port < 0 || port >= mp->n_port )
            ERROR( ( "Bad args" ) );
        return mp->sbuf[port];
    }

    inline void mp_size_recv_buffer( mp_t *mp, int port, int sz )
    {
        char *ALIGNED( 128 ) buf;

        if ( !mp || port < 0 || port > mp->n_port || sz < 1 )
            ERROR( ( "Bad args" ) );

        // If there already a large enough buffer, we are done
        if ( mp->rbuf_sz[port] >= sz )
            return;

        // Try to reduce the number of reallocs
        sz = (int)( sz * (double)RESIZE_FACTOR );

        // If no buffer allocated for this port, malloc it and return
        if ( !mp->rbuf[port] )
        {
            MALLOC_ALIGNED( mp->rbuf[port], sz, 128 );
            mp->rbuf_sz[port] = sz;
            return;
        }

        // Resize the existing buffer (preserving any data in it)
        // (FIXME: THIS IS PROBABLY SILLY!)
        MALLOC_ALIGNED( buf, sz, 128 );
        COPY( buf, mp->rbuf[port], mp->rbuf_sz[port] );
        FREE_ALIGNED( mp->rbuf[port] );
        mp->rbuf[port] = buf;
        mp->rbuf_sz[port] = sz;
    }

    inline void mp_size_send_buffer( mp_t *mp, int port, int sz )
    {
        char *ALIGNED( 128 ) buf;

        // Check input arguments
        if ( !mp || port < 0 || port > mp->n_port || sz < 1 )
            ERROR( ( "Bad args" ) );

        // Is there already a large enough buffer
        if ( mp->sbuf_sz[port] >= sz )
            return;

        // Try to reduce the number of reallocs
        sz = (int)( sz * (double)RESIZE_FACTOR );

        // If no buffer allocated for this port, malloc it and return
        if ( !mp->sbuf[port] )
        {
            MALLOC_ALIGNED( mp->sbuf[port], sz, 128 );
            mp->sbuf_sz[port] = sz;
            return;
        }

        // Resize the existing buffer (preserving any data in it)
        // (FIXME: THIS IS PROBABLY SILLY!)
        MALLOC_ALIGNED( buf, sz, 128 );
        COPY( buf, mp->sbuf[port], mp->sbuf_sz[port] );
        FREE_ALIGNED( mp->sbuf[port] );
        mp->sbuf[port] = buf;
        mp->sbuf_sz[port] = sz;
    }

    /* ---- END CUT-AND-PASTE JOB FROM DMPPOLICY ---- */

    inline void mp_begin_recv( mp_t *mp, int port, int sz, int src, int tag )
    {
        if ( !mp || port < 0 || port >= mp->n_port || sz < 1 ||
             sz > mp->rbuf_sz[port] || src < 0 || src >= world_size )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        mp->rreq_sz[port] = sz;
        MPRequest request( P2PTag::irecv, tag, sz, port, src );
        p2p.post( request );
        p2p.irecv( static_cast<char *>( mp->rbuf[port] ), sz, tag, port );
    }

    inline void mp_begin_send( mp_t *mp, int port, int sz, int dst, int tag )
    {
        if ( !mp || port < 0 || port >= mp->n_port || sz < 1 ||
             sz > mp->sbuf_sz[port] || dst < 0 || dst >= world_size )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        mp->sreq_sz[port] = sz;
        MPRequest request( P2PTag::isend, tag, sz, port, dst );
        p2p.post( request );
        p2p.isend( static_cast<char *>( mp->sbuf[port] ), sz, tag, port );
    }

    inline void mp_end_recv( mp_t *mp, int port )
    {
        int sz;
        if ( !mp || port < 0 || port >= mp->n_port )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::wait_recv, 0, 0, port );
        p2p.post( request );
        p2p.wait_recv( port );
        p2p.get_count<char>( port, sz );
        /* FIXME: SHOULDN'T WE CHECK SZ==RREQ_SZ HERE? */
    }

    inline void mp_end_send( mp_t *mp, int port )
    {
        if ( !mp || port < 0 || port >= mp->n_port )
            ERROR( ( "Bad args" ) );
        P2PConnection &p2p = P2PConnection::instance();
        MPRequest request( P2PTag::wait_send, 0, 0, port );
        p2p.post( request );
        p2p.wait_send( port );
    }

#undef RESIZE_FACTOR

}; // struct RelayPolicy

#endif // RelayPolicy_h
