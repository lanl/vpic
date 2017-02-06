#ifndef DMPPolicy_h
#define DMPPolicy_h

#include <mpi.h>
#include <cstdlib>
#include <cstdlib>

#include "../checkpt/checkpt.h"

/* Define this comm and mp opaque handles */
/* FIXME: PARENT, COLOR AND KEY ARE FOR FUTURE EXPANSION */

struct collective {
  collective_t * parent;
  int color;
  int key;
  MPI_Comm comm;
};

struct mp {
  int n_port;
  char * ALIGNED(128) * rbuf; char * ALIGNED(128) * sbuf;
  int * rbuf_sz;              int * sbuf_sz;
  int * rreq_sz;              int * sreq_sz;
  MPI_Request * rreq;         MPI_Request * sreq;
};

/* Create the world collective */

static collective_t __world = { NULL, 0, 0, MPI_COMM_SELF };
collective_t * _world = &__world;
int _world_rank = 0;
int _world_size = 1;

/* collective checkpointer */
/* FIXME: SINCE RIGHT NOW, THERE IS ONLY THE WORLD COLLECTIVE AND NO WAY
   TO CREATE CHILDREN COLLECTIVES, THIS IS BASICALLY A PLACEHOLDER. */

void
checkpt_collective( const collective_t * comm ) {
  CHECKPT_VAL( int, world_rank );
  CHECKPT_VAL( int, world_size );
}

collective_t *
restore_collective( void ) {
  int rank, size;
  RESTORE_VAL( int, rank );
  RESTORE_VAL( int, size );
  if( size!=world_size )
    ERROR(( "The number of nodes that made this checkpt (%i) is different "
            "from the number of nodes currently (%i)",
            size, world_size ));
  if( rank!=world_rank )
    ERROR(( "This node (%i) is reading a checkpoint previously written by "
            "a different node (%i).", rank, world_rank ));
  return world;
}

/* mp checkpointer */

void
checkpt_mp( mp_t * mp ) {
  int port;
  CHECKPT( mp, 1 );
  CHECKPT( mp->rbuf,    mp->n_port ); CHECKPT( mp->sbuf,    mp->n_port );
  CHECKPT( mp->rbuf_sz, mp->n_port ); CHECKPT( mp->sbuf_sz, mp->n_port );
  CHECKPT( mp->rreq_sz, mp->n_port ); CHECKPT( mp->sreq_sz, mp->n_port );
  CHECKPT( mp->rreq,    mp->n_port ); CHECKPT( mp->sreq,    mp->n_port );
  for( port=0; port<mp->n_port; port++ ) {
    CHECKPT_ALIGNED( mp->rbuf[port], mp->rbuf_sz[port], 128 );
    CHECKPT_ALIGNED( mp->sbuf[port], mp->sbuf_sz[port], 128 );
  }
}

mp_t *
restore_mp( void ) {
  mp_t * mp;
  int port;
  RESTORE( mp );
  RESTORE( mp->rbuf    ); RESTORE( mp->sbuf    );
  RESTORE( mp->rbuf_sz ); RESTORE( mp->sbuf_sz );
  RESTORE( mp->rreq_sz ); RESTORE( mp->sreq_sz );
  RESTORE( mp->rreq    ); RESTORE( mp->sreq    );
  for( port=0; port<mp->n_port; port++ ) {
    RESTORE_ALIGNED( mp->rbuf[port] );
    RESTORE_ALIGNED( mp->sbuf[port] );
  }
  return mp;
}

struct DMPPolicy {

  // FIXME-KJB: The whole sizing process in here is kinda silly and should
  // be removed in the long haul.

# define RESIZE_FACTOR 1.3125
# define TRAP( x ) do {                                                  \
     int ierr = (x);                                                     \
     if( ierr!=MPI_SUCCESS ) ERROR(( "MPI error %i on "#x, ierr ));      \
   } while(0)

  inline void
  boot_mp( int * pargc,
           char *** pargv ) {
    TRAP( MPI_Init( pargc, pargv ) );
    TRAP( MPI_Comm_dup( MPI_COMM_WORLD, &__world.comm ) );
    __world.parent = NULL, __world.color = 0, __world.key = 0;
    TRAP( MPI_Comm_rank( __world.comm, &_world_rank ) );
    TRAP( MPI_Comm_size( __world.comm, &_world_size ) );
    REGISTER_OBJECT( &__world, checkpt_collective, restore_collective, NULL );
  }

  inline void
  halt_mp( void ) {
    UNREGISTER_OBJECT( &__world );
    TRAP( MPI_Comm_free( &__world.comm ) );
    __world.parent = NULL, __world.color = 0, __world.key = 0;
    __world.comm = MPI_COMM_SELF;
    _world_size = 1;
    _world_rank = 0;
    TRAP( MPI_Finalize() );
  }

  inline void
  mp_abort( int reason ) {
    MPI_Abort( world->comm, reason );
  }

  inline void
  mp_barrier( void ) {
    TRAP( MPI_Barrier( world->comm ) );
  }

  inline void
  mp_allsum_d( double * local,
               double * global,
               int n ) {
    if( !local || !global || n<1 || std::abs(local-global)<n ) {
	 	ERROR(( "Bad args" ));
	 } // if
    TRAP( MPI_Allreduce( local, global, n, MPI_DOUBLE, MPI_SUM, world->comm ) );
  }

  inline void
  mp_allsum_i( int * local,
               int * global,
               int n ) {
    if( !local || !global || n<1 || std::abs(local-global)<n ) {
	 	ERROR(( "Bad args" ));
	 } // if
    TRAP( MPI_Allreduce( local, global, n, MPI_INT, MPI_SUM, world->comm ) );
  }

  inline void
  mp_allgather_i( int * sbuf,
                  int * rbuf,
                  int n ) {
    if( !sbuf || !rbuf || n<1 ) ERROR(( "Bad args" ));
    TRAP( MPI_Allgather( sbuf, n, MPI_INT, rbuf, n, MPI_INT, world->comm ) );
  }

  inline void
  mp_allgather_i64( int64_t * sbuf,
                    int64_t * rbuf,
                    int n ) {
    if( !sbuf || !rbuf || n<1 ) ERROR(( "Bad args" ));
    TRAP( MPI_Allgather( sbuf, n, MPI_LONG_LONG, rbuf, n, MPI_LONG_LONG, world->comm ) );
  }

  inline void
  mp_gather_uc( unsigned char * sbuf,
                unsigned char * rbuf,
                int n ) {
    if( !sbuf || (!rbuf && world_rank==0) || n<1 ) ERROR(( "Bad args" ));
    TRAP( MPI_Gather( sbuf, n, MPI_CHAR, rbuf, n, MPI_CHAR, 0, world->comm ) );
  }

  inline void
  mp_send_i( int * buf,
             int n,
             int dst ) {
    if( !buf || n<1 || dst<0 || dst>=world_size ) ERROR(( "Bad args" ));
    TRAP( MPI_Send( buf, n, MPI_INT, dst, 0, world->comm ) );
  }

  inline void
  mp_recv_i( int * buf,
             int n,
             int src ) {
    if( !buf || n<1 || src<0 || src>=world_size ) ERROR(( "Bad args" ));
    TRAP( MPI_Recv( buf, n, MPI_INT, src, 0, world->comm, MPI_STATUS_IGNORE ) );
  }

  inline mp_t *
  new_mp( int n_port ) {
    mp_t * mp;
    if( n_port<1 ) ERROR(( "Bad args" ));
    MALLOC( mp, 1 );
    mp->n_port = n_port;
    MALLOC( mp->rbuf,    n_port ); MALLOC( mp->sbuf,    n_port );
    MALLOC( mp->rbuf_sz, n_port ); MALLOC( mp->sbuf_sz, n_port );
    MALLOC( mp->rreq_sz, n_port ); MALLOC( mp->sreq_sz, n_port );
    MALLOC( mp->rreq,    n_port ); MALLOC( mp->sreq,    n_port );
    CLEAR(  mp->rbuf,    n_port ); CLEAR(  mp->sbuf,    n_port );
    CLEAR(  mp->rbuf_sz, n_port ); CLEAR(  mp->sbuf_sz, n_port );
    CLEAR(  mp->rreq_sz, n_port ); CLEAR(  mp->sreq_sz, n_port );
    CLEAR(  mp->rreq,    n_port ); CLEAR(  mp->sreq,    n_port );
    REGISTER_OBJECT( mp, checkpt_mp, restore_mp, NULL );
    return mp;
  }

  inline void
  delete_mp( mp_t * mp ) {
    int port;
    if( !mp ) return;
    UNREGISTER_OBJECT( mp );
    for( port=0; port<mp->n_port; port++ ) {
      FREE_ALIGNED( mp->rbuf[port] ); FREE_ALIGNED( mp->sbuf[port] );
    }
    FREE( mp->rreq    ); FREE( mp->sreq    );
    FREE( mp->rreq_sz ); FREE( mp->sreq_sz );
    FREE( mp->rbuf_sz ); FREE( mp->sbuf_sz );
    FREE( mp->rbuf    ); FREE( mp->sbuf    );
    FREE( mp );
  }

  inline void * ALIGNED(128)
  mp_recv_buffer( mp_t * mp,
                  int port ) {
    if( !mp || port<0 || port>=mp->n_port ) ERROR(( "Bad args" ));
    return mp->rbuf[port];
  }

  inline void * ALIGNED(128)
  mp_send_buffer( mp_t * mp,
                  int port ) {
    if( !mp || port<0 || port>=mp->n_port ) ERROR(( "Bad args" ));
    return mp->sbuf[port];
  }

  inline void
  mp_size_recv_buffer( mp_t * mp,
                       int port,
                       int sz ) {
    char * ALIGNED(128) buf;

    if( !mp || port<0 || port>mp->n_port || sz<1 ) ERROR(( "Bad args" ));

    // If there already a large enough buffer, we are done
    if( mp->rbuf_sz[port]>=sz ) return;

    // Try to reduce the number of reallocs
    sz = (int)( sz*(double)RESIZE_FACTOR );

    // If no buffer allocated for this port, malloc it and return
    if( !mp->rbuf[port] ) {
      MALLOC_ALIGNED( mp->rbuf[port], sz, 128 );
      mp->rbuf_sz[port] = sz;
      return;
    }

    // Resize the existing buffer (preserving any data in it)
    // (FIXME: THIS IS PROBABLY SILLY!)
    MALLOC_ALIGNED( buf, sz, 128 );
    COPY( buf, mp->rbuf[port], mp->rbuf_sz[port] );
    FREE_ALIGNED( mp->rbuf[port] );
    mp->rbuf[port]    = buf;
    mp->rbuf_sz[port] = sz;
  }

  inline void
  mp_size_send_buffer( mp_t * mp,
                       int port,
                       int sz ) {
    char * ALIGNED(128) buf;

    // Check input arguments
    if( !mp || port<0 || port>mp->n_port || sz<1 ) ERROR(( "Bad args" ));

    // Is there already a large enough buffer
    if( mp->sbuf_sz[port]>=sz ) return;

    // Try to reduce the number of reallocs
    sz = (int)( sz*(double)RESIZE_FACTOR );

    // If no buffer allocated for this port, malloc it and return
    if( !mp->sbuf[port] ) {
      MALLOC_ALIGNED( mp->sbuf[port], sz, 128 );
      mp->sbuf_sz[port] = sz;
      return;
    }

    // Resize the existing buffer (preserving any data in it)
    // (FIXME: THIS IS PROBABLY SILLY!)
    MALLOC_ALIGNED( buf, sz, 128 );
    COPY( buf, mp->sbuf[port], mp->sbuf_sz[port] );
    FREE_ALIGNED( mp->sbuf[port] );
    mp->sbuf[port]    = buf;
    mp->sbuf_sz[port] = sz;
  }

  inline void
  mp_begin_recv( mp_t * mp,
                 int port,
                 int sz,
                 int src,
                 int tag ) {
    if( !mp || port<0 || port>=mp->n_port || sz<1 || sz>mp->rbuf_sz[port] ||
        src<0 || src>=world_size ) ERROR(( "Bad args" ));
    mp->rreq_sz[port] = sz;
    TRAP(MPI_Irecv(mp->rbuf[port], sz, MPI_BYTE, src, tag, world->comm, &mp->rreq[port]));
  }

  inline void
  mp_begin_send( mp_t * mp,
                 int port,
                 int sz,
                 int dst,
                 int tag ) {
    if( !mp || port<0 || port>=mp->n_port || dst<0 || dst>=world_size ||
        sz<1 || mp->sbuf_sz[port]<sz ) ERROR(( "Bad args" ));
    mp->sreq_sz[port] = sz;
    TRAP(MPI_Issend(mp->sbuf[port],sz, MPI_BYTE, dst, tag, world->comm, &mp->sreq[port]));
  }

  inline void
  mp_end_recv( mp_t * mp,
               int port ) {
    MPI_Status status;
    int sz;
    if( !mp || port<0 || port>=mp->n_port ) ERROR(( "Bad args" ));
    TRAP( MPI_Wait( &mp->rreq[port], &status ) );
    TRAP( MPI_Get_count( &status, MPI_BYTE, &sz ) );
    if( mp->rreq_sz[port]!=sz ) ERROR(( "Sizes do not match" ));
  }

  inline void
  mp_end_send( mp_t * mp,
               int port ) {
    if( !mp || port<0 || port>=mp->n_port ) ERROR(( "Bad args" ));
    TRAP( MPI_Wait( &mp->sreq[port], MPI_STATUS_IGNORE ) );
  }

# undef RESIZE_FACTOR
# undef TRAP

}; // struct DMPPolicy


#endif // DMPPolicy_h
