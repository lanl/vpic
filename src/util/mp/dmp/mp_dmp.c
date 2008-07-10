/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "mp_dmp.h"

// FIXME: AN MP SHOULD USE ITS OWN COMMUNICATOR _NOT_ MPI_COMM_WORLD!
// FIXME: WORLD_RANK AND WORLD_SIZE ... SHOULD BE GLOBAL CONST VARS

#if !defined(USE_MPRELAY)

// BJA: We're getting crashbugs from excessive realloc on some
// problems, so this multiplier will help reduce the number of such
// send/receive buffer resizes, though at a modest cost of excess
// memory
//
// KJB: The whole sizing process in here is kinda silly and should be
// removed in the long haul.  In the meantime, 1.3125 is better resize
// factor theoretically (the "silver" ratio).

#define RESIZE_FACTOR 1.3125

#define TRAP( x ) do {                                                  \
    int ierr = (x);                                                     \
    if( ierr!=MPI_SUCCESS ) ERROR(( "MPI error %i on "#x, ierr ));      \
  } while(0)

void
mp_init_dmp( int argc,
             char ** argv ) {
  TRAP( MPI_Init( &argc, &argv ) );
}

void
mp_finalize_dmp( mp_handle h ) {
  TRAP( MPI_Finalize() ); 
}

mp_handle
new_mp_dmp( void ) {
  mp_t *mp;
  MALLOC( mp, 1 );
  CLEAR( mp, 1 );
  TRAP( MPI_Comm_rank( MPI_COMM_WORLD, &mp->rank  ) );
  TRAP( MPI_Comm_size( MPI_COMM_WORLD, &mp->nproc ) );
  mp->elapsed_ref = MPI_Wtime();
  return (mp_handle)mp;
}

void
delete_mp_dmp( mp_handle *h ) {
  mp_t *mp;
  if( h==NULL ) return;
  mp = (mp_t *)*h;
  if( mp==NULL ) return;
  CLEAR( mp, 1 );
  FREE( mp );
  *h = NULL;
}

int
mp_rank_dmp( mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  return mp==NULL ? -1 : mp->rank;
}

int
mp_nproc_dmp( mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  return mp==NULL ? 0 : mp->nproc;
}

void * ALIGNED(16)
mp_recv_buffer_dmp( int tag,
                    mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  if( mp==NULL || tag<0 || tag>=NUM_BUF ) return NULL;
  return mp->rbuf[tag];
}

void * ALIGNED(16)
mp_send_buffer_dmp( int tag,
                    mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  if( mp==NULL || tag<0 || tag>=NUM_BUF ) return NULL;
  return mp->sbuf[tag];
}

void
mp_abort_dmp( int reason,
              mp_handle h ) {  
  MPI_Abort( MPI_COMM_WORLD, reason );
}

void
mp_barrier_dmp( mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  if( mp==NULL ) return;
  TRAP( MPI_Barrier( MPI_COMM_WORLD ) );
}

// Returns the time elapsed since the communicator was created.  Every
// rank gets the same value.
double
mp_elapsed_dmp( mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  double local_t, global_t;

  if( mp==NULL ) return -1;
  local_t = MPI_Wtime() - mp->elapsed_ref;
  TRAP( MPI_Allreduce( &local_t, &global_t, 1, MPI_DOUBLE,
                       MPI_MAX, MPI_COMM_WORLD ) );
  
  return global_t;
}

// Stop watch. Different ranks may get different values.  First call
// to stop watch returns a measure of the overhead.
double
mp_time00_dmp( mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL ) return -1;
  mp->time00_toggle ^= 1;
  if( mp->time00_toggle ) mp->time00_ref = MPI_Wtime();

  return MPI_Wtime() - mp->time00_ref;
}

double
mp_wtime_dmp(void) {
  return MPI_Wtime();
}

void
mp_size_recv_buffer_dmp( int tag,
                         int size,
                         mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  char * ALIGNED(16) rbuf;

  // Check input arguments
  if( mp==NULL              ) ERROR(( "Bad handle" ));
  if( tag<0 || tag>=NUM_BUF ) ERROR(( "Bad tag" ));
  if( size<=0               ) ERROR(( "Bad size" ));

  // If no buffer allocated for this tag
  if( mp->rbuf[tag]==NULL ) {
    MALLOC_ALIGNED( mp->rbuf[tag], size, 16 );
    mp->rbuf_size[tag] = size;
    return;
  }

  // Is there already a large enough buffer
  if( mp->rbuf_size[tag]>=size ) return;

  // Try to reduce the number of realloc calls

  size*=RESIZE_FACTOR; 

  // Create the new recv buffer

  MALLOC_ALIGNED( rbuf, size, 16 );

  // Preserve the old recv buffer data
  // FIXME: THIS IS PROBABLY SILLY!
  COPY( rbuf, mp->rbuf[tag], mp->rbuf_size[tag] );

  // Free the old recv buffer

  FREE_ALIGNED( mp->rbuf[tag] );
  mp->rbuf[tag]      = rbuf;
  mp->rbuf_size[tag] = size;
}

void
mp_size_send_buffer_dmp( int tag,
                         int size,
                         mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  char * ALIGNED(16) sbuf;

  // Check input arguments
  if( mp==NULL              ) ERROR(( "Bad handle" ));
  if( tag<0 || tag>=NUM_BUF ) ERROR(( "Bad tag" ));
  if( size<=0               ) ERROR(( "Bad size" ));

  // If no buffer allocated for this tag
  if( mp->sbuf[tag]==NULL ) {
    MALLOC_ALIGNED( mp->sbuf[tag], size, 16 );
    mp->sbuf_size[tag] = size;
    return;
  }

  // Is there already a large enough buffer
  if( mp->sbuf_size[tag]>=size ) return;

  // Try to reduce the number of realloc calls

  size *= RESIZE_FACTOR; 

  // Create the new send buffer

  MALLOC_ALIGNED( sbuf, size, 16 );

  // Preserve the old send buffer data
  // FIXME: SILLY
  COPY( sbuf, mp->sbuf[tag], mp->sbuf_size[tag] );

  // Free the old recv buffer

  FREE_ALIGNED( mp->sbuf[tag] );
  mp->sbuf[tag]      = sbuf;
  mp->sbuf_size[tag] = size;
}

void
mp_begin_recv_dmp( int recv_buf,
                   int msg_size,
                   int sender,
                   int msg_tag,
                   mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL                         ) ERROR(( "Bad handle" ));
  if( recv_buf<0 || recv_buf>=NUM_BUF  ) ERROR(( "Bad recv_buf" )); 
  if( msg_size<=0                      ) ERROR(( "Bad msg_size" ));
  if( sender<0 || sender>=mp->nproc    ) ERROR(( "Bad sender" ));
  if( mp->rbuf[recv_buf]==NULL         ) ERROR(( "NULL recv_buf" ));
  if( mp->rbuf_size[recv_buf]<msg_size )
    ERROR(( "recv_buf too small" ));

  mp->rreq_size[recv_buf] = msg_size;

  TRAP( MPI_Irecv( mp->rbuf[recv_buf], msg_size, MPI_BYTE,
                   sender, msg_tag, MPI_COMM_WORLD,
                   &mp->rreq[recv_buf] ) );
}

void
mp_begin_send_dmp( int send_buf,
                   int msg_size,
                   int receiver,
                   int msg_tag,
                   mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL                          ) ERROR(( "Bad handle" ));
  if( send_buf<0 || send_buf>=NUM_BUF   ) ERROR(( "Bad send_buf" ));
  if( msg_size<=0                       ) ERROR(( "Bad msg_size" ));
  if( receiver<0 || receiver>=mp->nproc ) ERROR(( "Bad receiver" ));
  if( mp->sbuf[send_buf]==NULL          ) ERROR(( "NULL send_buf" ));
  if( mp->sbuf_size[send_buf]<msg_size  )
    ERROR(( "send_buf too small" ));

  mp->sreq_size[send_buf] = msg_size;

  TRAP( MPI_Issend( mp->sbuf[send_buf], msg_size, MPI_BYTE,
                    receiver,	msg_tag, MPI_COMM_WORLD,
                    &mp->sreq[send_buf] ) );
}

void
mp_end_recv_dmp( int recv_buf,
                 mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  MPI_Status status;
  int size;

  if( mp==NULL                        ) ERROR(( "Bad handle" ));
  if( recv_buf<0 || recv_buf>=NUM_BUF ) ERROR(( "Bad recv_buf" ));

  TRAP( MPI_Wait( &mp->rreq[recv_buf], &status ) );

  TRAP( MPI_Get_count( &status, MPI_BYTE, &size ) );
  if( mp->rreq_size[recv_buf] != size )
    ERROR(( "Sizes do not match" ));
}

void
mp_end_send_dmp( int send_buf,
                 mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL                        ) ERROR(( "Bad handle" ));
  if( send_buf<0 || send_buf>=NUM_BUF ) ERROR(( "Bad send_buf" ));
  TRAP( MPI_Wait( &mp->sreq[send_buf], MPI_STATUS_IGNORE ) );
}

void
mp_allsum_d_dmp( double *local,
                 double *global,
                 int n,
                 mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  
  if( mp==NULL            ) ERROR(( "Bad handle" ));
  if( local==NULL         ) ERROR(( "Bad local" ));
  if( global==NULL        ) ERROR(( "Bad global" ));
  if( abs(local-global)<n ) ERROR(( "Overlapping local and global" ));

  TRAP( MPI_Allreduce( local, global, n, MPI_DOUBLE,
                       MPI_SUM, MPI_COMM_WORLD ) );
}

void
mp_allsum_i_dmp( int *local,
                 int *global,
                 int n,
                 mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  
  if( mp==NULL            ) ERROR(( "Bad handle" ));
  if( local==NULL         ) ERROR(( "Bad local" ));
  if( global==NULL        ) ERROR(( "Bad global" ));
  if( abs(local-global)<n ) ERROR(( "Overlapping local and global" ));
  
  TRAP( MPI_Allreduce( local, global, n, MPI_INT,
                       MPI_SUM, MPI_COMM_WORLD ) );
}

void
mp_allgather_i_dmp( int *sbuf,
                    int *rbuf,
                    int n,
                    mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  
  if( mp==NULL   ) ERROR(( "Bad handle" ));
  if( sbuf==NULL ) ERROR(( "Bad send" ));
  if( rbuf==NULL ) ERROR(( "Bad recv" ));
  if( n<1        ) ERROR(( "Bad n" ));
  
  TRAP( MPI_Allgather( sbuf, n, MPI_INT,
                       rbuf, n, MPI_INT, MPI_COMM_WORLD ) );
}

void
mp_allgather_i64_dmp( int64_t *sbuf,
                      int64_t *rbuf,
                      int n,
                      mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL   ) ERROR(( "Bad handle" ));
  if( sbuf==NULL ) ERROR(( "Bad send" ));
  if( rbuf==NULL ) ERROR(( "Bad recv" ));
  if( n<1        ) ERROR(( "Bad n" ));

  // FIXME: THIS IS BROKEN
  TRAP( MPI_Allgather( sbuf, n, MPI_LONG_LONG,
                       rbuf, n, MPI_LONG_LONG, MPI_COMM_WORLD ) );
}

// We need blocking send/receive to implement turnstiles.

void
mp_send_i_dmp( int *buf,
               int n,
               int dst,
               mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  
  if( mp==NULL ) ERROR(( "Bad handle" )); 
  TRAP( MPI_Send( buf, n, MPI_INT, dst, 0, MPI_COMM_WORLD ) );
}

void
mp_recv_i_dmp( int *buf,
               int n,
               int src,
               mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  MPI_Status status;
  
  if( mp==NULL ) ERROR(( "Bad handle" )); 
  TRAP( MPI_Recv( buf, n, MPI_INT, src, 0, MPI_COMM_WORLD, &status ) );
}

#endif // USE_MPRELAY
