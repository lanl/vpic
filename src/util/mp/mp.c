/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#if defined USE_MPI_STUBS
#include <mpi_stubs.h>
#else
#include <mpi.h>
#endif /* USE_MPI_STUBS */

#include <mp.h>
#include <string.h>  /* For memcpy */ 
/* A mp can handle up to 54 simultaneous communications (27 sends and 27
   receives) ... This is based on every volume, face, edge and corner of a
   box simultaneously engaged in an overlapped send and an overlapped
   receive */

#define NUM_BUF 27

/* BJA:
   We're getting crashbugs from excessive realloc on some problems, so 
   this multiplier will help reduce the number of such send/receive 
   buffer resizes, though at a modest cost of excess memory */ 

#define RESIZE_FACTOR 1.1

typedef struct _mp_t {
  int rank, nproc;
  double elapsed_ref;
  double time00_ref;
  int time00_toggle;
  char * ALIGNED rbuf[NUM_BUF];
  char * ALIGNED sbuf[NUM_BUF];
  int rbuf_size[NUM_BUF];
  int sbuf_size[NUM_BUF];
  MPI_Request rreq[NUM_BUF];
  MPI_Request sreq[NUM_BUF];
  int rreq_size[NUM_BUF];
  int sreq_size[NUM_BUF];
} mp_t;

mp_handle new_mp(void) {
  mp_t *mp;
  int i;

  mp = (mp_t *)malloc(sizeof(mp_t));
  if( mp==NULL ) return (mp_handle)NULL;

  MPI_Comm_rank( MPI_COMM_WORLD, &mp->rank );
  MPI_Comm_size( MPI_COMM_WORLD, &mp->nproc );
  mp->elapsed_ref = MPI_Wtime();
  mp->time00_ref = 0;
  mp->time00_toggle = 0;
  for( i=0; i<NUM_BUF; i++ ) {
    mp->rbuf[i] = NULL;
    mp->sbuf[i] = NULL;
    mp->rbuf_size[i] = 0;
    mp->sbuf_size[i] = 0;
    /* FIXME: Init rreq and sreq? */
    mp->rreq_size[i] = 0;
    mp->sreq_size[i] = 0;
  }

  return (mp_handle)mp;
}

void delete_mp( mp_handle *h ) {
  mp_t *mp;
  int i;

  if( h==NULL ) return;
  mp = (mp_t *)*h;
  if( mp==NULL ) return;
  mp->rank = -1;
  mp->nproc = 0;
  mp->elapsed_ref = 0;
  mp->time00_ref = 0;
  mp->time00_toggle = 0;
  for( i=0; i<NUM_BUF; i++ ) {
    if( mp->rbuf[i]!=NULL ) free_aligned(mp->rbuf[i]);
    if( mp->sbuf[i]!=NULL ) free_aligned(mp->sbuf[i]);
    mp->rbuf[i] = NULL;
    mp->sbuf[i] = NULL;
    mp->rbuf_size[i] = 0;
    mp->sbuf_size[i] = 0;
    /* FIXME: Delete rreq and sreq? */
    mp->rreq_size[i] = 0;
    mp->sreq_size[i] = 0;
  }
  free(mp);
  *h = NULL;
}

int mp_rank( mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  return mp==NULL ? -1 : mp->rank;
}

int mp_nproc( mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  return mp==NULL ? 0 : mp->nproc;
}

void * ALIGNED mp_recv_buffer( int tag, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  if( mp==NULL || tag<0 || tag>=NUM_BUF ) return NULL;
  return mp->rbuf[tag];
}

void * ALIGNED mp_send_buffer( int tag, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  if( mp==NULL || tag<0 || tag>=NUM_BUF ) return NULL;
  return mp->sbuf[tag];
}

void mp_abort( int reason, mp_handle h ) {  
  MPI_Abort(MPI_COMM_WORLD,reason);
}

/* Adding for clean termination */ 
void mp_finalize( mp_handle h ) {
  MPI_Finalize(); 
}


void mp_barrier( mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  if( mp==NULL ) return;
  MPI_Barrier(MPI_COMM_WORLD);
}

/* Returns the time elapsed since the communicator was created.
   Every rank gets the same value */
double mp_elapsed( mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  double local_t, global_t;

  if( mp==NULL ) return -1;
  local_t = MPI_Wtime() - mp->elapsed_ref;
  MPI_Allreduce( &local_t, &global_t, 1, MPI_DOUBLE,
		 MPI_MAX, MPI_COMM_WORLD );
  
  return global_t;
}

/* Stop watch. Different ranks may get different values.
   First call to stop watch returns a measure of the overhead */
double mp_time00( mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL ) return -1;
  mp->time00_toggle ^= 1;
  if( mp->time00_toggle ) mp->time00_ref = MPI_Wtime();

  return MPI_Wtime() - mp->time00_ref;
}

error_code mp_size_recv_buffer( int tag, int size, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  char * ALIGNED rbuf;

  /* Check input arguments */
  if( mp==NULL            ) return ERROR_CODE("Bad handle");
  if( tag<0 || tag>=NUM_BUF ) return ERROR_CODE("Bad tag");
  if( size<=0               ) return ERROR_CODE("Bad size");

  /* If no buffer allocated for this tag */
  if( mp->rbuf[tag]==NULL ) {
    mp->rbuf[tag] = (char * ALIGNED)
      malloc_aligned( size, preferred_alignment );
    if( mp->rbuf[tag]==NULL ) return ERROR_CODE("malloc_aligned failed");
    mp->rbuf_size[tag] = size;
    return SUCCESS;
  }

  /* Is there already a large enough buffer */
  if( mp->rbuf_size[tag]>=size ) return SUCCESS;

  /* Try to reduce the number of realloc calls */

  size*=RESIZE_FACTOR; 

  /* Create the new recv buffer */

  rbuf = (char * ALIGNED)malloc_aligned( size, preferred_alignment );
  if( rbuf==NULL ) return ERROR_CODE("malloc_aligned failed");

  /* Preserve the old recv buffer data */

  memcpy( rbuf, mp->rbuf[tag], mp->rbuf_size[tag] );

  /* Free the old recv buffer */

  free_aligned( mp->rbuf[tag] );
  mp->rbuf[tag]      = rbuf;
  mp->rbuf_size[tag] = size;

  return SUCCESS;
}

error_code mp_size_send_buffer( int tag, int size, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  char * ALIGNED sbuf;

  /* Check input arguments */
  if( mp==NULL              ) return ERROR_CODE("Bad handle");
  if( tag<0 || tag>=NUM_BUF ) return ERROR_CODE("Bad tag");
  if( size<=0               ) return ERROR_CODE("Bad size");

  /* If no buffer allocated for this tag */
  if( mp->sbuf[tag]==NULL ) {
    mp->sbuf[tag] = (char * ALIGNED)
      malloc_aligned( size, preferred_alignment );
    if( mp->sbuf[tag]==NULL ) return ERROR_CODE("malloc_aligned failed");
    mp->sbuf_size[tag] = size;
    return SUCCESS;
  }

  /* Is there already a large enough buffer */
  if( mp->sbuf_size[tag]>=size ) return SUCCESS;

  /* Try to reduce the number of realloc calls */

  size *= RESIZE_FACTOR; 

  /* Create the new send buffer */

  sbuf = (char * ALIGNED)malloc_aligned( size, preferred_alignment );
  if( sbuf==NULL ) return ERROR_CODE("malloc_aligned failed");

  /* Preserve the old send buffer data */

  memcpy( sbuf, mp->sbuf[tag], mp->sbuf_size[tag] );

  /* Free the old recv buffer */

  free_aligned( mp->sbuf[tag] );
  mp->sbuf[tag]      = sbuf;
  mp->sbuf_size[tag] = size;

  return SUCCESS;
}

error_code mp_begin_recv( int recv_buf, int msg_size,
                          int sender, int msg_tag, mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL                         ) return ERROR_CODE("Bad handle");
  if( recv_buf<0 || recv_buf>=NUM_BUF  ) return ERROR_CODE("Bad recv_buf"); 
  if( msg_size<=0                      ) return ERROR_CODE("Bad msg_size");
  if( sender<0 || sender>=mp->nproc    ) return ERROR_CODE("Bad sender");
  if( mp->rbuf[recv_buf]==NULL         ) return ERROR_CODE("NULL recv_buf");
  if( mp->rbuf_size[recv_buf]<msg_size )
    return ERROR_CODE("recv_buf too small");

  mp->rreq_size[recv_buf] = msg_size;
  switch( MPI_Irecv( mp->rbuf[recv_buf], msg_size, MPI_BYTE,
                     sender, msg_tag, MPI_COMM_WORLD,
                     &mp->rreq[recv_buf] ) ) {
  case MPI_SUCCESS:   return SUCCESS;
  case MPI_ERR_COMM:  return ERROR_CODE("MPI_ERR_COMM");
  case MPI_ERR_TYPE:  return ERROR_CODE("MPI_ERR_TYPE");
  case MPI_ERR_COUNT: return ERROR_CODE("MPI_ERR_COUNT");
  case MPI_ERR_TAG:   return ERROR_CODE("MPI_ERR_TAG");
  case MPI_ERR_RANK:  return ERROR_CODE("MPI_ERR_RANK");
  default:            break;
  }
  return ERROR_CODE("Unknown MPI error");
}

error_code mp_begin_send( int send_buf, int msg_size, int receiver,
                          int msg_tag, mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL                          ) return ERROR_CODE("Bad handle");
  if( send_buf<0 || send_buf>=NUM_BUF   ) return ERROR_CODE("Bad send_buf");
  if( msg_size<=0                       ) return ERROR_CODE("Bad msg_size");
  if( receiver<0 || receiver>=mp->nproc ) return ERROR_CODE("Bad receiver");
  if( mp->sbuf[send_buf]==NULL          ) return ERROR_CODE("NULL send_buf");
  if( mp->sbuf_size[send_buf]<msg_size  )
    return ERROR_CODE("send_buf too small");

  mp->sreq_size[send_buf] = msg_size;
  switch( MPI_Issend( mp->sbuf[send_buf], msg_size, MPI_BYTE,
                      receiver,	msg_tag, MPI_COMM_WORLD,
                      &mp->sreq[send_buf] ) ) {
  case MPI_SUCCESS:   return SUCCESS;
  case MPI_ERR_COMM:  return ERROR_CODE("MPI_ERR_COMM");
  case MPI_ERR_COUNT: return ERROR_CODE("MPI_ERR_COUNT");
  case MPI_ERR_TYPE:  return ERROR_CODE("MPI_ERR_TYPE");
  case MPI_ERR_TAG:   return ERROR_CODE("MPI_ERR_TAG");
  case MPI_ERR_RANK:  return ERROR_CODE("MPI_ERR_RANK");
  case MPI_ERR_OTHER: return ERROR_CODE("MPI_ERR_OTHER");
  default:            break;
  }
  return ERROR_CODE("Unknown MPI error");
}

error_code mp_end_recv( int recv_buf, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  MPI_Status status;
  int size;

  if( mp==NULL                        ) return ERROR_CODE("Bad handle");
  if( recv_buf<0 || recv_buf>=NUM_BUF ) return ERROR_CODE("Bad recv_buf");

  switch( MPI_Wait( &mp->rreq[recv_buf], &status ) ) {
  case MPI_SUCCESS:     break;
  case MPI_ERR_REQUEST: return ERROR_CODE("MPI_Wait - MPI_ERR_REQUEST");
  case MPI_ERR_ARG:     return ERROR_CODE("MPI_Wait - MPI_ERR_ARG");
  default:              return ERROR_CODE("MPI_Wait - Unknown MPI error");
  }

  switch( MPI_Get_count( &status, MPI_BYTE, &size ) ) {
  case MPI_SUCCESS:     break;
  case MPI_ERR_ARG:     return ERROR_CODE("MPI_Get_count - MPI_ERR_ARG");
  case MPI_ERR_TYPE:    return ERROR_CODE("MPI_Get_count - MPI_ERR_TYPE");
  default:              return ERROR_CODE("MPI_Get_count - Unknown MPI error");
  }
  if( mp->rreq_size[recv_buf] != size )
    return ERROR_CODE("Sizes do not match");

  return SUCCESS;
}

error_code mp_end_send( int send_buf, mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL                        ) return ERROR_CODE("Bad handle");
  if( send_buf<0 || send_buf>=NUM_BUF ) return ERROR_CODE("Bad send_buf");

  switch( MPI_Wait( &mp->sreq[send_buf], MPI_STATUS_IGNORE ) ) {
  case MPI_SUCCESS:     return SUCCESS;
  case MPI_ERR_REQUEST: return ERROR_CODE("MPI_ERR_REQUEST");
  case MPI_ERR_ARG:     return ERROR_CODE("MPI_ERR_ARG");
  default:              break;
  }
  return ERROR_CODE("Unknown MPI error");
}

error_code mp_allsum_d( double *local, double *global, int n, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  
  if( mp==NULL            ) return ERROR_CODE("Bad handle");
  if( local==NULL         ) return ERROR_CODE("Bad local");
  if( global==NULL        ) return ERROR_CODE("Bad global");
  if( abs(local-global)<n ) return ERROR_CODE("Overlapping local and global");
  
  switch( MPI_Allreduce( local, global, n, MPI_DOUBLE,
                         MPI_SUM, MPI_COMM_WORLD ) ) {
  case MPI_SUCCESS:    return SUCCESS;
  case MPI_ERR_COMM:   return ERROR_CODE("MPI_ERR_COMM");
  case MPI_ERR_BUFFER: return ERROR_CODE("MPI_ERR_BUFFER");
  case MPI_ERR_COUNT:  return ERROR_CODE("MPI_ERR_COUNT");
  case MPI_ERR_TYPE:   return ERROR_CODE("MPI_ERR_TYPE");
  case MPI_ERR_OP:     return ERROR_CODE("MPI_ERR_OP");
  default:             break;
  }
  return ERROR_CODE("Unknown MPI error");
}

error_code mp_allsum_i( int *local, int *global, int n, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  
  if( mp==NULL            ) return ERROR_CODE("Bad handle");
  if( local==NULL         ) return ERROR_CODE("Bad local");
  if( global==NULL        ) return ERROR_CODE("Bad global");
  if( abs(local-global)<n ) return ERROR_CODE("Overlapping local and global");
  
  switch( MPI_Allreduce( local, global, n, MPI_INT,
                         MPI_SUM, MPI_COMM_WORLD ) ) {
  case MPI_SUCCESS:    return SUCCESS;
  case MPI_ERR_COMM:   return ERROR_CODE("MPI_ERR_COMM");
  case MPI_ERR_BUFFER: return ERROR_CODE("MPI_ERR_BUFFER");
  case MPI_ERR_COUNT:  return ERROR_CODE("MPI_ERR_COUNT");
  case MPI_ERR_TYPE:   return ERROR_CODE("MPI_ERR_TYPE");
  case MPI_ERR_OP:     return ERROR_CODE("MPI_ERR_OP");
  default:             break;
  }
  return ERROR_CODE("Unknown MPI error");
}

error_code mp_allgather_i( int *sbuf, int *rbuf, int n, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  
  if( mp==NULL   ) return ERROR_CODE("Bad handle");
  if( sbuf==NULL ) return ERROR_CODE("Bad send");
  if( rbuf==NULL ) return ERROR_CODE("Bad recv");
  if( n<1        ) return ERROR_CODE("Bad n");
  
  switch( MPI_Allgather( sbuf, n, MPI_INT,
                         rbuf, n, MPI_INT, MPI_COMM_WORLD ) ) {
  case MPI_SUCCESS:    return SUCCESS;
  case MPI_ERR_COMM:   return ERROR_CODE("MPI_ERR_COMM");
  case MPI_ERR_OTHER:  return ERROR_CODE("MPI_ERR_OTHER");
  case MPI_ERR_COUNT:  return ERROR_CODE("MPI_ERR_COUNT");
  case MPI_ERR_TYPE:   return ERROR_CODE("MPI_ERR_TYPE");
  case MPI_ERR_BUFFER: return ERROR_CODE("MPI_ERR_BUFFER");
  default:             break;
  }
  return ERROR_CODE("Unknown MPI error");
}

error_code mp_allgather_i64( INT64_TYPE *sbuf, INT64_TYPE *rbuf, int n, mp_handle h ) {
  mp_t *mp = (mp_t *)h;

  if( mp==NULL ) return ERROR_CODE("Bad handle");
  if( sbuf==NULL ) return ERROR_CODE("Bad send");
  if( rbuf==NULL ) return ERROR_CODE("Bad recv");
  if( n<1      ) return ERROR_CODE("Bad n");

  switch( MPI_Allgather( sbuf, n, MPI_INTEGER8,
                         rbuf, n, MPI_INTEGER8, MPI_COMM_WORLD ) ) {
  case MPI_SUCCESS:  return SUCCESS;
  case MPI_ERR_COMM: return ERROR_CODE("MPI_ERR_COMM");
  case MPI_ERR_OTHER: return ERROR_CODE("MPI_ERR_OTHER");
  case MPI_ERR_COUNT: return ERROR_CODE("MPI_ERR_COUNT");
  case MPI_ERR_TYPE: return ERROR_CODE("MPI_ERR_TYPE");
  case MPI_ERR_BUFFER: return ERROR_CODE("MPI_ERR_BUFFER");
  default:           break;
  }
  return ERROR_CODE("Unknown MPI error");
}

/* We need blocking send/receive to implement turnstiles. */

error_code mp_send_i( int *buf, int n, int dst, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  
  if( mp==NULL   ) return ERROR_CODE("Bad handle");
  
  switch( MPI_Send( buf, n, MPI_INT, dst, 0, MPI_COMM_WORLD ) ) {

  case MPI_SUCCESS:    return SUCCESS;
  case MPI_ERR_COMM:   return ERROR_CODE("MPI_ERR_COMM");
  case MPI_ERR_OTHER:  return ERROR_CODE("MPI_ERR_OTHER");
  case MPI_ERR_COUNT:  return ERROR_CODE("MPI_ERR_COUNT");
  case MPI_ERR_TYPE:   return ERROR_CODE("MPI_ERR_TYPE");
  case MPI_ERR_BUFFER: return ERROR_CODE("MPI_ERR_BUFFER");
  default:             break;
  }
  return ERROR_CODE("Unknown MPI error");
}

error_code mp_recv_i( int *buf, int n, int src, mp_handle h ) {
  mp_t *mp = (mp_t *)h;
  MPI_Status status;
  
  if( mp==NULL   ) return ERROR_CODE("Bad handle");
  
  switch( MPI_Recv( buf, n, MPI_INT, src, 0, MPI_COMM_WORLD, &status ) ) {

  case MPI_SUCCESS:    return SUCCESS;
  case MPI_ERR_COMM:   return ERROR_CODE("MPI_ERR_COMM");
  case MPI_ERR_OTHER:  return ERROR_CODE("MPI_ERR_OTHER");
  case MPI_ERR_COUNT:  return ERROR_CODE("MPI_ERR_COUNT");
  case MPI_ERR_TYPE:   return ERROR_CODE("MPI_ERR_TYPE");
  case MPI_ERR_BUFFER: return ERROR_CODE("MPI_ERR_BUFFER");
  default:             break;
  }
  return ERROR_CODE("Unknown MPI error");
}


