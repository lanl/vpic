/* 
*/
#ifndef mp_t_h
#define mp_t_h

#include <util_base.h>
#include <mpi.h>

#define BEGIN_TURNSTYLE do {                                          \
   int _rank, _size, _baton = 0;                                      \
   MPI_Comm_rank( MPI_COMM_WORLD, &_rank );                           \
   MPI_Comm_size( MPI_COMM_WORLD, &_size );                           \
   if( _rank!=0 ) MPI_Recv( &_baton, 1, MPI_INT, _rank-1, _rank-1,    \
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE );      \
   do

#define END_TURNSTYLE                                                    \
     while(0);                                                           \
     if( _rank!=_size-1 ) MPI_Send( &_baton, 1, MPI_INT, _rank+1, _rank, \
                                    MPI_COMM_WORLD );                    \
   } while(0)

/* A mp can handle up to 54 simultaneous communications (27 sends and 27
	receives) ... This is based on every volume, face, edge and corner of a
	box simultaneously engaged in an overlapped send and an overlapped
	receive */
#define NUM_BUF 27

typedef struct mp {
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

#endif // mp_t_h
