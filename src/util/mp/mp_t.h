/* 
*/
#ifndef mp_t_h
#define mp_t_h

#include <mpi.h>
#include <common.h>

/* A mp can handle up to 54 simultaneous communications (27 sends and 27
	receives) ... This is based on every volume, face, edge and corner of a
	box simultaneously engaged in an overlapped send and an overlapped
	receive */
#define NUM_BUF 27

/* Opaque data type so users do not have to see message passing internals */
typedef void * mp_handle;

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
