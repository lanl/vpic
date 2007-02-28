/*
 * Written by:
 * Dr.-Ing. Benjamin Karl Bergen
 * Computational Physics and Methods (CCS-2)
 * Computer, Computational, and Statistical Sciences Division
 * Los Alamos National Laboratory
 * February 2007 - Original Version
 */

#if defined USE_MPI_STUBS

#include <mpi_stubs.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

int MPI_Recv( void *buf, int count, MPI_Datatype datatype,
   int source, int tag, MPI_Comm comm, MPI_Status *request ) {
   return( MPI_SUCCESS );
} /* MPI_Send; */

int MPI_Irecv( void *buf, int count, MPI_Datatype datatype,
   int source, int tag, MPI_Comm comm, MPI_Request *request ) {
   return( MPI_SUCCESS );
} /* MPI_Irecv */

int MPI_Send( void *buf, int count, MPI_Datatype datatype,
   int dest, int tag, MPI_Comm comm ) {
   return( MPI_SUCCESS );
} /* MPI_Send; */

int MPI_Issend( void *buf, int count, MPI_Datatype datatype,
   int dest, int tag, MPI_Comm comm, MPI_Request *request ) {
   return( MPI_SUCCESS );
} /* MPI_Issend */

int MPI_Wait( MPI_Request *request, MPI_Status *status ) {
   return( MPI_SUCCESS );
} /* MPI_Wait */

int MPI_Init( int *argc, char ***argv ) {
   return( MPI_SUCCESS );
} /* MPI_Init */

int MPI_Comm_rank( MPI_Comm comm, int *rank ) {
   *rank = 0;
   return( MPI_SUCCESS );
} /* MPI_Comm_rank */

int MPI_Comm_size( MPI_Comm comm, int *size ) {
   *size = 1;
   return( MPI_SUCCESS );
} /* MPI_Comm_size */

int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count) {
	*count = status->count;
   return( MPI_SUCCESS );
} /* MPI_Get_count */

double MPI_Wtime( void ) {
   static int sec = -1;
   struct timeval tv;
   gettimeofday(&tv, (void *)0);
   if (sec < 0) sec = tv.tv_sec;
   return (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
} /* MPI_Wtime */

int MPI_Barrier( MPI_Comm comm ) {
   return( MPI_SUCCESS );
} /* MPI_Barrier */

int  MPI_Finalize( void ) {
   return( MPI_SUCCESS );
} /* MPI_Finalize */

int MPI_Abort( MPI_Comm comm, int errorcode ) {
   exit( errorcode );
   return( MPI_SUCCESS );
} /* MPI_Abort */

int MPI_Allgather( void *sendbuf, int sendcount, MPI_Datatype sendtype,
   void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm ) {
   return( MPI_SUCCESS );
} /* MPI_Allgather */

int  MPI_Allreduce( void *sendbuf, void *recvbuf, int nitems,
   MPI_Datatype type, MPI_Op op, MPI_Comm comm ) {
   int i;
   if( type == MPI_INT ) {
      int *pd_sendbuf, *pd_recvbuf;
      pd_sendbuf = (int *) sendbuf;    
      pd_recvbuf = (int *) recvbuf;    

      for( i=0; i<nitems; i++ ) {
         *(pd_recvbuf+i) = *(pd_sendbuf+i);
      } /* for */
   } /* if */

   if( type == MPI_DOUBLE ) {
      double *pd_sendbuf, *pd_recvbuf;
      pd_sendbuf = (double *) sendbuf;    
      pd_recvbuf = (double *) recvbuf;    

      for( i=0; i<nitems; i++ ) {
         *(pd_recvbuf+i) = *(pd_sendbuf+i);
      } /* for */
   } /* if */

   return( MPI_SUCCESS );
} /* MPI_Allreduce */
  
#endif /* USE_MPI_STUBS */
