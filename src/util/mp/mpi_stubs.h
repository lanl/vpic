/*
 * Written by:
 * Dr.-Ing. Benjamin Karl Bergen
 * Computational Physics and Methods (CCS-2)
 * Computer, Computational, and Statistical Sciences Division
 * Los Alamos National Laboratory
 * February 2007 - Original Version
 */

#ifndef mpi_dummy_h
#define mpi_dummy_h

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

/* constants */

#define MPI_COMM_WORLD 0

#define MPI_MAX 1
#define MPI_SUM 2
#define MPI_MIN 3

#define MPI_DATATYPE_NULL 0
#define MPI_CHAR 1
#define MPI_SHORT 2
#define MPI_INT 3
#define MPI_LONG 4
#define MPI_UNSIGNED_CHAR 5
#define MPI_UNSIGNED_SHORT 6
#define MPI_UNSIGNED 7
#define MPI_UNSIGNED_LONG 8
#define MPI_FLOAT 9
#define MPI_DOUBLE 10
#define MPI_LONG_DOUBLE 11
#define MPI_LONG_LONG 12
#define MPI_INTEGER 13
#define MPI_REAL 14
#define MPI_DOUBLE_PRECISION 15
#define MPI_COMPLEX 16
#define MPI_DOUBLE_COMPLEX 17
#define MPI_LOGICAL 18
#define MPI_CHARACTER 19
#define MPI_INTEGER1 20
#define MPI_INTEGER2 21
#define MPI_INTEGER4 22
#define MPI_INTEGER8 23
#define MPI_REAL4 24
#define MPI_REAL8 25
#define MPI_REAL16 26
#define MPI_BYTE 27
#define MPI_PACKED 28
#define MPI_UB 29
#define MPI_LB 30
#define MPI_FLOAT_INT 31
#define MPI_DOUBLE_INT 32
#define MPI_LONG_INT 33
#define MPI_2INT 34
#define MPI_SHORT_INT 35
#define MPI_LONG_DOUBLE_INT 36
#define MPI_2REAL 37
#define MPI_2DOUBLE_PRECISION 38
#define MPI_2INTEGER 39
#define MPI_LONG_LONG_INT 40

#define MPI_SUCCESS                   0
#define MPI_ERR_BUFFER                1
#define MPI_ERR_COUNT                 2
#define MPI_ERR_TYPE                  3
#define MPI_ERR_TAG                   4
#define MPI_ERR_COMM                  5
#define MPI_ERR_RANK                  6
#define MPI_ERR_REQUEST               7
#define MPI_ERR_ROOT                  8
#define MPI_ERR_GROUP                 9
#define MPI_ERR_OP                    10
#define MPI_ERR_TOPOLOGY              11
#define MPI_ERR_DIMS                  12
#define MPI_ERR_ARG                   13
#define MPI_ERR_UNKNOWN               14
#define MPI_ERR_TRUNCATE              15
#define MPI_ERR_OTHER                 16
#define MPI_ERR_INTERN                17
#define MPI_ERR_IN_STATUS             18
#define MPI_ERR_PENDING               19
#define MPI_ERR_ACCESS                20
#define MPI_ERR_AMODE                 21
#define MPI_ERR_ASSERT                22
#define MPI_ERR_BAD_FILE              23
#define MPI_ERR_BASE                  24
#define MPI_ERR_CONVERSION            25
#define MPI_ERR_DISP                  26
#define MPI_ERR_DUP_DATAREP           27
#define MPI_ERR_FILE_EXISTS           28
#define MPI_ERR_FILE_IN_USE           29
#define MPI_ERR_FILE                  30
#define MPI_ERR_INFO_KEY              31
#define MPI_ERR_INFO_NOKEY            32
#define MPI_ERR_INFO_VALUE            33
#define MPI_ERR_INFO                  34
#define MPI_ERR_IO                    35
#define MPI_ERR_KEYVAL                36
#define MPI_ERR_LOCKTYPE              37
#define MPI_ERR_NAME                  38
#define MPI_ERR_NO_MEM                39
#define MPI_ERR_NOT_SAME              40
#define MPI_ERR_NO_SPACE              41
#define MPI_ERR_NO_SUCH_FILE          42
#define MPI_ERR_PORT                  43
#define MPI_ERR_QUOTA                 44
#define MPI_ERR_READ_ONLY             45
#define MPI_ERR_RMA_CONFLICT          46
#define MPI_ERR_RMA_SYNC              47
#define MPI_ERR_SERVICE               48
#define MPI_ERR_SIZE                  49
#define MPI_ERR_SPAWN                 50
#define MPI_ERR_UNSUPPORTED_DATAREP   51
#define MPI_ERR_UNSUPPORTED_OPERATION 52
#define MPI_ERR_WIN                   53
#define MPI_ERR_LASTCODE              54

#define MPI_ERR_SYSRESOURCE          -2

/* data structures */

typedef struct {
   int count;
   int MPI_SOURCE;
   int MPI_TAG;
   int MPI_ERROR;
} MPI_Status;

#define MPI_STATUS_IGNORE ((MPI_Status *) 0)

typedef int MPI_Request;

typedef int MPI_Datatype;

typedef int MPI_Comm;

typedef int MPI_Op;

/* prototypes */

int  MPI_Init( int *argc, char ***argv );

int MPI_Finalize( void );

int MPI_Comm_rank( MPI_Comm comm, int *rank );

int MPI_Comm_size( MPI_Comm comm, int *size );

int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);

double MPI_Wtime();

int MPI_Send( void *buf, int count, MPI_Datatype datatype, int dest,
   int tag, MPI_Comm comm );

int MPI_Issend( void *buf, int count, MPI_Datatype datatype, int dest,
   int tag, MPI_Comm comm, MPI_Request *request );

int MPI_Recv( void *buf, int count, MPI_Datatype datatype,
   int source, int tag, MPI_Comm comm, MPI_Status *request );

int MPI_Irecv( void *buf, int count, MPI_Datatype datatype, int source,
   int tag, MPI_Comm comm, MPI_Request *request );

int MPI_Allgather( void *sendbuf, int sendcount, MPI_Datatype sendtype,
   void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm );

int  MPI_Allreduce( void *sendbuf, void *recvbuf, int nitems,
   MPI_Datatype type, MPI_Op op, MPI_Comm comm );

int MPI_Barrier( MPI_Comm comm );

int MPI_Wait(MPI_Request *request, MPI_Status *status);

int MPI_Abort( MPI_Comm comm, int errorcode );

#if defined(__cplusplus)
} /* extern */
#endif /* __cplusplus */

#endif /* mpi_dummy_h */
