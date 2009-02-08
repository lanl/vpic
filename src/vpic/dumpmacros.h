/*------------------------------------------------------------------------------
 * MACROS for I/O
 *
 * From origianl dump.cxx by Kevin Bowers
 * Created 10/2008
------------------------------------------------------------------------------*/
#ifndef dumpmacros_h
#define dumpmacros_h

#define WRITE_HEADER_V0(dump_type,sp_id,q_m,fileIO) BEGIN_PRIMITIVE {      \
  /* Binary compatibility information */                                   \
  WRITE( char,      CHAR_BIT,               fileIO );                      \
  WRITE( char,      sizeof(short int),      fileIO );                      \
  WRITE( char,      sizeof(int),            fileIO );                      \
  WRITE( char,      sizeof(float),          fileIO );                      \
  WRITE( char,      sizeof(double),         fileIO );                      \
  WRITE( short int, 0xcafe,                 fileIO );                      \
  WRITE( int,       0xdeadbeef,             fileIO );                      \
  WRITE( float,     1.0,                    fileIO );                      \
  WRITE( double,    1.0,                    fileIO );                      \
  /* Dump type and header format version */                                \
  WRITE( int,       0 /* Version */,        fileIO );                      \
  WRITE( int,       dump_type,              fileIO );                      \
  /* High level information */                                             \
  WRITE( int,       step,                   fileIO );                      \
  WRITE( int,       nxout,                  fileIO );                      \
  WRITE( int,       nyout,                  fileIO );                      \
  WRITE( int,       nzout,                  fileIO );                      \
  WRITE( float,     grid->dt,               fileIO );                      \
  WRITE( float,     dxout,                  fileIO );                      \
  WRITE( float,     dyout,                  fileIO );                      \
  WRITE( float,     dzout,                  fileIO );                      \
  WRITE( float,     grid->x0,               fileIO );                      \
  WRITE( float,     grid->y0,               fileIO );                      \
  WRITE( float,     grid->z0,               fileIO );                      \
  WRITE( float,     grid->cvac,             fileIO );                      \
  WRITE( float,     grid->eps0,             fileIO );                      \
  WRITE( float,     grid->damp,             fileIO );                      \
  WRITE( int,       mp_rank(grid->mp),      fileIO );                      \
  WRITE( int,       mp_nproc(grid->mp),     fileIO );                      \
  /* Species parameters */                                                 \
  WRITE( int,       sp_id,                  fileIO );                      \
  WRITE( float,     q_m,                    fileIO );                      \
} END_PRIMITIVE
 
// Note dim _MUST_ be a pointer to an int
 
#define WRITE_ARRAY_HEADER(p,ndim,dim,fileIO) BEGIN_PRIMITIVE { \
  WRITE( int, sizeof(p[0]), fileIO );                           \
  WRITE( int, ndim,         fileIO );                           \
  fileIO.write( dim, ndim );                                    \
} END_PRIMITIVE
 
// The WRITE macro copies the output "value" into a temporary variable
// of the requested output "type" so that the write to the "file"
// occurs from a known binary data type. For example, if grid.dx were
// changed from a float to a double, routines which say
// WRITE(float,grid.dx,out) will still work fine without modification
// (you will lose some precision in the output file obviously). Note:
// No effort is made to convert the raw binary data from the root
// node's native format. Usually, integer data will be little endian
// 32-bit values and floating data with be little endian 32-bit IEEE
// single precision write copies. However, specialty types could be
// created so that the type cast __WRITE_tmp = (type)(value)
// automatically does the underlying conversion in C++
 
#define WRITE(type,value,fileIO) BEGIN_PRIMITIVE { \
  type __WRITE_tmp = (type)(value);                \
  fileIO.write( &__WRITE_tmp, 1 );                 \
} END_PRIMITIVE
 
// Note: strlen does not include the terminating NULL
#define WRITE_STRING(string,fileIO) BEGIN_PRIMITIVE {     \
  int __WRITE_STRING_len = 0;                             \
  if( string!=NULL ) __WRITE_STRING_len = strlen(string); \
  fileIO.write( &__WRITE_STRING_len, 1 );                 \
  if( __WRITE_STRING_len>0 )                              \
    fileIO.write( string, __WRITE_STRING_len );           \
} END_PRIMITIVE
 
#define READ(type,value,fileIO) BEGIN_PRIMITIVE { \
  type __READ_tmp;                                \
  fileIO.read(&__READ_tmp, 1 );                   \
  (value) = __READ_tmp;                           \
} END_PRIMITIVE

#define F_WRITE_HEADER_V0(dump_type,sp_id,q_m,fileIO) BEGIN_PRIMITIVE { \
  /* Binary compatibility information */                    \
  F_WRITE( char,      CHAR_BIT,               fileIO );     \
  F_WRITE( char,      sizeof(short int),      fileIO );     \
  F_WRITE( char,      sizeof(int),            fileIO );     \
  F_WRITE( char,      sizeof(float),          fileIO );     \
  F_WRITE( char,      sizeof(double),         fileIO );     \
  F_WRITE( short int, 0xcafe,                 fileIO );     \
  F_WRITE( int,       0xdeadbeef,             fileIO );     \
  F_WRITE( float,     1.0,                    fileIO );     \
  F_WRITE( double,    1.0,                    fileIO );     \
  /* Dump type and header format version */                 \
  F_WRITE( int,       0 /* Version */,        fileIO );     \
  F_WRITE( int,       dump_type,              fileIO );     \
  /* High level information */                              \
  F_WRITE( int,       step,                   fileIO );     \
  F_WRITE( int,       imxstr-2,               fileIO );     \
  F_WRITE( int,       jmxstr-2,               fileIO );     \
  F_WRITE( int,       kmxstr-2,               fileIO );     \
  F_WRITE( float,     grid->dt,               fileIO );     \
  F_WRITE( float,     dxstr,                  fileIO );     \
  F_WRITE( float,     dystr,                  fileIO );     \
  F_WRITE( float,     dzstr,                  fileIO );     \
  F_WRITE( float,     grid->x0,               fileIO );     \
  F_WRITE( float,     grid->y0,               fileIO );     \
  F_WRITE( float,     grid->z0,               fileIO );     \
  F_WRITE( float,     grid->cvac,             fileIO );     \
  F_WRITE( float,     grid->eps0,             fileIO );     \
  F_WRITE( float,     grid->damp,             fileIO );     \
  F_WRITE( int,       mp_rank(grid->mp),      fileIO );     \
  F_WRITE( int,       mp_nproc(grid->mp),     fileIO );     \
  /* Species parameters */                                  \
  F_WRITE( int,       sp_id,                  fileIO );     \
  F_WRITE( float,     q_m,                    fileIO );     \
} END_PRIMITIVE
 
#define F_WRITE_HEADER_PAR(dump_type,sp_id,q_m,fileIO) BEGIN_PRIMITIVE { \
  /* Binary compatibility information */                    \
  F_WRITE( char,      CHAR_BIT,               fileIO );     \
  F_WRITE( char,      sizeof(short int),      fileIO );     \
  F_WRITE( char,      sizeof(int),            fileIO );     \
  F_WRITE( char,      sizeof(float),          fileIO );     \
  F_WRITE( char,      sizeof(double),         fileIO );     \
  F_WRITE( short int, 0xcafe,                 fileIO );     \
  F_WRITE( int,       0xdeadbeef,             fileIO );     \
  F_WRITE( float,     1.0,                    fileIO );     \
  F_WRITE( double,    1.0,                    fileIO );     \
  /* Dump type and header format version */                 \
  F_WRITE( int,       0 /* Version */,        fileIO );     \
  F_WRITE( int,       dump_type,              fileIO );     \
  /* High level information */                              \
  F_WRITE( int,       step,                   fileIO );     \
  F_WRITE( int,       grid->nx,               fileIO );     \
  F_WRITE( int,       grid->ny,               fileIO );     \
  F_WRITE( int,       grid->nz,               fileIO );     \
  F_WRITE( float,     grid->dt,               fileIO );     \
  F_WRITE( float,     grid->dx,               fileIO );     \
  F_WRITE( float,     grid->dy,               fileIO );     \
  F_WRITE( float,     grid->dz,               fileIO );     \
  F_WRITE( float,     grid->x0,               fileIO );     \
  F_WRITE( float,     grid->y0,               fileIO );     \
  F_WRITE( float,     grid->z0,               fileIO );     \
  F_WRITE( float,     grid->cvac,             fileIO );     \
  F_WRITE( float,     grid->eps0,             fileIO );     \
  F_WRITE( float,     grid->damp,             fileIO );     \
  F_WRITE( int,       mp_rank(grid->mp),      fileIO );     \
  F_WRITE( int,       mp_nproc(grid->mp),     fileIO );     \
  /* Species parameters */                                  \
  F_WRITE( int,       sp_id,                  fileIO );     \
  F_WRITE( float,     q_m,                    fileIO );     \
} END_PRIMITIVE
 
// Note dim _MUST_ be a pointer to an int
 
#define F_WRITE_ARRAY_HEADER(psiz,ndim,dim,fileIO) BEGIN_PRIMITIVE { \
  F_WRITE( int, psiz, fileIO );                         \
  F_WRITE( int, ndim, fileIO );                         \
  fileIO.write( dim, ndim );                            \
} END_PRIMITIVE
 
#define F_WRITE(type,value,fileIO) BEGIN_PRIMITIVE {    \
  type __F_WRITE_tmp = (type)(value);                   \
  fileIO.write( &__F_WRITE_tmp, 1 );                    \
} END_PRIMITIVE
 
#define F_READ(type,value,fileIO) BEGIN_PRIMITIVE {     \
  type __F_READ_tmp;                                    \
  fileIO.read( &__F_READ_tmp, 1 );                      \
  (value) = __F_READ_tmp;                               \
} END_PRIMITIVE

# define ABORT(cond) if( cond ) ERROR(( #cond ))

# define SETIVAR( V, A, S )                             \
	{                                                   \
	V=(A);                                              \
	if ( mp_rank(grid->mp)==0 )                         \
		MESSAGE(("* Modifying %s to value %d", S, A));  \
	}

# define SETDVAR( V, A, S )                             \
	{                                                   \
	V=(A);                                              \
	if ( mp_rank(grid->mp)==0 )                         \
		MESSAGE(("* Modifying %s to value %le", S, A)); \
	}

# define ITEST( V, N, A ) if ( sscanf( line, N " %d", &iarg )==1 ) SETIVAR( V, A, N );

# define DTEST( V, N, A ) if ( sscanf( line, N " %le", &darg )==1 ) SETDVAR( V, A, N );

#define MAX0(A,B) ((A) > (B) ? (A) : (B))
#define MIN0(A,B) ((A) < (B) ? (A) : (B))
 
#endif // dumpmacros_h
