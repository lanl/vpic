#ifndef dumpmacros_h
#define dumpmacros_h

/* FIXME: WHEN THESE MACROS WERE HOISTED AND VARIOUS HACKS DONE TO THEM
   THEY BECAME _VERY_ _DANGEROUS. */

#define WRITE_HEADER_V0(dump_type,sp_id,q_m,cstep,fileIO) do { \
    /* Binary compatibility information */                     \
    WRITE( char,           CHAR_BIT,               fileIO );   \
    WRITE( char,           sizeof(short int),      fileIO );   \
    WRITE( char,           sizeof(int),            fileIO );   \
    WRITE( char,           sizeof(float),          fileIO );   \
    WRITE( char,           sizeof(double),         fileIO );   \
    WRITE( unsigned short, 0xcafe,                 fileIO );   \
    WRITE( unsigned int,   0xdeadbeef,             fileIO );   \
    WRITE( float,          1.0,                    fileIO );   \
    WRITE( double,         1.0,                    fileIO );   \
    /* Dump type and header format version */                  \
    WRITE( int,            0 /* Version */,        fileIO );   \
    WRITE( int,            dump_type,              fileIO );   \
    /* High level information */                               \
    WRITE( int,            cstep,                  fileIO );   \
    WRITE( int,            nxout,                  fileIO );   \
    WRITE( int,            nyout,                  fileIO );   \
    WRITE( int,            nzout,                  fileIO );   \
    WRITE( float,          grid->dt,               fileIO );   \
    WRITE( float,          dxout,                  fileIO );   \
    WRITE( float,          dyout,                  fileIO );   \
    WRITE( float,          dzout,                  fileIO );   \
    WRITE( float,          grid->x0,               fileIO );   \
    WRITE( float,          grid->y0,               fileIO );   \
    WRITE( float,          grid->z0,               fileIO );   \
    WRITE( float,          grid->cvac,             fileIO );   \
    WRITE( float,          grid->eps0,             fileIO );   \
    WRITE( float,          0 /* damp */,           fileIO );   \
    WRITE( int,            rank(),                 fileIO );   \
    WRITE( int,            nproc(),                fileIO );   \
    /* Species parameters */                                   \
    WRITE( int,            sp_id,                  fileIO );   \
    WRITE( float,          q_m,                    fileIO );   \
  } while(0)
 
// Note dim _MUST_ be a pointer to an int
 
#define WRITE_ARRAY_HEADER(p,ndim,dim,fileIO) do { \
    WRITE( int, sizeof(p[0]), fileIO );            \
    WRITE( int, ndim,         fileIO );            \
    fileIO.write( dim, ndim );                     \
  } while(0)
 
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
 
#define WRITE(type,value,fileIO) do { \
    type __WRITE_tmp = (type)(value); \
    fileIO.write( &__WRITE_tmp, 1 );  \
  } while(0)
 
// Note: strlen does not include the terminating \0
#define WRITE_STRING(string,fileIO) do {                    \
    int __WRITE_STRING_len = 0;                             \
    if( string ) __WRITE_STRING_len = strlen(string);       \
    fileIO.write( &__WRITE_STRING_len, 1 );                 \
    if( __WRITE_STRING_len>0 )                              \
      fileIO.write( string, __WRITE_STRING_len );           \
  } while(0)
 
#define READ(type,value,fileIO) do { \
    type __READ_tmp;                 \
    fileIO.read(&__READ_tmp, 1 );    \
    (value) = __READ_tmp;            \
  } while(0)

#define F_WRITE_HEADER_V0(dump_type,sp_id,q_m,fileIO) do { \
    /* Binary compatibility information */                 \
    F_WRITE( char,      CHAR_BIT,               fileIO );  \
    F_WRITE( char,      sizeof(short int),      fileIO );  \
    F_WRITE( char,      sizeof(int),            fileIO );  \
    F_WRITE( char,      sizeof(float),          fileIO );  \
    F_WRITE( char,      sizeof(double),         fileIO );  \
    F_WRITE( short int, 0xcafe,                 fileIO );  \
    F_WRITE( int,       0xdeadbeef,             fileIO );  \
    F_WRITE( float,     1.0,                    fileIO );  \
    F_WRITE( double,    1.0,                    fileIO );  \
    /* Dump type and header format version */              \
    F_WRITE( int,       0 /* Version */,        fileIO );  \
    F_WRITE( int,       dump_type,              fileIO );  \
    /* High level information */                           \
    F_WRITE( int,       step(),                 fileIO );  \
    F_WRITE( int,       imxstr-2,               fileIO );  \
    F_WRITE( int,       jmxstr-2,               fileIO );  \
    F_WRITE( int,       kmxstr-2,               fileIO );  \
    F_WRITE( float,     grid->dt,               fileIO );  \
    F_WRITE( float,     dxstr,                  fileIO );  \
    F_WRITE( float,     dystr,                  fileIO );  \
    F_WRITE( float,     dzstr,                  fileIO );  \
    F_WRITE( float,     grid->x0,               fileIO );  \
    F_WRITE( float,     grid->y0,               fileIO );  \
    F_WRITE( float,     grid->z0,               fileIO );  \
    F_WRITE( float,     grid->cvac,             fileIO );  \
    F_WRITE( float,     grid->eps0,             fileIO );  \
    F_WRITE( float,     0 /*damp*/,             fileIO );  \
    F_WRITE( int,       rank(),                 fileIO );  \
    F_WRITE( int,       nproc(),                fileIO );  \
    /* Species parameters */                               \
    F_WRITE( int,       sp_id,                  fileIO );  \
    F_WRITE( float,     q_m,                    fileIO );  \
  } while(0)
 
#define F_WRITE_HEADER_PAR(dump_type,sp_id,q_m,fileIO) do { \
    /* Binary compatibility information */                  \
    F_WRITE( char,      CHAR_BIT,               fileIO );   \
    F_WRITE( char,      sizeof(short int),      fileIO );   \
    F_WRITE( char,      sizeof(int),            fileIO );   \
    F_WRITE( char,      sizeof(float),          fileIO );   \
    F_WRITE( char,      sizeof(double),         fileIO );   \
    F_WRITE( short int, 0xcafe,                 fileIO );   \
    F_WRITE( int,       0xdeadbeef,             fileIO );   \
    F_WRITE( float,     1.0,                    fileIO );   \
    F_WRITE( double,    1.0,                    fileIO );   \
    /* Dump type and header format version */               \
    F_WRITE( int,       0 /* Version */,        fileIO );   \
    F_WRITE( int,       dump_type,              fileIO );   \
    /* High level information */                            \
    F_WRITE( int,       step(),                 fileIO );   \
    F_WRITE( int,       grid->nx,               fileIO );   \
    F_WRITE( int,       grid->ny,               fileIO );   \
    F_WRITE( int,       grid->nz,               fileIO );   \
    F_WRITE( float,     grid->dt,               fileIO );   \
    F_WRITE( float,     grid->dx,               fileIO );   \
    F_WRITE( float,     grid->dy,               fileIO );   \
    F_WRITE( float,     grid->dz,               fileIO );   \
    F_WRITE( float,     grid->x0,               fileIO );   \
    F_WRITE( float,     grid->y0,               fileIO );   \
    F_WRITE( float,     grid->z0,               fileIO );   \
    F_WRITE( float,     grid->cvac,             fileIO );   \
    F_WRITE( float,     grid->eps0,             fileIO );   \
    F_WRITE( float,     0 /*damp*/,             fileIO );   \
    F_WRITE( int,       rank(),                 fileIO );   \
    F_WRITE( int,       nproc(),                fileIO );   \
    /* Species parameters */                                \
    F_WRITE( int,       sp_id,                  fileIO );   \
    F_WRITE( float,     q_m,                    fileIO );   \
  } while(0)
 
// Note dim _MUST_ be a pointer to an int
 
#define F_WRITE_ARRAY_HEADER(psiz,ndim,dim,fileIO) do { \
    F_WRITE( int, psiz, fileIO );                       \
    F_WRITE( int, ndim, fileIO );                       \
    fileIO.write( dim, ndim );                          \
  } while(0)
 
#define F_WRITE(type,value,fileIO) do { \
    type __F_WRITE_tmp = (type)(value); \
    fileIO.write( &__F_WRITE_tmp, 1 );  \
  } while(0)
 
#define F_READ(type,value,fileIO) do { \
    type __F_READ_tmp;                 \
    fileIO.read( &__F_READ_tmp, 1 );   \
    (value) = __F_READ_tmp;            \
  } while(0)

#define ABORT(cond) if( cond ) ERROR(( #cond ))

#endif // dumpmacros_h
