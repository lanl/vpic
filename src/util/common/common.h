/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#ifndef _common_h_
#define _common_h_

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif

#include <stdlib.h> /* For malloc, realloc, free, size_t, NULL */
#include <limits.h> /* For LONG_MAX */

/* 16-bit ints */
#ifndef INT16_TYPE
#define INT16_TYPE short int
#endif

/* 32-bit ints */
#ifndef INT32_TYPE
#define INT32_TYPE int
#endif

#ifndef INT64_TYPE
#if LONG_MAX==2147483647L
/* Long int is 32-bits ... assume long long int is 64-bits */
/* Warning: This is not ANSI-C89 compliant */
__extension__
typedef long long int int64;
#define INT64_TYPE int64
#else
/* Long int is not 32-bits ... assume long int is a 64-bits */
#define INT64_TYPE long int
#endif
#endif

/* PREFERRED_ALIGNMENT is the default alignment */
#ifndef PREFERRED_ALIGNMENT
#define PREFERRED_ALIGNMENT 16
#endif

/* RESTRICT indicates memory accessed through a pointer will not be accessed
   by other means in the scope of the pointer */
#ifndef RESTRICT
#define RESTRICT
#endif

/* ALIGNED indicates a pointer can be assumed to aligned */
#ifndef ALIGNED
#define ALIGNED
#endif

/* These macros facillitate doing evil tricks */
#define BEGIN_PRIMITIVE do
#define END_PRIMITIVE   while(0)
#define STRINGIFY(s)#s
#define EXPAND_AND_STRINGIFY(s)STRINGIFY(s)
#define CONCAT3(a,b,c)a/**/b/**/c

/* INDEX_FORTRAN_x and INDEX_C_x give macros for accessing multi-dimensional
   arrays with different conventions. To eliminate potential side effects and
   maximize optimization possibilites, xl, xh, yl, yh, zl, zh should be
   local constant ints */
#define INDEX_FORTRAN_1(x,xl,xh)                 \
 ((x)-(xl))
#define INDEX_FORTRAN_2(x,y,xl,xh,yl,yh)         \
 ((x)-(xl) + ((xh)-(xl)+1)*((y)-(yl)))
#define INDEX_FORTRAN_3(x,y,z,xl,xh,yl,yh,zl,zh) \
 ((x)-(xl) + ((xh)-(xl)+1)*((y)-(yl) + ((yh)-(yl)+1)*((z)-(zl))))

#define INDEX_C_1(x,xl,xh)                 \
 ((x)-(xl))
#define INDEX_C_2(x,y,xl,xh,yl,yh)         \
 ((y)-(yl) + ((yh)-(yl)+1)*((x)-(xl)))
#define INDEX_C_3(x,y,z,xl,xh,yl,yh,zl,zh) \
 ((z)-(zl) + ((zh)-(zl)+1)*((y)-(yl) + ((yh)-(yl)+1)*((x)-(xl))))

/* The following macros deal with linked lists */
#define LIST_FOR_EACH(node,list)        \
  for((node)=(list); (node)!=NULL; (node)=(node)->next)
#define LIST_FIND_FIRST(node,list,cond) \
  for((node)=(list); (node)!=NULL; (node)=(node)->next) if(cond) break

/* The following macros give a provide a simple logging capabilty. Due to the
   way they work, usage needs double parenthesis. That is:

     ERROR(("Could not allocate %i bytes", req));

   Will print the following message to the log:

   Error at src/module/file.c(34):
           Could not allocate 45 bytes */

#define CHECKPOINT() BEGIN_PRIMITIVE {		           \
  print_log( "%s(%i): Checkpoint\n", __FILE__, __LINE__ ); \
} END_PRIMITIVE

#define MESSAGE(args) BEGIN_PRIMITIVE {		\
  print_log( "%s(%i): ", __FILE__, __LINE__ );	\
  print_log args;			        \
  print_log( "\n" );                            \
} END_PRIMITIVE

#define WARNING(args) BEGIN_PRIMITIVE {			     \
  print_log( "Warning at %s(%i):\n\t", __FILE__, __LINE__ ); \
  print_log args;                                            \
  print_log( "\n" );                                         \
} END_PRIMITIVE

#define ERROR(args) BEGIN_PRIMITIVE {			   \
  print_log( "Error at %s(%i):\n\t", __FILE__, __LINE__ ); \
  print_log args;					   \
  print_log( "\n" );                                       \
} END_PRIMITIVE

enum common_enums {
  preferred_alignment = PREFERRED_ALIGNMENT
};

typedef const char *error_code;
#define ERROR_CODE(s) \
  ((error_code)(__FILE__"("EXPAND_AND_STRINGIFY(__LINE__)"): "s))
#define SUCCESS ((error_code)NULL)

BEGIN_C_DECLS

/* In common.c */
extern void * ALIGNED malloc_aligned( size_t n, size_t a );
extern void free_aligned( void * ALIGNED mem );
extern void print_log( const char *fmt, ... );

END_C_DECLS

#endif
