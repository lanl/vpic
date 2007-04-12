/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#ifndef _util_base_h_
#define _util_base_h_

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif

#include <stdlib.h> /* For malloc, realloc, free, size_t, NULL */
#include <string.h> /* For string and memory manipulation */
#include <math.h>   /* For math prototypes */
#include <stdint.h> /* For fixed width integer types */
#include <limits.h> /* For integer limits */
#include <float.h>  /* For floating point limits */

/* PREFERRED_ALIGNMENT is the default alignment */

#ifndef PREFERRED_ALIGNMENT
#define PREFERRED_ALIGNMENT 16
#endif

/* ALIGNED indicates a pointer can be assumed to aligned */

#ifndef ALIGNED
#define ALIGNED
#endif

/* These macros facilitate doing evil tricks */

#define BEGIN_PRIMITIVE do
#define END_PRIMITIVE   while(0)
#define STRINGIFY(s)#s
#define EXPAND_AND_STRINGIFY(s)STRINGIFY(s)

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
           Could not allocate 45 bytes

   Note: Error messages are abortive but MESSAGE and WARNING are not

   FIXME: SHOULD PRINT THE WORLD RANK! */

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
  /* FIXME: SHOULD WAIT A FEW SECONDS HERE! */             \
  exit(1);                                                 \
} END_PRIMITIVE

/* FIXME: DEPRECATE THIS */
enum {
  preferred_alignment = PREFERRED_ALIGNMENT
};

/* FIXME: DEPRECATE THIS */
typedef const char *error_code;
#define ERROR_CODE(s) \
  ((error_code)(__FILE__"("EXPAND_AND_STRINGIFY(__LINE__)"): "s))
#define NO_ERROR ((error_code)NULL)

BEGIN_C_DECLS

/* In util.c */

void * ALIGNED
malloc_aligned( size_t n, size_t a );

void
free_aligned( void * ALIGNED mem );

void
print_log( const char *fmt, ... );

/* This function returns a value to prevent the compiler from optimizing
   it away the function body.  The caller should not use it though so the
   declaration casts away the return. */

#define nanodelay(i) ((void)_nanodelay(i))
uint32_t
_nanodelay( uint32_t i );

END_C_DECLS

#endif /* _util_base_h_ */
