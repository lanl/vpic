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

// C99 does requires some key macros of stdint to only be defined in
// C++ implementations if explicitly requested. 

#define __STDC_LIMIT_MACROS

// In a stunning language flaw, C++98 does not support 64-bit data
// types (including no standard way to initialize a non-trivial 64-bit
// type to a compile time constant).  Hilary ensures when stdint gives
// you 64-bit primitive data types.  C99's stdint defines some macros
// that fix this if explicitly requested.

#define __STDC_CONSTANT_MACROS

#include <stdlib.h> // For exit, size_t, NULL
#include <string.h> // For string and memory manipulation
#include <stdint.h> // For fixed width integer types
#include <math.h>   // For math prototypes
#include <limits.h> // For integer limits
#include <float.h>  // For floating point limits

// Opaque handle of a  set of communicating processes

struct collective;
typedef struct collective collective_t;

// These macros facilitate doing evil tricks

#define BEGIN_PRIMITIVE do
#define END_PRIMITIVE   while(0)

#define _UTIL_STRINGIFY(s)#s
#define EXPAND_AND_STRINGIFY(s)_UTIL_STRINGIFY(s)

// Function inlining

// STATIC_INLINE declares a function as having static inline semantics
// in the C99 sense.  Define this to nothing at compile time to disable
// or to the appropriate value supported by your compiler.

#ifndef STATIC_INLINE
#define STATIC_INLINE static __inline__
#endif

// Conditional annotations

// LIKELY and UNLIKELY allow the programmer to annotate the likelihood
// of a given conditional being true.  If NO_BRANCH_HINTS is passed
// at compile time, this optimization will be disabled.

#ifdef NO_BRANCH_HINTS
#define   LIKLEY(_c) _c
#define UNLIKLEY(_c) _c
#else
#ifndef LIKELY
#define LIKELY(_c)   __builtin_expect((_c),1)
#endif
#ifndef UNLIKELY
#define UNLIKELY(_c) __builtin_expect((_c),0)
#endif
#endif

// Pointer qualifiers

// This pointer modifier indicates that a pointer can be assumed to
// have at least the given power-of-two alignment.
//
// May want to use compiler intrinsics to provide compiler
// the aligment information as well.
//
// NOTE: This is *currently* intentionally empty, but can be specialized on a
// per-platform basis
#ifndef ALIGNED
#define ALIGNED(a)
#endif

// This pointer modifier indicates that a pointer is restricted in
// the C99 sense.  This is added to allow C++ code (which technically
// does not have restricted pointers) to use this optimization.
// Define this to nothing at compile time to disable this or to the
// appropriate value supported by your compiler.

#ifndef RESTRICT
#define RESTRICT __restrict
#endif 

// Normal pointers (e.g. a *) are in whatever address space the given
// compile unit uses.  However, sometimes it is necessary to declare
// pointers that are understandable in multiple address spaces.  The
// MEM_PTR macro declares a pointer to external memory with the
// specified alignment.  This allows declaractions to be compiled on
// both the SPU and PPU with appropriate annotations to necessary
// write the appropriate DMA transfers.

# define MEM_PTR(type,align) type * ALIGNED(align)

// The SIZEOF_MEM_PTR macro gives the number of bytes taken by a MEM_PTR.

#define SIZEOF_MEM_PTR sizeof(MEM_PTR(void,1))

// DECLARE_ALIGNED_ARRAY declares an array containing count elements
// with the given alignment in memory.  The scope of the array is the
// scope of context in which it is declared.  Note: This macro is
// really two statements (there is no way to bundle the macro into one
// semantic statement linguistically without defeating the purpose of
// the macro).  Thus, it is not as robust as it could be.  Thus, sure
// any usage of this macro occurs in contexts where two back-to-back
// statments in the same context would be in the same scope.  That is:
//
//   if(...) { DECLARE_ALIGNED_ARRAY(type,align,name,count); ... } // OKAY!
//   else      DECLARE_ALIGNED_ARRAY(type,align,name,count);     // NOT OKAY!
//
// For 99.9% of the expected usage, this should not matter.
//
// align should be a power of two.

#if 0 // C99 has (dubious) issues with this
#define DECLARE_ALIGNED_ARRAY(type,align,name,count)                    \
  char _aa_##name[(count)*sizeof(type)+(align)];                        \
  type * ALIGNED(align) const name = (type * ALIGNED(align))            \
    ( ( (size_t)_aa_##name + (align) - 1 ) & (~((align)-1)) )
#else // Sigh ... this is technically not portable
#define DECLARE_ALIGNED_ARRAY(type,align,name,count)    \
  type name[(count)] __attribute__ ((aligned (align)))
#endif

// PAD(s,a) computes the amount of bytes necessary to add to "s" bytes
// to make "s" evenly divisible by "a" (a power of two).  Note: PAD is
// more wasteful than necessary if no padding needed.  Ideally PAD
// should be:
//   ( ( (a) - ( (s) & ( (a)-1 ) ) ) & ( (a)-1 ) )
// Unfortunately, C++98 does not support zero sized structure members
// and the preprocessor cannot evaluate the typical inputs to size to
// allow correct autogeneration when no alignment necessary ... sigh
// ...

#define PAD(s,a) ( (a) - ( (s) & ( (a)-1 ) ) ) 

// POW2_CEIL rounds "u" up to the nearest multiple of the power of two
// "a".  If u is a multiple of "a", its value is unchanged.  "a" should
// be safe against multiple dereferencing and the same type as "u".

#define POW2_CEIL(u,a) ( ((u)+(a)-1) & (~((a)-1)) )

// ALIGN_PTR rounds "p" up to the nearest multiple of the power of two
// "a".  If p is a multiple of "a", its value is unchanged.  "a" should
// be safe against multiple dereferencing.  The result is cast to a
// pointer of type "t".

#define ALIGN_PTR(t,p,a) ((t *)POW2_CEIL( (size_t)(p), (size_t)(a) ))

// Workload distribution macros

// Let items be enumerated 0,1, ... "N"-1 and pipelines be
// enumerated 0,1, ... "P" (0:P-1 refer to pipeline processes
// and pipeline P refers to the dispatching process).
// DISTRIBUTE determines the first item "i" and the number
// of items "n" that should be assigned to pipeline "p".
// The items are assigned such that each pipeline process
// gets an approximate equal share of items and that the
// number of items each pipeline process gets is a
// multiple of the block size "b".  The dispatching process
// is assigned all remaining items.  The items are assigned
// in monotonically increasing order.
//
// This macro is robust.  (All arguments only evaluated once;
// inputs can be same as output.)  Any compiler worth its
// salt will replace the divison and modulo with bit shifts
// and masks for power-of-two block sizes.

#define DISTRIBUTE( N, b, p, P, i, n ) BEGIN_PRIMITIVE {             \
    int _N = (N), _b = (b), _p = (p), _P = (P);                      \
    double _t = (double)(_N/_b)/(double)_P;                          \
    int _i =                    _b*(int)(_t*(double) _p   +0.5);     \
    (n) = (_p==_P) ? (_N%_b) : (_b*(int)(_t*(double)(_p+1)+0.5)-_i); \
    (i) = _i;                                                        \
  } END_PRIMITIVE

// INDEX_FORTRAN_x and INDEX_C_x give macros for accessing
// multi-dimensional arrays with different conventions. To eliminate
// potential side effects and maximize optimization possibilites, xl,
// xh, yl, yh, zl, zh should be local constant ints

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

// The following macros deal with linked lists

#define LIST_FOR_EACH(node,list)        \
  for((node)=(list); (node); (node)=(node)->next)

#define LIST_FIND_FIRST(node,list,cond) do {        \
    for((node)=(list); (node); (node)=(node)->next) \
      if(cond) break;                               \
  } while(0)

// Given an integer data type "type", MASK_BIT_RANGE returns a bit
// field of that type for which bits [f,l] inclusive are 1 and all
// other bits are zero.  Note that: 0<=f<=l<CHAR_BIT*sizeof(type) and
// f must be safe against multiple dereferencing.

#define MASK_BIT_RANGE(type,f,l) \
  ( ( ( ((type)1) << ( (l) - (f) + 1 ) ) - 1 ) << (f) )

// The following macros give a provide a simple logging capabilty. Due
// to the way they work, usage needs double parenthesis. That is:
//
//   ERROR(("Could not allocate %i bytes", req));
//
// will print the following message to the log:
//
//   Error at src/module/file.c(34):
//           Could not allocate 45 bytes
//
// Note: Error messages are abortive but MESSAGE and WARNING are not

#define _LOG_HDR __FILE__ "(" EXPAND_AND_STRINGIFY(__LINE__) ")"

#define CHECKPOINT() log_printf( _LOG_HDR"[%i]: Checkpoint\n", world_rank )

// WARNING: USE ONLY ONE LOG_PRINTF PER MESSAGE TO TRY AND MAKE THIS
// ATOMIC WHEN MULTIPLE RANKS USE SIMULTANEOUSLY

#define MESSAGE(args) do {                      \
    log_printf( _LOG_HDR "[%i]: ", world_rank ); \
    log_printf args;                            \
    log_printf( "\n" );                         \
  } while(0)

#define WARNING(args) do {                                      \
    log_printf( "Warning at " _LOG_HDR "[%i]:\n\t", world_rank ); \
    log_printf args;                                            \
    log_printf( "\n" );                                         \
  } while(0)

#define ERROR(args) do {                                      \
    log_printf( "Error at " _LOG_HDR "[%i]:\n\t", world_rank ); \
    log_printf args;                                          \
    log_printf( "\n" );                                       \
    nanodelay( 1000000000 ); /* Let the message out */        \
    exit(1);                                                  \
  } while(0)

// Element wise (rather than byte wise) mem{cpy,move,set} semantics

#define COPY(  d, s, n ) do { size_t _sz = (n)*sizeof(*(d)); if( _sz>0 ) memcpy(  (d), (s), _sz ); } while(0)
#define MOVE(  d, s, n ) do { size_t _sz = (n)*sizeof(*(d)); if( _sz>0 ) memmove( (d), (s), _sz ); } while(0)
#define CLEAR( d,    n ) do { size_t _sz = (n)*sizeof(*(d)); if( _sz>0 ) memset(  (d),   0, _sz ); } while(0)

BEGIN_C_DECLS

// These variables indicate the world communicator, the number of
// processes in it and the rank of this process in it.
// (The macros turn these into rvals that can't be modified
// by users accidentically).

#define world      ((collective_t *)_world)
extern collective_t * _world;

#define world_size ((int)_world_size)
extern int _world_size;

#define world_rank ((int)_world_rank)
extern int _world_rank;

// Strip all instances of key from the command line. Returns the
// number of times key was found.

int
strip_cmdline( int * pargc,
               char *** pargv,
               const char * key );

// Strip all instances of "key val" from the command line.  Returns
// val as an int of the last complete "key val" pair (if the last
// argument on the command line is "key", it too stripped but
// otherwise ignored).  If there are no instances of "key val" on
// the command line, returns default_val.

int
strip_cmdline_int( int * pargc,
                   char *** pargv,
                   const char * key,
                   int default_val );

// Same as strip_cmdline_int, but for doubles

double
strip_cmdline_double( int * pargc,
                      char *** pargv,
                      const char * key,
                      double default_val );

// Same as strip_cmdline_int, but for strings.  The lifetime of the
// returned '\0'-terminated string is the shorter of the lifetime of
// default_val or pargv.

const char *
strip_cmdline_string( int * pargc,
                      char *** pargv,
                      const char * key,
                      const char * default_val );

// In util.c
void detect_old_style_arguments(int* pargc, char *** pargv);

// MALLOC is guaranteed to succeed from the caller's point of view
// (thus, _no_ NULL checking the pointer is necessary).  n is the
// number of elements of the type of x to allocate (_not_ the number
// of bytes to allocate).  n==0 is a request for no elements and x is
// set NULL as a result.

#define MALLOC(x,n)                                                    \
  util_malloc( "MALLOC( "#x", "#n" (%lu bytes) ) at "                  \
               __FILE__ "(" EXPAND_AND_STRINGIFY(__LINE__) ") failed", \
               &(x), (n)*sizeof(*(x)) ) 

void
util_malloc( const char * err_fmt, // Has exactly one %lu in it
             void * mem_ref,
             size_t n );

// FREE frees memory allocated via MALLOC above.  It is safe to pass
// any value returned by MALLOC to FREE (_including_ a null pointer).
// The pointer to the memory will be set to NULL to indicate that it
// no longer points to anything.

#define FREE(x) util_free(&(x))

void
util_free( void * mem_ref );

// MALLOC_ALIGNED behaves equivalently to MALLOC.  The alignment must
// be a power of two.  Alignments smaller than 16 will be rounded up
// to 16.

#define MALLOC_ALIGNED(x,n,a)                                                  \
  util_malloc_aligned( "MALLOC_ALIGNED( "#x", "                                \
                                         #n" (%lu bytes), "                    \
                                         #a" (%lu bytes) ) at "                \
                       __FILE__ "(" EXPAND_AND_STRINGIFY(__LINE__) ") failed", \
                       &(x), (n)*sizeof(*(x)), (a) ) 


void
util_malloc_aligned( const char * err_fmt, // Has exactly two %lu in it
                     void * mem_ref,
                     size_t n,
                     size_t a );

// FREE_ALIGNED behaves equivalently to FREE.

#define FREE_ALIGNED(x) util_free_aligned(&(x))

void
util_free_aligned( void * mem_ref );

void
log_printf( const char *fmt, ... );

// This function returns a value to prevent the compiler from
// optimizing it away the function body.  The caller should not use it
// though so the declaration casts away the return.

#define nanodelay(i) ((void)_nanodelay(i))
uint32_t
_nanodelay( uint32_t i );

END_C_DECLS

#endif // _util_base_h_
