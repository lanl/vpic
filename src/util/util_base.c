/*
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "util_base.h" // Declarations
#include <stdio.h>     // For vfprintf
#include <stdarg.h>    // For va_list, va_start, va_end

/****************************************************************************/

#define STRIP_CMDLINE( what, T, convert )                         \
T                                                                 \
strip_cmdline_##what( int * pargc,                                \
                      char *** pargv,                             \
                      const char * key,                           \
                      T val ) {                                   \
  int i, n = 0;                                                   \
  for( i=0; i<(*pargc); i++ )                                     \
    if( strcmp( (*pargv)[i], key ) ) (*pargv)[n++] = (*pargv)[i]; \
    else if( (++i)<(*pargc) ) val = convert( (*pargv)[i] );       \
  (*pargv)[n] = NULL, (*pargc) = n; /* ANSI -argv is NULL terminated */ \
  return val;                                                     \
}

int
strip_cmdline( int * pargc,
               char *** pargv,
               const char * key ) {
  int i, n = 0, val = 0;
  for( i=0; i<(*pargc); i++ )
    if( strcmp( (*pargv)[i], key ) ) (*pargv)[n++] = (*pargv)[i];
    else val++;
  (*pargv)[n] = NULL, (*pargc) = n; /* ANSI - argv is NULL terminated */
  return val;
}

STRIP_CMDLINE( int,    int,          atoi )
STRIP_CMDLINE( double, double,       atof )
STRIP_CMDLINE( string, const char *,      )

#undef STRIP_CMDLINE

/****************************************************************************/

void
util_malloc( const char * err,
             void * mem_ref,
             size_t n ) {
  char * mem;

  // If no err given, use a default error
  if( !err ) err = "malloc failed (n=%lu)";

  // Check that mem_ref is valid
  if( !mem_ref ) ERROR(( err, (unsigned long)n ));

  // A do nothing request
  if( n==0 ) { *(char **)mem_ref = NULL; return; }

  // Allocate the memory ... abort if the allocation fails
  mem = (char *)malloc(n);
  if( !mem ) ERROR(( err, (unsigned long)n ));
  *(char **)mem_ref = mem;
}

void
util_free( void * mem_ref ) {
  char * mem;
  if( !mem_ref ) return;
  mem = *(char **)mem_ref;
  if( mem ) free( mem );
  *(char **)mem_ref = NULL;
}

void
util_malloc_aligned( const char * err,
                     void * mem_ref,
                     size_t n,
                     size_t a ) {
  char *mem_u, *mem_a, **mem_p;

  // If no err given, use a default error
  if( !err ) err = "malloc aligned failed (n=%lu, a=%lu)";

  // Check that mem_ref is valid and a is a power of two
  if( !mem_ref || a==0 || (a&(a-1))!=0 )
    ERROR(( err, (unsigned long)n, (unsigned long)a ));

  // A do nothing request
  if( n==0 ) { *(char **)mem_ref = NULL; return; }

  // Adjust small alignments to a minimal valid alignment
  // and convert a into a mask of the address LSB
  if( a<16 ) a = 16;
  a--;

  // Allocate the raw unaligned memory ... abort if the allocation fails
  mem_u = (char *)malloc( n + a + sizeof(char *) );
  if( !mem_u ) ERROR(( err, (unsigned long)n, (unsigned long)a ));

  // Compute the pointer to the aligned memory and save a pointer to the
  // raw unaligned memory for use on free_aligned
  mem_a = (char *)(((unsigned long int)(mem_u + a + sizeof(char *)))&(~a));
  mem_p = (char **)(mem_a - sizeof(char *));
  mem_p[0] = mem_u;

  *(char **)mem_ref = mem_a;
}

void
util_free_aligned( void * mem_ref ) {
  char *mem_u, *mem_a, **mem_p;
  if( !mem_ref ) return;
  mem_a = *(char **)mem_ref;
  if( mem_a ) {
    mem_p = (char **)(mem_a - sizeof(char *));
    mem_u = mem_p[0];
    free( mem_u );
  }
  *(char **)mem_ref = NULL;
}

/*****************************************************************************/

void
log_printf( const char *fmt, ... ) {
  va_list ap;
  va_start( ap, fmt );
  vfprintf( stderr, fmt, ap );
  va_end( ap );
  fflush( stderr );
}

uint32_t
_nanodelay( uint32_t i ) {
  uint32_t a = 0;
  for( ; i; i-- ) a^=0xdeadbeef, a>>=1;
  return a;
}

