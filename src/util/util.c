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
#include <unistd.h>

void
util_malloc( const char * err,
             void * mem_ref, 
             size_t n ) {
  char * mem;

  // If no err given, use a default error 
  if( err==NULL ) err = "malloc failed (n=%lu)";

  // Check that mem_ref is valid
  if( mem_ref==NULL ) ERROR(( err, (unsigned long)n ));

  // A do nothing request
  if( n==0 ) { *(char **)mem_ref = NULL; return; }

  // Allocate the memory ... abort if the allocation fails
  mem = (char *)malloc(n);
  if( mem==NULL ) ERROR(( err, (unsigned long)n ));
  *(char **)mem_ref = mem;
}

void
util_free( void * mem_ref ) {
  char * mem;
  if( mem_ref==NULL ) return;
  mem = *(char **)mem_ref;
  if( mem!=NULL ) free( mem );
  *(char **)mem_ref = NULL;
}

// FIXME TEMPORARY HACK
// FIXME: This is a hack to make the current processor rank available to the
// Error output messages during memory allocation.  This should never
// appear in production code.
int err_rank;

void
util_malloc_aligned( const char * err,
                     void * mem_ref, 
                     size_t n,
                     size_t a ) {
  char *mem_u, *mem_a, **mem_p;

  // If no err given, use a default error 
  if( err==NULL ) err = "malloc aligned failed (n=%lu, a=%lu)";

  // Check that mem_ref is valid and a is a power of two 
  if( mem_ref==NULL || a==0 || (a&(a-1))!=0 )
    ERROR(( err, (unsigned long)n, (unsigned long)a ));

  // A do nothing request 
  if( n==0 ) { *(char **)mem_ref = NULL; return; }

  // Adjust small alignments to a minimal valid alignment 
  // and convert a into a mask of the address LSB
  if( a<16 ) a = 16;
  a--;

  // Allocate the raw unaligned memory ... abort if the allocation fails 
  mem_u = (char *)malloc( n + a + sizeof(char *) );
  if( mem_u==NULL ) {

// This is a temporary diagnostic for use in debugging Roadrunner memory
// allocation errors.  This ensures that, in the case of memory allocation
// failure, the code will hang without exiting, so that node state can be
// determined.
#define TEMPORARY_HACK 1
#if TEMPORARY_HACK
    char hostname[256];
	gethostname(hostname, 256);
	WARNING(("Memory allocation failed on rank %d (host %s)", err_rank,
	  hostname));
    while(1) {
	  sleep(1);
	} // while
#endif

    ERROR(( err, (unsigned long)n, (unsigned long)a ));
  } // if

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
  if( mem_ref==NULL ) return;
  mem_a = *(char **)mem_ref;
  if( mem_a!=NULL ) {
    mem_p = (char **)(mem_a - sizeof(char *));
    mem_u = mem_p[0];
    free( mem_u );
  }
  *(char **)mem_ref = NULL;
}

void print_log( const char *fmt, ... ) {
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

