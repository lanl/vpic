/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <stdio.h>    /* For vfprintf */
#include <stdarg.h>   /* For va_list, va_start, va_end */
#include <stdlib.h>   /* For malloc, realloc, free */

#include <common.h>   /* Declarations */

void * ALIGNED malloc_aligned( size_t n, size_t a ) {
  char *mem_u, *mem_a, **mem_p;

  a--;
  mem_u = (char *)malloc( n + a + sizeof(char *) );
  if( mem_u==NULL ) {
    mem_a = NULL;
  } else {
    mem_a = (char *)(((unsigned long int)(mem_u + a + sizeof(char *)))&(~a));
    mem_p = (char **)(mem_a - sizeof(char *));
    mem_p[0] = mem_u;
  }

  return (void * ALIGNED)mem_a;
}

void free_aligned( void * ALIGNED mem ) {
  char *mem_u, *mem_a, **mem_p;

  if( mem!=NULL ) {
    mem_a = (char *)mem;
    mem_p = (char **)(mem_a - sizeof(char *));
    mem_u = mem_p[0];
    free(mem_u);
  }
  return;
}

void print_log( const char *fmt, ... ) {
  va_list ap;
  va_start( ap, fmt );
  vfprintf( stderr, fmt, ap );
  va_end( ap );
  fflush( stderr );
}
