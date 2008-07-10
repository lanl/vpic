/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised and extened from earlier V4PIC versions
 *
 */

#define IN_spa
#include "spa_private.h"

// FIXME: V4 ACCELERATE THIS!

void
sort_p( species_t * sp,
        const grid_t * g ) {
  particle_t * ALIGNED(128) p = sp->p;
  const int32_t * __restrict ALIGNED(128) sfc = g->sfc;
  const int np                = sp->np; 
  const int nc                = (g->nx+2)*(g->ny+2)*(g->nz+2);
  const int nc1               = nc+1;
  int * __restrict ALIGNED(128) partition = sp->partition;

  static int * __restrict ALIGNED(128) next = NULL;
  static int max_nc1 = 0;

  int i, j;

  // FIXME: TEMPORARY HACK UNTIL SPECIES_ADVANCE API INSTALLED
  if( sp->partition==NULL ) MALLOC_ALIGNED( sp->partition, nc1, 128 );
  partition = sp->partition;

  if( np==0 ) return; // Do not need to sort

  // Allocate the sorting intermediate
  // Making this into a static is done to avoid heap shredding
 
  if( max_nc1<nc1 ) {
    int * tmp = next; // Hack around __restrict__ issues
    FREE_ALIGNED(   tmp );
    MALLOC_ALIGNED( tmp, nc1, 128 );
    next    = tmp;
    max_nc1 = nc1;
  }

  // Count particles in each cell
  CLEAR( next, nc1 );
  for( i=0; i<np; i++ ) next[ sfc[ p[i].i ] ]++;

  // Convert the count to a partitioning (and save a copy in next)
  j=0;
  for( i=0; i<nc1; i++ ) {
    partition[i] = j;
    j += next[i];
    next[i] = partition[i];
  }

  if( sp->sort_out_of_place ) {

    // Throw down the particle array in order

    particle_t * ALIGNED(128) new_p;
    const particle_t * __restrict ALIGNED(32) in_p;
    /**/  particle_t * __restrict ALIGNED(32) out_p;

    MALLOC_ALIGNED( new_p, sp->max_np, 128 );

    in_p  = sp->p;
    out_p = new_p;
    for( i=0; i<np; i++ ) out_p[ next[ sfc[ in_p[i].i ] ]++ ] = in_p[i];

    FREE_ALIGNED( sp->p );
    sp->p = new_p;

  } else {

    // Run sort cycles until the list is sorted

    particle_t save_p, * ALIGNED(32) src, * ALIGNED(32) dest;

    i=0;
    while( i<nc ) {
      if( next[i]>=partition[i+1] ) i++;
      else {
        src = &p[ next[i] ];
        for(;;) {
          dest = &p[ next[ sfc[ src->i ] ]++ ];
          if( src==dest ) break;
          save_p = *dest;
          *dest  = *src;
          *src   = save_p;
        }
      }

    }
  }
}

