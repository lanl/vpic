//============================================================================//
// Written by:
//   Kevin J. Bowers, Ph.D.
//   Plasma Physics Group (X-1)
//   Applied Physics Division
//   Los Alamos National Lab
// March/April 2004 - Revised and extened from earlier V4PIC versions.
//============================================================================//

#define IN_spa

#include "../species_advance.h"

//----------------------------------------------------------------------------//
// This is the legacy thread serial version of the particle sort.
//----------------------------------------------------------------------------//

#if defined(VPIC_USE_LEGACY_SORT)

//----------------------------------------------------------------------------//
//
//----------------------------------------------------------------------------//

void
sort_p( species_t * sp )
{
  if ( !sp )
    ERROR( ( "Bad args" ) );

  sp->last_sorted = sp->g->step;

  particle_t * ALIGNED(128) p = sp->p;
  size_t * ALIGNED(128) p_id = sp->p_id;

  const int np                = sp->np;
  const int nc                = sp->g->nv;
  const int nc1               = nc + 1;

  int * RESTRICT ALIGNED(128) partition = sp->partition;

  static int * RESTRICT ALIGNED(128) next = NULL;

  static int max_nc1 = 0;

  int i, j;

  // Do not need to sort.
  if ( np == 0 )
    return;

  // Allocate the sorting intermediate. Making this into a static is done to
  // avoid heap shredding.

  if ( max_nc1 < nc1 )
  {
    // Hack around RESTRICT issues.
    int *tmp = next;

    FREE_ALIGNED( tmp );

    MALLOC_ALIGNED( tmp, nc1, 128 );

    next    = tmp;
    max_nc1 = nc1;
  }

  // Count particles in each cell.
  CLEAR( next, nc1 );

  for( i = 0; i < np; i++ )
  {
    next[ p[i].i ]++;
  }

  // Convert the count to a partitioning and save a copy in next.
  j = 0;
  for( i = 0; i < nc1; i++ )
  {
    partition[i]  = j;
    j            += next[i];
    next[i]       = partition[i];
  }

  if ( sp->sort_out_of_place )
  {
    // Throw down the particle array in order.

    /**/  particle_t *          ALIGNED(128) new_p;
    const particle_t * RESTRICT ALIGNED( 32)  in_p;
    /**/  particle_t * RESTRICT ALIGNED( 32) out_p;

    const size_t* RESTRICT ALIGNED( 32)  in_p_id;
    /**/  size_t* RESTRICT ALIGNED( 32) out_p_id;

    MALLOC_ALIGNED( new_p, sp->max_np, 128 );
    MALLOC_ALIGNED( new_p_id, sp->max_np, 128 );

    in_p  = sp->p;
    in_p_id  = sp->p_id;

    out_p = new_p;
    out_p_id = new_p_id;

    for( i = 0; i < np; i++ )
    {
      out_p[ next[ in_p[i].i ] ] = in_p[i];
      out_p_id[ next[ in_p[i].i ] ] = in_p_id[i];
      next[ in_p[i].i ]++;  /* advance to next free slot for this cell */
    }

    FREE_ALIGNED( sp->p );
    FREE_ALIGNED( sp->p_id );

    sp->p = new_p;
    sp->p_id = new_p_id;
  }

  else
  {
    // Run sort cycles until the list is sorted.

    particle_t               save_p;
    particle_t * ALIGNED(32) src;
    particle_t * ALIGNED(32) dest;
    size_t save_pid;
    size_t * ALIGNED(32) srcid;
    size_t * ALIGNED(32) destid;

    i = 0;
    while( i < nc )
    {
      if ( next[i] >= partition[i+1] )
      {
        i++;
      }

      else
      {
        src = &p[ next[i] ];
        srcid = &p_id[ next[i] ];

        for( ; ; )
        {
          dest = &p[ next[ src->i ] ];
          destid = &p_id[ next[ src->i ] ];
          next[ src->i ]++;  /* advance to next free slot for this cell */

          if ( src == dest ) break;

          save_p = *dest;
          *dest  = *src;
          *src   = save_p;
          save_pid = *destid;
          *destid = *srcid;
          *srcid = save_pid;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------//
// This is the new thread parallel version of the particle sort.
//----------------------------------------------------------------------------//

#else

//----------------------------------------------------------------------------//
// Top level function to select and call the proper sort_p function using the
// desired particle sort abstraction.  Currently, the only abstraction
// available is the pipeline abstraction.
//----------------------------------------------------------------------------//

void
sort_p( species_t * sp )
{
  if ( ! sp )
  {
    ERROR( ( "Bad args." ) );
  }

  // Conditionally execute this when more abstractions are available.
  sort_p_pipeline( sp );
}

#endif
