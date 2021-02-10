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
  #ifdef VPIC_GLOBAL_PARTICLE_ID
  const int sp_has_ids       = sp->has_ids;
  size_t * ALIGNED(128) p_id = sp->p_id;
  #endif
  #ifdef VPIC_PARTICLE_ANNOTATION
  typedef VPIC_PARTICLE_ANNOTATION annotation_t;
  const int sp_has_annotation = sp->has_annotation;
  annotation_t* ALIGNED(128) p_annotation = sp->p_annotation;
  #endif

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

    MALLOC_ALIGNED( new_p, sp->max_np, 128 );
    in_p  = sp->p;
    out_p = new_p;

    #ifdef VPIC_GLOBAL_PARTICLE_ID
    /**/  size_t*          ALIGNED(128) new_p_id;
    const size_t* RESTRICT ALIGNED( 32)  in_p_id;
    /**/  size_t* RESTRICT ALIGNED( 32) out_p_id;

    if(sp_has_ids) {
      MALLOC_ALIGNED( new_p_id, sp->max_np, 128 );
      in_p_id  = sp->p_id;
      out_p_id = new_p_id;
    }
    #endif
    #ifdef VPIC_PARTICLE_ANNOTATION
    /**/  annotation_t*          ALIGNED(128) new_p_annotation;
    const annotation_t* RESTRICT ALIGNED( 32)  in_p_annotation;
    /**/  annotation_t* RESTRICT ALIGNED( 32) out_p_annotation;

    if(sp_has_annotation) {
      MALLOC_ALIGNED( new_p_annotation, sp->max_np*sp_has_annotation, 128 );
      in_p_annotation  = sp->p_annotation;
      out_p_annotation = new_p_annotation;
    }
    #endif

    for( i = 0; i < np; i++ )
    {
      out_p[ next[ in_p[i].i ] ] = in_p[i];
      #ifdef VPIC_GLOBAL_PARTICLE_ID
      if(sp_has_ids) {
        out_p_id[ next[ in_p[i].i ] ] = in_p_id[i];
      }
      #endif
      #ifdef VPIC_PARTICLE_ANNOTATION
      if(sp_has_annotation) {
        for(int a = 0; a < sp_has_annotation; a++) {
         out_p_annotation[ next[ in_p[i].i ]*sp_has_annotation + a ] = in_p_annotation[i*sp_has_annotation+a];
        }
      }
      #endif
      next[ in_p[i].i ]++;  /* advance to next free slot for this cell */
    }

    FREE_ALIGNED( sp->p );
    sp->p = new_p;

    #ifdef VPIC_GLOBAL_PARTICLE_ID
    if(sp_has_ids) {
      FREE_ALIGNED( sp->p_id );
      sp->p_id = new_p_id;
    }
    #endif
    #ifdef VPIC_PARTICLE_ANNOTATION
    if(sp_has_annotation) {
      FREE_ALIGNED( sp->p_annotation );
      sp->p_annotation  = new_p_annotation;
    }
    #endif
  }

  else
  {
    // Run sort cycles until the list is sorted.

    particle_t               save_p;
    particle_t * ALIGNED(32) src;
    particle_t * ALIGNED(32) dest;
    #ifdef VPIC_GLOBAL_PARTICLE_ID
    size_t save_pid;
    size_t * ALIGNED(32) srcid;
    size_t * ALIGNED(32) destid;
    #endif
    #ifdef VPIC_PARTICLE_ANNOTATION
    annotation_t save_p_annotation;
    annotation_t* ALIGNED(32) src_annotation;
    annotation_t* ALIGNED(32) dest_annotation;
    #endif

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
        #ifdef VPIC_GLOBAL_PARTICLE_ID
        if(sp_has_ids) {
          srcid = &p_id[ next[i] ];
        }
        #endif
        #ifdef VPIC_PARTICLE_ANNOTATION
        if(sp_has_annotation) {
          src_annotation = & (p_annotation[ next[i]*sp_has_annotation] );
        }
        #endif

        for( ; ; )
        {
          dest = &p[ next[ src->i ] ];
          #ifdef VPIC_GLOBAL_PARTICLE_ID
          if(sp_has_ids) {
            destid = &p_id[ next[ src->i ] ];
          }
          #endif
          #ifdef VPIC_PARTICLE_ANNOTATION
          if(sp_has_annotation) {
            for(int a = 0; a < sp_has_annotation; a++) {
              dest_annotation = & (p_annotation[ next[ src->i ]*sp_has_annotation] );
           }
          }
          #endif
          next[ src->i ]++;  /* advance to next free slot for this cell */

          if ( src == dest ) break;

          save_p = *dest;
          *dest  = *src;
          *src   = save_p;
          #ifdef VPIC_GLOBAL_PARTICLE_ID
          if(sp_has_ids) {
            save_pid = *destid;
            *destid = *srcid;
            *srcid = save_pid;
          }
          #endif
          #ifdef VPIC_PARTICLE_ANNOTATION
          if(sp_has_annotation) {
            for(int a = 0; a < sp_has_annotation; a++) {
              save_p_annotation  = dest_annotation[a];
              dest_annotation[a] = src_annotation[a];
              src_annotation[a]  = save_p_annotation;
           }
          }
          #endif
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
