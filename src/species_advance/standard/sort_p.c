//============================================================================//
// Written by:
//   Kevin J. Bowers, Ph.D.
//   Plasma Physics Group (X-1)
//   Applied Physics Division
//   Los Alamos National Lab
// March/April 2004 - Revised and extened from earlier V4PIC versions.
//============================================================================//

#define IN_spa
#include "spa_private.h"

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

    for( i = 0; i < np; i++ )
    {
      out_p[ next[ in_p[i].i ]++ ] = in_p[i];
    }

    FREE_ALIGNED( sp->p );

    sp->p = new_p;
  }

  else
  {
    // Run sort cycles until the list is sorted.

    particle_t               save_p;
    particle_t * ALIGNED(32) src;
    particle_t * ALIGNED(32) dest;

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

        for( ; ; )
        {
          dest = &p[ next[ src->i ]++ ];

          if ( src == dest ) break;

          save_p = *dest;
          *dest  = *src;
          *src   = save_p;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------//
// This is the new thread parallel version of the particle sort.
//----------------------------------------------------------------------------//

#else

#if defined( __SSE__ )
#include "xmmintrin.h"
#endif

// FIXME: HOOK UP IN-PLACE / OUT-PLACE OPTIONS AGAIN.

//----------------------------------------------------------------------------//
// 
//----------------------------------------------------------------------------//

void
coarse_count_pipeline( sort_p_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline )
{
  const particle_t * RESTRICT ALIGNED(128) p_src = args->p;

  int i, i1;

  int n_subsort = args->n_subsort;
  int vl        = args->vl;
  int vh        = args->vh;
  int cp_stride = POW2_CEIL( n_subsort, 4 );

  // On pipeline stack to avoid cache hot spots.
  int count[256];

  // No straggler cleanup needed.
  if ( pipeline_rank == n_pipeline )
  {
    return;
  }

  if ( n_subsort > 256 )
  {
    ERROR( ( "n_subsort too large." ) );
  }

  DISTRIBUTE( args->n, 1, pipeline_rank, n_pipeline, i, i1 );

  i1 += i;

  // Clear the local coarse count.
  CLEAR( count, n_subsort );

  // Local coarse count the input particles.
  for( ; i < i1; i++ )
  {
    count[ V2P( p_src[i].i, n_subsort, vl, vh ) ]++;
  }

  // Copy local coarse count to output.
  COPY( args->coarse_partition + cp_stride*pipeline_rank,
	count,
	n_subsort );
}

//----------------------------------------------------------------------------//
// 
//----------------------------------------------------------------------------//

void
coarse_sort_pipeline( sort_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline )
{
  const particle_t * RESTRICT ALIGNED(128) p_src = args->p;
  /**/  particle_t * RESTRICT ALIGNED(128) p_dst = args->aux_p;

  int i, i1;
  int n_subsort = args->n_subsort;
  int vl        = args->vl;
  int vh        = args->vh;
  int cp_stride = POW2_CEIL( n_subsort, 4 );
  int j;

  // On pipeline stack to avoid cache hot spots and to allow reuse of coarse
  // partitioning for fine sort stage.
  int next[ 256 ];

  // No straggler cleanup needed.
  if ( pipeline_rank == n_pipeline )
  {
    return;
  }

  if ( n_subsort > 256 )
  {
    ERROR( ( "n_subsort too large." ) );
  }

  DISTRIBUTE( args->n, 1, pipeline_rank, n_pipeline, i, i1 );

  i1 += i;

  // Load the local coarse partitioning into next.
  COPY( next,
	args->coarse_partition + cp_stride*pipeline_rank,
	n_subsort );

  // Copy particles into aux array in coarse sorted order.
  for( ; i < i1; i++ )
  {
    j = next[ V2P( p_src[i].i, n_subsort, vl, vh ) ]++;

#   if defined( __SSE__ )

    _mm_store_ps( &p_dst[j].dx, _mm_load_ps( &p_src[i].dx ) );
    _mm_store_ps( &p_dst[j].ux, _mm_load_ps( &p_src[i].ux ) );

#   else

    p_dst[j] = p_src[i];

#   endif
  }
}

//----------------------------------------------------------------------------//
// 
//----------------------------------------------------------------------------//

void
subsort_pipeline( sort_p_pipeline_args_t * args,
                  int pipeline_rank,
                  int n_pipeline )
{
  const particle_t * RESTRICT ALIGNED(128) p_src = args->aux_p;
  /**/  particle_t * RESTRICT ALIGNED(128) p_dst = args->p;

  int i0, i1, v0, v1, i, j, v, sum, count;

  int subsort;

  int n_subsort = args->n_subsort;

  int * RESTRICT ALIGNED(128) partition = args->partition;
  int * RESTRICT ALIGNED(128) next      = args->next;

  // No straggler cleanup needed.
  if ( pipeline_rank == n_pipeline )
  {
    return;
  }

  for( subsort = pipeline_rank; subsort < n_subsort; subsort += n_pipeline )
  {
    // This subsort sorts particles in [i0,i1) in the aux array. These
    // particles are in voxels [v0,v1).
    i0 = args->coarse_partition[ subsort   ];
    i1 = args->coarse_partition[ subsort+1 ];

    v0 = P2V( subsort,   n_subsort, args->vl, args->vh );
    v1 = P2V( subsort+1, n_subsort, args->vl, args->vh );

    // Clear fine grained count.
    CLEAR( &next[v0], v1 - v0 );

    // Fine grained count.
    for( i = i0; i < i1; i++ )
    {
      next[ p_src[i].i ]++;
    }

    // Compute the partitioning.
    sum = i0;
    for( v = v0; v < v1; v++ )
    {
      count         = next[v];
      next[v]       = sum;
      partition[v]  = sum;
      sum          += count;
    }
    // All subsorts who write this agree.
    partition[v1] = sum;

    // Local fine grained sort.
    for( i = i0; i < i1; i++ )
    {
      v = p_src[i].i;
      j = next[v]++;

#     if defined( __SSE__ )

      _mm_store_ps( &p_dst[j].dx, _mm_load_ps( &p_src[i].dx ) );
      _mm_store_ps( &p_dst[j].ux, _mm_load_ps( &p_src[i].ux ) );

#     else

      p_dst[j] = p_src[i];

#     endif
    }
  }
}

//----------------------------------------------------------------------------//
// 
//----------------------------------------------------------------------------//

void
sort_p( species_t * sp )
{
  if ( !sp )
    ERROR( ( "Bad args" ) );

  sp->last_sorted = sp->g->step;

  static char * ALIGNED(128)     scratch = NULL;
  static size_t              max_scratch = 0;

  size_t sz_scratch;

  particle_t * RESTRICT ALIGNED(128) p = sp->p;
  particle_t * RESTRICT ALIGNED(128) aux_p;

  int n_particle = sp->np;

  int * RESTRICT ALIGNED(128) partition = sp->partition;
  int * RESTRICT ALIGNED(128) next;

  int vl = VOXEL( 1,
		  1,
		  1,
		  sp->g->nx,
		  sp->g->ny,
		  sp->g->nz );

  int vh = VOXEL( sp->g->nx,
		  sp->g->ny,
		  sp->g->nz,
		  sp->g->nx,
		  sp->g->ny,
		  sp->g->nz );

  int n_voxel = sp->g->nv;

  int * RESTRICT ALIGNED(128) coarse_partition;

  int n_pipeline = N_PIPELINE;
  int n_subsort  = N_PIPELINE;

  int cp_stride = POW2_CEIL( n_subsort, 4 );

  int i, pipeline_rank, subsort, count, sum;

  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );

  // Ensure enough scratch space is allocated for the sorting.
  sz_scratch = ( sizeof( *p ) * n_particle      +
		 128                            +
                 sizeof( *partition ) * n_voxel +
		 128                            +
                 sizeof( *coarse_partition ) * ( cp_stride * n_pipeline + 1 ) );

  if ( sz_scratch > max_scratch )
  {
    FREE_ALIGNED( scratch );

    MALLOC_ALIGNED( scratch, sz_scratch, 128 );

    max_scratch = sz_scratch;
  }

  aux_p            = ALIGN_PTR( particle_t, scratch,            128 );
  next             = ALIGN_PTR( int,        aux_p + n_particle, 128 );
  coarse_partition = ALIGN_PTR( int,        next  + n_voxel,    128 );

  // Setup pipeline arguments.
  args->p                = p;
  args->aux_p            = aux_p;
  args->coarse_partition = coarse_partition;
  args->next             = next;
  args->partition        = partition;
  args->n                = n_particle;
  args->n_subsort        = n_subsort;
  args->vl               = vl;
  args->vh               = vh;
  args->n_voxel          = n_voxel;

  if ( n_subsort != 1 )
  {
    // Do the coarse count.
    EXEC_PIPELINES( coarse_count, args, 0 );

    WAIT_PIPELINES();

    // Convert the coarse count into a coarse partitioning.
    sum = 0;
    for( subsort = 0; subsort < n_subsort; subsort++ )
    {
      for( pipeline_rank = 0; pipeline_rank < n_pipeline; pipeline_rank++ )
      {
        i                   = subsort + cp_stride * pipeline_rank;
        count               = coarse_partition[i];
        coarse_partition[i] = sum;
        sum                += count;
      }
    }

    // Do the coarse sort.
    EXEC_PIPELINES( coarse_sort, args, 0 );

    WAIT_PIPELINES();

    // Convert the coarse partitioning used during the coarse sort into the
    // partitioning of the particle list by subsort pipelines.
    coarse_partition[ n_subsort ] = n_particle;

    // Do fine grained subsorts.  While the fine grained subsorts are
    // executing, clear the ghost parts of the partitioning array.
    EXEC_PIPELINES( subsort, args, 0 );

    CLEAR( partition, vl );

    for( i = vh + 1; i < n_voxel; i++ )
    {
      partition[i] = n_particle;
    }

    WAIT_PIPELINES();

  }

  else
  {
    // Just do the subsort when single threaded.  We need to hack the aux
    // arrays and what not to make it look like coarse sorting was done to
    // the subsort pipeline.
    coarse_partition[0] = 0;
    coarse_partition[1] = n_particle;

    args->p     = aux_p;
    args->aux_p = p;

    subsort_pipeline( args, 0, 1 );

    CLEAR( partition, vl );

    for( i = vh + 1; i < n_voxel; i++ )
    {
      partition[i] = n_particle;
    }

    // Results ended up in the wrong place as a result of the ugly hack above.
    // Copy it to the right place and undo the above hack. FIXME: IF WILLING
    // TO MOVE SP->P AROUND AND DO MORE MALLOCS PER STEP I.E. HEAP
    // FRAGMENTATION, COULD AVOID THIS COPY.
    COPY( p, aux_p, n_particle );
  }
}

#endif
