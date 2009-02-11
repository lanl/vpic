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

#if 0

void
sort_p( species_t * sp,
        const grid_t * g ) {
  particle_t * ALIGNED(128) p = sp->p;
  //const int32_t * __restrict ALIGNED(128) sfc = g->sfc;
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
  //for( i=0; i<np; i++ ) next[ sfc[ p[i].i ] ]++;
  for( i=0; i<np; i++ ) next[ p[i].i ]++;

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
    //for( i=0; i<np; i++ ) out_p[ next[ sfc[ in_p[i].i ] ]++ ] = in_p[i];
    for( i=0; i<np; i++ ) out_p[ next[ in_p[i].i ]++ ] = in_p[i];

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
          //dest = &p[ next[ sfc[ src->i ] ]++ ];
          dest = &p[ next[ src->i ]++ ];
          if( src==dest ) break;
          save_p = *dest;
          *dest  = *src;
          *src   = save_p;
        }
      }

    }
  }
}

#else 

#if defined(__SSE__)
#include <xmmintrin.h>
#endif

/* FIXME: ELIMINATE IN-PLACE / OUT-PLACE OPTIONS FROM SP FIELD */
/* FIXME: ADD RESTRICT TO UTIL_BASE.H */
/* FIXME: Add N_VOXEL convenience field to grid */
#ifndef RESTRICT
#define RESTRICT __restrict
#endif

/* Given the voxel index, compute the pipeline responsible for sorting
   particles within that voxel.  This takes into account that p*V
   might overflow 32-bits.  This macro is robust. */
#define V2P( v, P, V ) ( ((int64_t)(v)*(int64_t)(P)) /  \
                         (int64_t)(V) )

/* Given the pipeline rank, compute the first voxel the pipeline is
   responsible for sorting.  This is based on:
     p = floor(vP/V) =>
     p <= vP/V < p+1 =>
     pV/P <= v < (p+1)V/P
   The range of voxels which satisfy this inequality is then:
     [ ceil(pV/P), ceil((p+1)V/P) )
   In integer math, the lower bound is thus:
     v = (p*V + P-1)/P
   This takes into account that p*V might overflow 32-bits.
   This macro is mostly robust. */
#define P2V( p, P, V ) ( ((int64_t)(p)*(int64_t)(V) + (int64_t)((P)-1)) / \
                         (int64_t)(P) )

void
coarse_count_pipeline( sort_p_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {
  /*const*/ double n_target = (double)args->np / (double)n_pipeline;
  /**/      int i  = (int)( n_target*(double) pipeline_rank    + 0.5 );
  /*const*/ int i1 = (int)( n_target*(double)(pipeline_rank+1) + 0.5 );
  /*const*/ int n_voxel = args->n_voxel;
  const particle_t * RESTRICT ALIGNED(128) p_src = args->p;
  int p;
  int count[ MAX_PIPELINE ]; /* On pipeline stack to avoid cache hot
                                spots */

  if( pipeline_rank==n_pipeline ) return; /* No straggler cleanup needed */
  
  /* Clear local coarse count */
  for( p=0; p<n_pipeline; p++ ) count[p] = 0;
  
  /* Do local coarse count */
  for( ; i<i1; i++ ) count[ V2P( p_src[i].i, n_pipeline, n_voxel ) ]++;
  
  /* Copy local coarse count to output */
  for( p=0; p<n_pipeline; p++ )
    args->coarse_partition[ pipeline_rank*n_pipeline + p ] = count[p];
}

void
coarse_sort_pipeline( sort_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline ) {
  /*const*/ double n_target = (double)args->np / (double)n_pipeline;
  /**/      int i  = (int)( n_target*(double) pipeline_rank    + 0.5 );
  /*const*/ int i1 = (int)( n_target*(double)(pipeline_rank+1) + 0.5 );
  /*const*/ int n_voxel = args->n_voxel;
  const particle_t * RESTRICT ALIGNED(128) p_src = args->p;
  /**/  particle_t * RESTRICT ALIGNED(128) p_dst = args->aux_p;
  int j, p;
  int next[ MAX_PIPELINE ]; /* On pipeline stack to avoid cache hot
                               spots and to allow reuse of coarse
                               partition for fine sort stage. */

  if( pipeline_rank==n_pipeline ) return; /* No straggler cleanup needed */

  /* Copy local partitioning to next array */
  for( p=0; p<n_pipeline; p++ )
    next[p] = args->coarse_partition[ pipeline_rank*n_pipeline + p ];
  
  /* Copy particles into aux array in coarse sorted order */
  for( ; i<i1; i++ ) {
    p = V2P( p_src[i].i, n_pipeline, n_voxel );
    j = next[p]++;
#   if defined(__SSE__)
    _mm_store_ps( &p_dst[j].dx, _mm_load_ps( &p_src[i].dx ) );
    _mm_store_ps( &p_dst[j].ux, _mm_load_ps( &p_src[i].ux ) );
#   else
    p_dst[j] = p_src[i];
#   endif
  }
}

void
subsort_pipeline( sort_p_pipeline_args_t * args,
                  int pipeline_rank,
                  int n_pipeline ) {
  /* This pipeline is to sort particles in voxels [v0,v1). */
  /*const*/ int n_voxel = args->n_voxel;
  /*const*/ int v0 = P2V( pipeline_rank,   n_pipeline, n_voxel );
  /*const*/ int v1 = P2V( pipeline_rank+1, n_pipeline, n_voxel );
  
  /* Particles in this voxel range in [i0,i1) in the aux array */
  /*const*/ int i0 = args->coarse_partition[pipeline_rank  ];
  /*const*/ int i1 = args->coarse_partition[pipeline_rank+1];    
  const particle_t * RESTRICT ALIGNED(128) p_src = args->aux_p;
  /**/  particle_t * RESTRICT ALIGNED(128) p_dst = args->p;
  
  int * RESTRICT ALIGNED(128) partition = args->partition;
  int * RESTRICT ALIGNED(128) next      = args->next;
  
  int i, j, v, sum, count;

  if( pipeline_rank==n_pipeline ) return; /* No straggler cleanup needed */

  /* Clear fine grained count */
  CLEAR( &next[v0], v1-v0 );

  /* Fine grained count */
  for( i=i0; i<i1; i++ ) next[ p_src[i].i ]++;

  /* Compute the partitioning */
  sum = i0;
  for( v=v0; v<v1; v++ ) {
    count = next[v];
    next[v] = sum;
    partition[v] = sum;
    sum += count;
  }
  partition[v1] = sum; /* All threads who write this agree */

  /* Local fine grained sort */
  for( i=i0; i<i1; i++ ) {
    v = p_src[i].i;
    j = next[v]++;
#   if defined(__SSE__)
    _mm_store_ps( &p_dst[j].dx, _mm_load_ps( &p_src[i].dx ) );
    _mm_store_ps( &p_dst[j].ux, _mm_load_ps( &p_src[i].ux ) );
#   else
    p_dst[j] = p_src[i];
#   endif
  }
}

void
sort_p( species_t * sp,
        const grid_t * g ) {

  static char * ALIGNED(128) scratch = NULL;
  static size_t max_scratch = 0;
  size_t sz_scratch;

  int n_pipeline = N_PIPELINE;
  int n_voxel    = (g->nx+2)*(g->ny+2)*(g->nz+2);
  int p, q, sum, count;

  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );

  /* FIXME: TEMPORARY HACK UNTIL SPECIES_ADVANCE API INSTALLED */
  if( sp->partition==NULL ) MALLOC_ALIGNED( sp->partition, n_voxel+1, 128 );

  /* Insure enough scratch space is allocated for the sorting */
  sz_scratch = ( sizeof(args->aux_p[0])*sp->np +            /* For aux_p */
                 128 + sizeof(args->next[0])*(n_voxel+1) ); /* For next */
  if( sz_scratch > max_scratch ) {
    FREE_ALIGNED( scratch );
    MALLOC_ALIGNED( scratch, sz_scratch, 128 );
    max_scratch = sz_scratch;
  }

  /* Setup pipeline arguments */
  args->p         = sp->p;
  args->aux_p     = (particle_t *)scratch;
  args->partition = sp->partition;
  args->next      = (int *)( args->aux_p + sp->np );//ALIGN_PTR( int, args->aux_p + sp->np, 128 );
  args->np        = sp->np;
  args->n_voxel   = n_voxel;

  if( n_pipeline>1 ) {

    /* Thread parallel coarse count */
    EXEC_PIPELINES( coarse_count, args, 0 );
    WAIT_PIPELINES();

    /* Convert the coarse count into a coarse partitioning */
    sum = 0;
    for( p=0; p<n_pipeline; p++ )
      for( q=0; q<n_pipeline; q++ ) {
        count = args->coarse_partition[ p + q*n_pipeline ];
        args->coarse_partition[ p+ q*n_pipeline ] = sum;
        sum += count;
      }
    
    /* Copy input into aux array in coarse sorted order */
    EXEC_PIPELINES( coarse_sort, args, 0 );
    WAIT_PIPELINES();
    
    /* Do fine grained subsorts */
    args->coarse_partition[n_pipeline] = sum; /* Convert coarse partition into
                                           the particle ranges assigned
                                           to each pipeline */
    EXEC_PIPELINES( subsort, args, 0 );
    WAIT_PIPELINES();

  } else {

    /* Just do the subsort when single threaded.  We need to hack the
       aux arrays and what not to make it look like coarse sorting was
       done to the subsort pipeline. */
    args->p                            = (particle_t *)scratch;
    args->aux_p                        = sp->p;
    args->coarse_partition[0]          = 0;
    args->coarse_partition[n_pipeline] = sp->np;
    subsort_pipeline( args, 0, 1 );

    /* Results ended up in the wrong place as a result of the ugly
       hack above.  Copy it to the right place and undo the above
       hack.  FIXME: IF WILLING TO MOVE SP->P AROUND AND DO MORE
       MALLOCS PER STEP (I.E. HEAP FRAGMENTATION), COULD AVOID THIS
       COPY. */
    COPY( args->aux_p, args->p, args->np );

  }
}

#endif
