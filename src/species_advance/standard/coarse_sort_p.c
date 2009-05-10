#define IN_spa
#define HAS_SPU_PIPELINE
#include "spa_private.h"

// FIXME: ALTIVEC ACCELERATE!
#if defined(__SSE__)
#include <xmmintrin.h>
#endif

// FIXME: ADD RESTRICT TO UTIL_BASE.H
#ifndef RESTRICT
#define RESTRICT __restrict
#endif

// See sort_p.c for details.
#define V2B( v, B, V ) ( ((int64_t)(v)*(int64_t)(B)) / (int64_t)(V) )

#define BS coarse_block_size
#define MB max_coarse_bucket

void
coarse_count_pipeline( sort_p_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {
  const particle_t * RESTRICT ALIGNED(128) p_src;
  int i, i1, n_voxel, n_coarse_bucket;

  int count[ MB ]; // On pipe stack to avoid cache hot spots

  // Get pipeline args
  p_src           = args->p;
  DISTRIBUTE( args->n, BS, pipeline_rank, n_pipeline, i, i1 ); i1 += i;
  n_voxel         = args->n_voxel;
  n_coarse_bucket = args->n_coarse_bucket;
  
  // Clear local coarse count
  CLEAR( count, n_coarse_bucket );
  
  // Local coarse count input particles
  for( ; i<i1; i++ ) count[ V2B( p_src[i].i, n_coarse_bucket, n_voxel ) ]++;
  
  // Copy local coarse count to output
  COPY( args->coarse_partition + MB*pipeline_rank, count, n_coarse_bucket );
}

void
coarse_sort_pipeline( sort_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline ) {
  const particle_t * RESTRICT ALIGNED(128) p_src;
  /**/  particle_t * RESTRICT ALIGNED(128) p_dst;
  int i, i1, n_voxel, n_coarse_bucket;
  int j;

  int next[ MB ]; // On pipeline stack to avoid cache hot spots and to
                  // allow reuse of coarse partitioning for fine sort
                  // stage.

  // Get pipeline args
  p_src           = args->p;
  p_dst           = args->aux_p;
  DISTRIBUTE( args->n, BS, pipeline_rank, n_pipeline, i, i1 ); i1 += i;
  n_voxel         = args->n_voxel;
  n_coarse_bucket = args->n_coarse_bucket;
  
  // Load local coarse partitioning to next
  COPY( next, args->coarse_partition + MB*pipeline_rank, n_coarse_bucket );

  // Copy particles into aux array in coarse sorted order
  for( ; i<i1; i++ ) {
    j = next[ V2B( p_src[i].i, n_coarse_bucket, n_voxel ) ]++;
#   if defined(__SSE__)
    _mm_store_ps( &p_dst[j].dx, _mm_load_ps( &p_src[i].dx ) );
    _mm_store_ps( &p_dst[j].ux, _mm_load_ps( &p_src[i].ux ) );
#   else
    p_dst[j] = p_src[i];
#   endif
  }
}

void
coarse_sort_p( sort_p_pipeline_args_t * ALIGNED(128) args ) { 
  int * ALIGNED(128) coarse_partition = args->coarse_partition;
  int b, n_coarse_bucket = args->n_coarse_bucket;
  int p, n_pipeline = N_PIPELINE + 1; // Include straggler cleanup pipeline
  int count, sum = 0;

  // Do the coarse count
  EXEC_PIPELINES( coarse_count, args, 0 );
  WAIT_PIPELINES();
  
  // Convert the coarse count into a coarse partitioning
  for( b=0; b<n_coarse_bucket; b++ )
    for( p=0; p<n_pipeline; p++ ) {
      count = coarse_partition[ b + MB*p ];
      coarse_partition[ b + MB*p ] = sum;
      sum += count;
    }
  
  // Do the coarse sort
  EXEC_PIPELINES( coarse_sort, args, 0 );
  WAIT_PIPELINES();  
}
