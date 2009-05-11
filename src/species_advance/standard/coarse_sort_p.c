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

#define BS sort_block_size
#define MP max_subsort_pipeline

void
coarse_count_pipeline( sort_p_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {
  const particle_t * RESTRICT ALIGNED(128) p_src;
  int i, i1, n_subsort_pipeline, vl, vh;

  int count[ MP ]; // On pipe stack to avoid cache hot spots

  // Get pipeline args
  p_src              = args->p;
  DISTRIBUTE( args->n, BS, pipeline_rank, n_pipeline, i, i1 ); i1 += i;
  n_subsort_pipeline = args->n_subsort_pipeline;
  vl                 = args->vl;
  vh                 = args->vh;
  
  // Clear local coarse count
  CLEAR( count, n_subsort_pipeline );
  
  // Local coarse count input particles
  for( ; i<i1; i++ ) count[ V2P( p_src[i].i, n_subsort_pipeline, vl, vh ) ]++;
  
  // Copy local coarse count to output
  COPY( args->coarse_partition + MP*pipeline_rank, count, n_subsort_pipeline );
}

void
coarse_sort_pipeline( sort_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline ) {
  const particle_t * RESTRICT ALIGNED(128) p_src;
  /**/  particle_t * RESTRICT ALIGNED(128) p_dst;
  int i, i1, n_subsort_pipeline, vl, vh;
  int j;

  int next[ MP ]; // On pipeline stack to avoid cache hot spots and to
                  // allow reuse of coarse partitioning for fine sort
                  // stage.

  // Get pipeline args
  p_src              = args->p;
  p_dst              = args->aux_p;
  DISTRIBUTE( args->n, BS, pipeline_rank, n_pipeline, i, i1 ); i1 += i;
  n_subsort_pipeline = args->n_subsort_pipeline;
  vl                 = args->vl;
  vh                 = args->vh;
  
  // Load local coarse partitioning to next
  COPY( next, args->coarse_partition + MP*pipeline_rank, n_subsort_pipeline );

  // Copy particles into aux array in coarse sorted order
  for( ; i<i1; i++ ) {
    j = next[ V2P( p_src[i].i, n_subsort_pipeline, vl, vh ) ]++;
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
  int b, n_subsort_pipeline = args->n_subsort_pipeline;
  int p, n_pipeline = N_PIPELINE + 1; // Include straggler cleanup pipeline
  int count, sum = 0;

  // Do the coarse count
  EXEC_PIPELINES( coarse_count, args, 0 );
  WAIT_PIPELINES();
  
  // Convert the coarse count into a coarse partitioning
  for( b=0; b<n_subsort_pipeline; b++ )
    for( p=0; p<n_pipeline; p++ ) {
      count = coarse_partition[ b + MP*p ];
      coarse_partition[ b + MP*p ] = sum;
      sum += count;
    }
  
  // Do the coarse sort
  EXEC_PIPELINES( coarse_sort, args, 0 );
  WAIT_PIPELINES();  
}
