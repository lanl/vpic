#define IN_spa
#define HAS_SPU_PIPELINE
#include "../spa_private.h"

#if 0

#include <spu_mfcio.h>

#ifdef IN_HARNESS
#include <profile.h>
#endif

/* See ../sort_p.c for details. */
#define V2P( v, P, V ) ( ((int64_t)(v)*(int64_t)(P)) / (int64_t)(V) )
#define P2V( p, P, V ) ( ((int64_t)(p)*(int64_t)(V) + (int64_t)((P)-1)) / (int64_t)(P) )

void
_SPUEAR_coarse_count_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                   int pipeline_rank,
                                   int n_pipeline ) {
  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  mfc_write_tag_mask( 1 );
  mfc_read_tag_status_all();

  /* Determine our workload */

  /*const*/ double n_target = (double)args->np / (double)n_pipeline;
  /**/      int i  = (int)( n_target*(double) pipeline_rank    + 0.5 );
  /*const*/ int i1 = (int)( n_target*(double)(pipeline_rank+1) + 0.5 );
  /*const*/ int n_voxel = args->n_voxel;
  MEM_PTR( const particle_t, 128 ) p_src = args->p;
  int p;

  /* _coarse_partition mirrors coarse_partition in memory (including
     alignment restrictions for support of small DMA transfers).
     count points to where this pipeline should accumulate counts in
     the mirror (which has the same offset locally as where this
     pipeline should store counts in main memory). */
  DECLARE_ALIGNED_ARRAY( int, 128, _coarse_partition,
                         MAX_PIPELINE*MAX_PIPELINE+1 );
  int * __restrict count =
    _coarse_partition + sizeof(int)*pipeline_rank*n_pipeline;

  int channel, n_channel;

  /* _voxel holds voxel data loaded from particles assigned to this
     pipeline.  It is overallocated to support restrictions on small
     DMA transfers.  voxel is syntactic sugar to access this. */
  DECLARE_ALIGNED_ARRAY( int32_t, 128, _voxel, 128 );
# define voxel(n) _voxel[4*(n)+3]

  /* Clear local coarse count */
  for( p=0; p<n_pipeline; p++ ) count[p] = 0;

  /* Do local coarse count */
  while( i<i1 ) {

    /* Fetch voxel indexes of the next 32 (or so) particles */
    /* FIXME: WOULD IT BE BETTER TO DOUBLE BUFFER THIS? */
    /* FIXME: WOULD IT BE BETTER TO WASTE BANDWIDTH AND USE FEWER DMA
       TRANSFERS? */

    n_channel = i1 - i;
    if( n_channel>32 ) n_channel = 32;
    for( channel=0; channel<n_channel; channel++ ) /* FIXME: Duff-able? */
      mfc_get( &voxel(channel),
               (p_src+sizeof(particle_t)*i+12) + sizeof(particle_t)*channel,
               sizeof(int32_t), channel, 0, 0 ); /* FIXME: Hoist invariant */
    mfc_write_tag_mask( 0xffffffff ); /* FIXME: Need to do this every time? */
    mfc_read_tag_status_all();

    /* Update counts */

    for( channel=0; channel<n_channel; channel++ )
      count[ V2P( voxel(channel), n_pipeline, n_voxel ) ]++;

    /* Move to the next batch of particles */
    i += n_channel;
  }

  /* Copy local coarse count to output */
  /* FIXME: WOULD IT BE BETTER TO WASTE BANDWIDTH AND USE FEWER DMA
     TRANSFERS? */

  for( p=0; p<n_pipeline; p++ ) {
    /* Wait for any pending transactions on this channel to complete
       and then write out the count */
    /* FIXME: IS THIS WAITING ON PENDING NECESSARY? */
    channel = p & 31;
    mfc_write_tag_mask( (1<<channel) );
    mfc_read_tag_status_all();
    mfc_put( &count[p],
             ( args->coarse_partition +
               sizeof(int)*pipeline_rank*n_pipeline ) + sizeof(int)*p,
             sizeof(int), channel, 0, 0 ); /* FIXME: Host invariant */
  }
}

/* FIXME: Assign in blocks of 4 particles */

void
_SPUEAR_coarse_sort_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                  int pipeline_rank,
                                  int n_pipeline ) {
  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  mfc_write_tag_mask( 1 );
  mfc_read_tag_status_all();

  /* Determine our workload */

  /*const*/ double n_target = (double)args->np / (double)n_pipeline;
  /**/      int i  = (int)( n_target*(double) pipeline_rank    + 0.5 );
  /*const*/ int i1 = (int)( n_target*(double)(pipeline_rank+1) + 0.5 );
  /*const*/ int n_voxel = args->n_voxel;
  MEM_PTR( const particle_t, 128 ) p_src = args->p;
  MEM_PTR(       particle_t, 128 ) p_dst = args->aux_p;
  int j, p;

  /* _coarse_partition mirrors coarse_partition in memory (including
     alignment restrictions for support of small DMA transfers).  next
     points to where this pipeline should load initial next values in
     the mirror (which has the same offset locally as where this
     pipeline load initial next values in main memory). */
  DECLARE_ALIGNED_ARRAY( int, 128, _coarse_partition,
                         MAX_PIPELINE*MAX_PIPELINE+1 );
  int * __restrict next =
    _coarse_partition + sizeof(int)*pipeline_rank*n_pipeline;

  int channel, n_channel;

  /* Copy this pipeline's partitioning into next */
  /* FIXME: WOULD IT BE BETTER TO WASTE BANDWIDTH AND USE FEWER DMA
     TRANSFERS? */

  for( p=0; p<n_pipeline; p++ ) {
    /* Wait for any pending transactions on this channel to complete
       and then write out the count */
    /* FIXME: IS THIS WAITING ON PENDING NECESSARY? */
    channel = p & 31;
    mfc_write_tag_mask( (1<<channel) );
    mfc_read_tag_status_all();
    mfc_get( &next[p],
             ( args->coarse_partition +
               sizeof(int)*pipeline_rank*n_pipeline ) + sizeof(int)*p,
             sizeof(int), channel, 0, 0 ); /* FIXME: Host invariant */
  }

  /* Copy particles into aux array in coarse sorted order */

  channel = 2;

# define READ_PBLOCK( b ) do {                                          \
    N = i1 - i;                                                         \
    if( N>512 ) N = 512;                                                \
    if( N ) mfc_get( p_buffer + 512*(b),                                \
                     p_src + i*sizeof(particle_t),                      \
                     N*sizeof(particle_t),                              \
                     (b), 0, 0 );                                       \
    i += N;                                                             \
    n_block[b] = N;                                                     \
  } while(0)
  
# define PROCESS_PBLOCK(b) do {                                         \
    N = n_block[b];                                                     \
    if( N ) {                                                           \
      mfc_write_tag_mask( 1<<b );                                       \
      mfc_read_tag_status_all();                                        \
      p_block = p_buffer + 512*(b);                                     \
      for( n=0; n<N; n++ ) {                                            \
        mfc_write_task_mask( 1<<channel );                              \
        mfc_read_tag_status_all();                                      \
        mfc_put( &p_block[n],                                           \
                 p_dst + sizeof(particle_t)*(next[ V2P( p_block[n].i ) ]++), \
                 sizeof(particle_t),                                    \
                 channel, 0, 0 );                                       \
        channel++; if( channel>31 ) channel = 2;                        \
      }                                                                 \
    }                                                                   \
  } while(0)

  READ_PBLOCK( 0 );
  for( b=0; n_block[b]; b=next[b] ) {
    READ_PBLOCK( next[b] );
    PROCESS_PBLOCK( b );
  }

  /* Make sure all DMA transfers have completed before returning */
  mfc_write_task_mask( 0xffffffff );
  mfc_read_tag_status_all();
}

#else

#include <stdio.h>

void
_SPUEAR_coarse_count_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                   int pipeline_rank,
                                   int n_pipeline ) {
  fprintf( stdout, "In coarse_count_pipeline\n" ); fflush( stdout );
}

void
_SPUEAR_coarse_sort_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                  int pipeline_rank,
                                  int n_pipeline ) {
  fprintf( stdout, "In coarse_sort_pipeline\n" ); fflush( stdout );
}

void
_SPUEAR_subsort_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                              int pipeline_rank,
                              int n_pipeline ) {
  fprintf( stdout, "In subsort_pipeline\n" ); fflush( stdout );
}

#endif
