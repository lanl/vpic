#define IN_spa
#define HAS_SPU_PIPELINE
#include "../spa_private.h"

#include <spu_mfcio.h>

// See ../sort_p.c for details.
#define V2B( v, B, V ) ( ((int64_t)(v)*(int64_t)(B)) / (int64_t)(V) )

#define BS coarse_block_size
#define MB max_coarse_bucket

void
_SPUEAR_coarse_count_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                   int pipeline_rank,
                                   int n_pipeline ) {
  MEM_PTR( const particle_t, 32 ) p_src;
  int i, i1, n_voxel, n_coarse_bucket;
  int b;

  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  DECLARE_ALIGNED_ARRAY( int, 128, voxel, BS );
  DECLARE_ALIGNED_ARRAY( int, 128, count, MB );

  // Get pipeline args and clear local coarse counts
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  CLEAR( count, n_coarse_bucket );
  mfc_write_tag_mask( 1<<0 );
  mfc_read_tag_status_all();
  p_src           = args->p;
  DISTRIBUTE( args->n, BS, pipeline_rank, n_pipeline, i, i1 ); i1 += i;
  n_voxel         = args->n_voxel;
  n_coarse_bucket = args->n_coarse_bucket;

  // Begin getting the voxels for the initial particle block
  if( i<i1 ) 
    for( b=0; b<BS; b++ )
      mfc_get( voxel+b, p_src+sizeof(particle_t)*(i+b)+3*sizeof(float),
               sizeof(int32_t), b, 0, 0 );

  // For all particle blocks
  for( ; i<i1; i+=BS ) {

    // For each particle in this block, wait for the particle's voxel
    // data to arrive, update the local count and starting getting the
    // corresponding particle voxel data in the next block
    for( b=0; b<BS; b++ ) {
      mfc_write_tag_mask( 1 << b );
      mfc_reag_tag_status_all();
      count[ V2B( voxel[b], n_coarse_bucket, n_voxel ) ]++;
      if( i+BS<i1 )
        mfc_get( voxel+b, p_src+sizeof(particle_t)*(i+BS+b)+3*sizeof(float),
                 sizeof(int32_t), b, 0, 0 );
    }
  }

  // Copy local count to the output and wait for all memory
  // transactions to complete
  mfc_put( count, args->coarse_partition + sizeof(int)*MB*pipeline_rank,
           sizeof(int)*MB, 31, 0, 0 );
  mfc_write_tag_mask( 0xffffffff );
  mfc_read_tag_status_all();
}

void
_SPUEAR_coarse_sort_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                  int pipeline_rank,
                                  int n_pipeline ) {
  MEM_PTR( const particle_t, 32 ) p_src;
  MEM_PTR(       particle_t, 32 ) p_dst;
  int i, i1, n_voxel, n_coarse_bucket;
  int b;

  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  DECLARE_ALIGNED_ARRAY( particle_t, 128, p_in, BS );
  DECLARE_ALIGNED_ARRAY( int, 128, next, MB );

  // Get pipeline args
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  mfc_write_tag_mask( 1<<0 );
  mfc_read_tag_status_all();
  p_src           = args->p;
  p_dst           = args->aux_p;
  DISTRIBUTE( args->n, BS, pipeline_rank, n_pipeline, i, i1 ); i1 += i;
  n_voxel         = args->n_voxel;
  n_coarse_bucket = args->n_coarse_bucket;

  // Get the coarse partitioning
  mfc_get( next, args->coarse_partition + sizeof(int)*MB*pipeline_rank,
           sizeof(int)*MB, 0, 0, 0 );
  mfc_write_tag_mask( 1<<0 );
  mfc_read_tag_status_all();
  
  // Begin getting the voxels for the initial particle block
  if( i<i1 ) 
    for( b=0; b<BS; b++ )
      mfc_get( p_in+b, p_src+sizeof(particle_t)*(i+b), sizeof(particle_t),
               b, 0, 0 );

  // For all particle blocks
  for( ; i<i1; i+=BS ) {

    // For each particle in this block, wait for the particle's load
    // to complete, wait for the store of the corresponding particle
    // in the previous block to compelete, determine where to store
    // the arrived particle, copy it into a local writeback buffer,
    // begin the load of the corresponding particle in the next block
    // and the store of the particle in this block.

    for( b=0; b<BS; b++ ) {
      mfc_write_tag_mask( 1<<b | 1<<(b+BS) );
      mfc_reag_tag_status_all();     
      j = next[ V2B( p_in[b].i, n_coarse_bucket, n_voxel ) ]++;     
      p_out[b] = p_in[b]; // FIXME: SPU accelerate?
      mfc_put( p_out+b, p_dst+sizeof(particle_t)*j, sizeof(particle_t),
               b+BS, 0, 0 );
      if( i+BS<i1 )
        mfc_get( p_in+b, p_src+sizeof(particle_t)*(i+BS+b), sizeof(particle_t),
                 b, 0, 0 );
    }
  }

  // Wait for all memory transactions to complete
  mfc_write_tag_mask( 0xffffffff );
  mfc_read_tag_status_all();
}
























#if 0

#define CLEAR( n ) do {                                                 \
    mfc_write_tag_mask( (1<<channel) );                                 \
    mfc_read_tag_status_all();                                          \
    mfc_put( zeros, d + i*sizeof(int), n*sizeof(int), channel, 0, 0 );  \
    channel++;                                                          \
    if( channel>31 ) channel = 0;                                       \
    i += n;                                                             \
  } while(0)
#define MAX_LS_ZEROS 4096 /* 16K chunks, maximum */

static void
clear_range( MEM_PTR( int, 128 ) d,
             int i,
             int i1 ) {
  int n, n_bulk, channel = 0;
  DECLARE_ALIGNED_ARRAY( int, 128, zeros, MAX_LS_ZEROS );

  /* Create a buffer of zeros on the local store for clearing the memory */
  memset( zeros, 0,
          sizeof(int)*( (i1-i)<MAX_LS_ZEROS ? (i1-i) : MAX_LS_ZEROS ) );

  /* Clear until i is aligned at 16 bytes */
  while( i&3 && i<i1 ) CLEAR(1);
 
  /* Clear the bulk as fast as possible */
  n_bulk = (i1-i)&(~3); /* Number of integers we can clear in 16 byte
                           aligned chunks */
  while( n_bulk ) {
    n = n_bulk;
    if( n>N_LS_ZEROS ) n = N_LS_ZEROS;
    CLEAR(n);
    n_bulk -= n;
  }

  /* Clear non aligned stragglers */
  while( i<i1 ) CLEAR(1);

  /* Make sure all DMA transfers have completed before returning */
  mfc_write_task_mask( 0xffffffff );
  mfc_read_tag_status_all();
}

#undef MAX_LS_ZEROS
#undef CLEAR

// sizeof(what) better be 1,2,4,8 or a multiple of 16.
// FIXME: WILL THIS RUN AFOUL OF ALIGNMENT RESTRICTIONS?

#define LOAD( what, ls, mem, index )                                    \
  mfc_get( (ls), (mem) + (index)*sizeof(what), sizeof(what), 0, 0 );    \
  mfc_write_tag_mask( 1 );                                              \
  mfc_read_tag_status_all()

#define STORE( what, ls, pmem, index )                                  \
  mfc_put( (ls), (mem) + (index)*sizeof(what), sizeof(what), 0, 0 );    \
  mfc_write_tag_mask( 1 );                                              \
  mfc_read_tag_status_all()

void
_SPUEAR_subsort_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                              int pipeline_rank,
                              int n_pipeline ) {
  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  mfc_write_tag_mask( 1 );
  mfc_read_tag_status_all();

  /* This pipeline is to sort particles in voxels [v0,v1). */
  /*const*/ int n_voxel = args->n_voxel;
  /*const*/ int v0 = P2V( pipeline_rank,   n_pipeline, n_voxel );
  /*const*/ int v1 = P2V( pipeline_rank+1, n_pipeline, n_voxel );
  
  /* Particles in this voxel range in [i0,i1) in the aux array */
  /*const*/ int i0 = args->coarse_partition[pipeline_rank  ];
  /*const*/ int i1 = args->coarse_partition[pipeline_rank+1];
  MEM_PTR( const particle_t, 128 ) p_src = args->aux_p;
  MEM_PTR(       particle_t, 128 ) p_dst = args->p;

  MEM_PTR( int, 128 ) partition = args->partition;
  MEM_PTR( int, 128 ) next      = args->next;
  
  int i, j, v, sum, count;

  particle_t p_i; /* FIXME: ALIGNMENT? */
  int next_v; /* FIXME: ALIGNMENT? */
  int partition_v;

  /* Clear fine grained count */
  clear_range( next, v0, v1 );

  /* Fine grained count
     FIXME: THIS DESPARATELY NEEDS CACHING AND BUFFERING */
  for( i=i0; i<i1; i++ ) {
    LOAD( particle_t, &p_i, p_src, i ); v = p_i.i;
    LOAD( int, &next_v, next, v ); next_v++; STORE( int, &next_v, next, v );
  }

  /* Convert the fine grained count into a partitioning
     FIXME: THIS DESPARATELY NEEDS BUFFERING */
  sum = i0;
  for( v=v0; v<v1; v++ ) {
    LOAD( int, &next_v, next, v ); count = next_v;
    next_v      = sum; STORE( int, &next_v,      next,      v );
    partition_v = sum; STORE( int, &partition_v, partition, v );
    sum += count;
  }
  partition_v = sum; STORE( int, &partition_v, partition, v );

  /* Local fine grained sort */
  /* FIXME: THIS DESPARATELY NEEDS CACHING AND BUFFERING */
  for( i=i0; i<i1; i++ ) {
    LOAD( particle_t, &p_i, p_src, i ); v = p_i.i;
    LOAD( int, &next_v, next, v ); j=next_v++; STORE( int, &next_v, next, v );
    STORE( particle_t, &p_i, p_dst, j );
  }
}

#undef STORE
#undef LOAD

#endif
