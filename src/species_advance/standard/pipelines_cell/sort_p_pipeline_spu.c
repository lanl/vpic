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
#define P2V( p, P, V ) ( ((int64_t)(p)*(int64_t)(V) + (int64_t)((P)-1)) / \
                         (int64_t)(P) )

#define N_BLOCK 512

void
_SPUEAR_coarse_count_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                   int pipeline_rank,
                                   int n_pipeline ) {
  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  mfc_write_tag_mask( 1 );
  mfc_read_tag_status_all();

  /*const*/ double n_target = (double)(args->np/4) / (double)n_pipeline;
  /**/      int i  = 4*(int)( n_target*(double) pipeline_rank    + 0.5 );
  /*const*/ int i1 = 4*(int)( n_target*(double)(pipeline_rank+1) + 0.5 );
  /*const*/ int n_voxel = args->n_voxel;
  MEM_PTR( const particle_t, 128 ) p_src = args->p;
  int p;
  int * __restrict count = args->coarse_partition + pipeline_rank*n_pipeline;

  int b, n, nb, n_block[2];
  const particle_t * __restrict ALIGNED(128) p_block;
  DECLARE_ALIGNED_ARRAY( particle_t, 128, p_buffer, 2*N_BLOCK );

  /* Clear local count */
  for( p=0; p<n_pipeline; p++ ) count[p] = 0;

  /* Do local coarse count */
  /* FIXME: THIS IS 8X THEORETICALLY BANDWIDTH INEFFICIENT.  BUT,
     UNKNOWN IF PRACTICALLY IT COULD BE IMPROVED MUCH WITH FINE
     GRAINED DMA TRANSFERS. */

  /* Particles are assigned to pipelines in blocks of 4.  The last
     pipeline handles stragglers. */
  if( pipeline_rank==n_pipeline-1 ) i1 = args->np;

# define READ_PBLOCK( b ) do {                                          \
    nb = i1 - i;                                                        \
    if( nb>N_BLOCK ) nb = N_BLOCK;                                      \
    if( nb ) mfc_get( p_buffer + N_BLOCK*(b),                           \
                      p_src + i*sizeof(particle_t),                     \
                      nb*sizeof(particle_t),                            \
                      (b), 0, 0 );                                      \
    i += nb;                                                            \
    n_block[b] = nb;                                                    \
  } while(0)
  
# define PROCESS_PBLOCK(b) do {                                         \
    nb = n_block[b];                                                    \
    if( nb ) {                                                          \
      mfc_write_tag_mask( 1<<b );                                       \
      mfc_read_tag_status_all();                                        \
      p_block = p_buffer + N_BLOCK*(b);                                 \
      for( n=0; n<nb; n++ )                                             \
        count[ V2P( p_block[n].i, n_pipeline, n_voxel ) ]++;            \
    }                                                                   \
  } while(0)

  READ_PBLOCK( 0 );
  for( b=0; n_block[b]; b=!b ) {
    READ_PBLOCK( !b );
    PROCESS_PBLOCK( b );
  }

  /* Copy local coarse counts to the output */
  argp += ( (uint64_t)args->coarse_partition - (uint64_t)args +
            sizeof(int)*pipeline_rank*n_pipeline );
  for( p=0; p<n_pipeline; p++) {
    mfc_write_tag_mask( 1 << (p&31) );
    mfc_read_tag_status_all();
    mfc_put( &count[p], argp + sizeof(int)*p, sizeof(int), p&31, 0, 0 );
  }

  /* Make sure all DMA transfers have completed before returning */
  mfc_write_task_mask( 0xffffffff );
  mfc_read_tag_status_all();
}

void
_SPUEAR_coarse_sort_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                  int pipeline_rank,
                                  int n_pipeline ) {
  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  mfc_write_tag_mask( 1 );
  mfc_read_tag_status_all();

  /*const*/ double n_target = (double)(args->np/4) / (double)n_pipeline;
  /**/      int i  = 4*(int)( n_target*(double) pipeline_rank    + 0.5 );
  /*const*/ int i1 = 4*(int)( n_target*(double)(pipeline_rank+1) + 0.5 );
  /*const*/ int n_voxel = args->n_voxel;
  MEM_PTR( const particle_t, 128 ) p_src = args->p;
  MEM_PTR(       particle_t, 128 ) p_dst = args->aux_p;
  int j, p;
  int * __restrict next;

  int channel, b, n, nb, n_block[2];
  const particle_t * __restrict ALIGNED(128) p_block;
  DECLARE_ALIGNED_ARRAY( particle_t, 128, p_buffer, 2*N_BLOCK );

  /* Copy local partitioning to next array */
  next = args->coarse_partition + pipeline_rank*n_pipeline;

  /* Copy particles into aux array in coarse sorted order */

  /* Particles are assigned to pipelines in blocks of 4.  The last
     pipeline handles stragglers. */
  if( pipeline_rank==n_pipeline-1 ) i1 = args->np;

# define READ_PBLOCK( b ) do {                                          \
    nb = i1 - i;                                                        \
    if( nb>N_BLOCK ) nb = N_BLOCK;                                      \
    if( nb ) mfc_get( p_buffer + N_BLOCK*(b),                           \
                      p_src + i*sizeof(particle_t),                     \
                      nb*sizeof(particle_t),                            \
                      (b), 0, 0 );                                      \
    i += nb;                                                            \
    n_block[b] = nb;                                                    \
  } while(0)
  
# define PROCESS_PBLOCK(b) do {                                         \
    nb = n_block[b];                                                    \
    if( nb ) {                                                          \
      mfc_write_tag_mask( 1<<b );                                       \
      mfc_read_tag_status_all();                                        \
      p_block = p_buffer + N_BLOCK*(b);                                 \
      for( n=0; n<nb; n++ ) {                                           \
        p = V2P( p_block[n].i, n_pipeline, n_voxel );                   \
        j = next[p]++;                                                  \
        /* FIXME: IS THE WAIT NECESSARY? */                             \
        mfc_write_task_mask( 1<<channel );                              \
        mfc_read_tag_status_all();                                      \
        mfc_put( &p_block[n], p_dst + j*sizeof(particle_t),             \
                 sizeof(particle_t), channel, 0, 0 );                   \
        channel++; if( channel>31 ) channel = 2;                        \
      }                                                                 \
    }                                                                   \
  } while(0)

  channel = 2;
  READ_PBLOCK( 0 );
  for( b=0; n_block[b]; b=!b ) {
    READ_PBLOCK( !b );
    PROCESS_PBLOCK( b );
  }

  /* Make sure all DMA transfers have completed before returning */
  mfc_write_task_mask( 0xffffffff );
  mfc_read_tag_status_all();
}

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
