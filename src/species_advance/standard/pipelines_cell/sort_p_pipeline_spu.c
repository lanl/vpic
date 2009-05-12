#define IN_spa
#define HAS_SPU_PIPELINE
#include "../spa_private.h"

#include <spu_mfcio.h>

#define BS sort_block_size
#define MP max_subsort
#define MV max_subsort_voxel

void
_SPUEAR_coarse_count_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                   int pipeline_rank,
                                   int n_pipeline ) {
  MEM_PTR( const particle_t, 32 ) p_src;
  int i, i1, n_subsort, vl, vh;
  int b;

  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  DECLARE_ALIGNED_ARRAY( int, 128, voxel, 4*BS );
  DECLARE_ALIGNED_ARRAY( int, 128, count, MP );

  // Get pipeline args
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  mfc_write_tag_mask( 1<<0 );
  mfc_read_tag_status_all();
  p_src     = args->p;
  DISTRIBUTE( args->n, BS, pipeline_rank, n_pipeline, i, i1 ); i1 += i;
  n_subsort = args->n_subsort;
  vl        = args->vl;
  vh        = args->vh;

  // Begin getting the voxels for the initial particle block
  if( i<i1 )
    for( b=0; b<BS; b++ )
      mfc_get( voxel+4*b+3, p_src+sizeof(particle_t)*(i+b)+3*sizeof(int32_t),
               sizeof(int32_t), b, 0, 0 );

  // Clear the local coarse counts
  CLEAR( count, n_subsort );

  // For all particle blocks
  for( ; i<i1; i+=BS ) {

    // For each particle in this block, wait for the particle's voxel
    // data to arrive, update the local count and starting getting the
    // corresponding particle voxel data in the next block
    for( b=0; b<BS; b++ ) {
      mfc_write_tag_mask( 1<<b );
      mfc_read_tag_status_all();
      count[ V2P( voxel[4*b+3], n_subsort, vl, vh ) ]++;
      if( i+BS<i1 )
        mfc_get( voxel+4*b+3,
                 p_src+sizeof(particle_t)*(i+BS+b)+3*sizeof(float),
                 sizeof(int32_t), b, 0, 0 );
    }
  }

  // Copy local count to the output and wait for all memory
  // transactions to complete
  mfc_put( count, args->coarse_partition + sizeof(int)*MP*pipeline_rank,
           sizeof(int)*MP, 31, 0, 0 );
  mfc_write_tag_mask( 0xffffffff );
  mfc_read_tag_status_all();
}

void
_SPUEAR_coarse_sort_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                                  int pipeline_rank,
                                  int n_pipeline ) {
  MEM_PTR( const particle_t, 32 ) p_src;
  MEM_PTR(       particle_t, 32 ) p_dst;
  int i, i1, n_subsort, vl, vh;
  int b, j;

  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  DECLARE_ALIGNED_ARRAY( particle_t, 128, p_in,  BS );
  DECLARE_ALIGNED_ARRAY( particle_t, 128, p_out, BS );
  DECLARE_ALIGNED_ARRAY( int, 128, next, MP );

  // Get pipeline args
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  mfc_write_tag_mask( 1<<0 );
  mfc_read_tag_status_all();
  p_src     = args->p;
  p_dst     = args->aux_p;
  DISTRIBUTE( args->n, BS, pipeline_rank, n_pipeline, i, i1 ); i1 += i;
  n_subsort = args->n_subsort;
  vl        = args->vl;
  vh        = args->vh;

  // Get the local coarse partitioning
  mfc_get( next, args->coarse_partition + sizeof(int)*MP*pipeline_rank,
           sizeof(int)*MP, 0, 0, 0 );
  mfc_write_tag_mask( 1<<0 );
  mfc_read_tag_status_all();
  
  // Begin getting particles for the initial block
  if( i<i1 ) 
    for( b=0; b<BS; b++ )
      mfc_get( p_in+b, p_src+sizeof(particle_t)*(i+b), sizeof(particle_t),
               b, 0, 0 );
  
  // For all blocks
  for( ; i<i1; i+=BS ) {

    // For each particle in this block, wait for the particle's load
    // to complete, wait for the store of the corresponding particle
    // in the previous block to complete, determine where to store the
    // arrived particle, copy it into a local writeback buffer, begin
    // the load of the corresponding particle in the next block and
    // the store of the particle in this block.

    for( b=0; b<BS; b++ ) {
      mfc_write_tag_mask( (1<<b) | (1<<(b+BS)) );
      mfc_read_tag_status_all();     
      j = next[ V2P( p_in[b].i, n_subsort, vl, vh ) ]++;     
      ((vec_float4 *)(p_out+b))[0] = ((vec_float4 *)(p_in+b))[0];
      ((vec_float4 *)(p_out+b))[1] = ((vec_float4 *)(p_in+b))[1];
      //p_out[b] = p_in[b]; // FIXME: SPU accelerate?
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

void
_SPUEAR_subsort_pipeline_spu( MEM_PTR( sort_p_pipeline_args_t, 128 ) argp,
                              int pipeline_rank,
                              int n_pipeline ) {
  MEM_PTR( const particle_t, 128 ) p_src;
  MEM_PTR(       particle_t, 128 ) p_dst;
  int i0, i1, v0, v1, i, j, v, sum, count;
  int b, nb;

  int * RESTRICT partition;

  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args,      1    );
  DECLARE_ALIGNED_ARRAY( int,                    128, _partition, MV+1+3 );
  DECLARE_ALIGNED_ARRAY( int,                    128, next,      MV   );
  DECLARE_ALIGNED_ARRAY( particle_t,             128, p_in,      BS   );
  DECLARE_ALIGNED_ARRAY( particle_t,             128, p_out,     BS   );

  // Get pipeline args
  mfc_get( args, argp, sizeof(*args), 0, 0, 0 );
  mfc_write_tag_mask( 1<<0 );
  mfc_read_tag_status_all();

  // This pipeline sorts particles in this voxel range in [i0,i1) in
  // the aux array.  These particles are in voxels [v0,v1).
  p_src = args->aux_p;
  p_dst = args->p;
  mfc_get( next+pipeline_rank,
           args->coarse_partition+sizeof(int)*pipeline_rank,
           sizeof(int), 0, 0, 0 );
  mfc_get( next+pipeline_rank+1,
           args->coarse_partition+sizeof(int)*(pipeline_rank+1),
           sizeof(int), 1, 0, 0 );
  mfc_write_tag_mask( (1<<0)|(1<<1) );
  mfc_read_tag_status_all();          
  i0 = next[pipeline_rank  ];
  i1 = next[pipeline_rank+1];
  v0 = P2V( pipeline_rank,   n_pipeline, args->vl, args->vh );
  v1 = P2V( pipeline_rank+1, n_pipeline, args->vl, args->vh );
  if( pipeline_rank==0            ) v0 = 0;
  if( pipeline_rank==n_pipeline-1 ) v1 = args->n_voxel;
  if( v1-v0 > MV ) return; // FIXME: DIAGNOSTIC ABORT

  // Begin getting the voxels for the initial particle block
  partition = _partition;
  nb = i1-i0; if( nb>BS ) nb = BS;
  for( b=0; b<nb; b++ )
    mfc_get( partition+4*b+3,
             p_src+sizeof(particle_t)*(i0+b)+3*sizeof(int32_t),
             sizeof(int32_t), b, 0, 0 );

  // Clear fine grained count
  CLEAR( next, v1-v0 );

  // Fine grained count
  // For all particle blocks
  for( i=i0; i<i1; i+=BS ) {
    // For each particle in this block, wait for the particle's voxel
    // data to arrive, update the fine grained count and start getting
    // the corresponding particle voxel data in the next block
    nb = i1-i; if( nb>BS ) nb = BS;
    for( b=0; b<nb; b++ ) {
      mfc_write_tag_mask( 1<<b );
      mfc_read_tag_status_all();
      next[ partition[4*b+3]-v0 ]++;
      if( i+BS+b<i1 )
        mfc_get( partition+4*b+3,
                 p_src+sizeof(particle_t)*(i+BS+b)+3*sizeof(int32_t),
                 sizeof(int32_t), b, 0, 0 );
    }
  }

  // Compute the partitioning
  partition += (v0&3); // Align LS partition with EA partition
  sum = i0;
  for( v=v0; v<v1; v++ ) {
    count = next[v-v0];
    next[v-v0] = sum;
    partition[v-v0] = sum;
    sum += count;
  }
  partition[v1-v0] = sum; // All subsorts who write this agree
  b = 0;
  for( v=v0; v<=v1; v+=nb ) {
    nb = v1+1-v; // Number of ints remaining to transfer
    if(      (v &  1) || (nb< 2) ) nb =  1;  //  1-int txfr ( 1-int aligned)
    else if( (v &  3) || (nb< 4) ) nb =  2;  //  2-int txfr ( 2-int aligned)
    else if( (v &  7) || (nb< 8) ) nb =  4;  //  4-int txfr ( 4-int aligned)
    else if( (v & 15) || (nb<16) ) nb =  8;  //  8-int txfr ( 8-int aligned)
    else if( (v & 31) || (nb<32) ) nb = 16;  // 16-int txfr (16-int aligned)
    else { // Large txfr (32-int aligned, 4-int multiple, under 16KB)
      nb &= ~3;
      if( nb>4096 ) nb = 4096;
    }
    mfc_put( partition+v-v0, args->partition+sizeof(int)*v, sizeof(int)*nb,
             b, 0, 0 ); b = (b+1) & 31;
  }

  // Begin getting particles for the initial particle block
  nb = i1-i0; if( nb>BS ) nb = BS;
  for( b=0; b<nb; b++ )
    mfc_get( p_in+b, p_src+sizeof(particle_t)*(i0+b), sizeof(particle_t),
             b, 0, 0 );

  // Fine grained sort
  // For all particle blocks
  for( i=i0; i<i1; i+=BS ) {
    // For each particle in this block, wait for the particle's load
    // to complete, wait for the store of the corresponding particle
    // in the previous block to complete, determine where to store the
    // arrived particle, copy it into a local writeback buffer, begin
    // the load of the corresponding particle in the next block and
    // the store of the particle in this block.
    nb = i1-i; if( nb>BS ) nb = BS;
    for( b=0; b<nb; b++ ) {
      mfc_write_tag_mask( (1<<b) | (1<<(b+BS)) );
      mfc_read_tag_status_all();     
      j = next[ p_in[b].i-v0 ]++;     
      ((vec_float4 *)(p_out+b))[0] = ((vec_float4 *)(p_in+b))[0];
      ((vec_float4 *)(p_out+b))[1] = ((vec_float4 *)(p_in+b))[1];
      mfc_put( p_out+b, p_dst+sizeof(particle_t)*j, sizeof(particle_t),
               b+BS, 0, 0 );
      if( i+BS+b<i1 )
        mfc_get( p_in+b, p_src+sizeof(particle_t)*(i+BS+b), sizeof(particle_t),
                 b, 0, 0 );
    }
  }

  // Complete all memory transactions
  mfc_write_tag_mask( 0xffffffff );
  mfc_read_tag_status_all();
}
