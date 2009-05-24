#define IN_spa
#define HAS_SPU_PIPELINE
#include "../spa_private.h"

#include <spu_mfcio.h>

// Lookhead must be a power-of-two.  8 saturates MFC queues, 16 saturates
// channels

#define LOOKAHEAD 8
#define MP max_subsort

void
_SPUEAR_coarse_count_pipeline_spu( sort_p_pipeline_args_t * args,
                                   int pipeline_rank,
                                   int n_pipeline ) {
  MEM_PTR( const particle_t, 32 ) p_src = args->p;
  int i, i1, n_subsort = args->n_subsort, vl = args->vl, vh = args->vh, b;

  int        * RESTRICT ALIGNED(128) count; SPU_MALLOC(count, MP,        128);
  particle_t * RESTRICT ALIGNED(128) p_in;  SPU_MALLOC(p_in,  LOOKAHEAD, 128);

  DISTRIBUTE( args->n, 1, pipeline_rank, n_pipeline, i, i1 ); i1 += i;

  // Begin getting the voxels for the initial particle block
  for( b=0; b<LOOKAHEAD; b++ )
    if( LIKELY( i+b<i1 ) )
      mfc_get( &p_in[b].i, p_src+sizeof(particle_t)*(i+b)+3*sizeof(int32_t),
               sizeof(int32_t), b, 0, 0 );

  // Clear the local coarse count
  CLEAR( count, n_subsort );

  // Local coarse count the input particles
  b = 0;
  for( ; i<i1; i++ ) {

    // Wait for the load of this particle's voxel index to complete,
    // update the local count and begin the load of the corresponding
    // particle in the next block's voxel index index
    mfc_write_tag_mask( 1<<b ); mfc_read_tag_status_all();
    count[ V2P( p_in[b].i, n_subsort, vl, vh ) ]++;
    if( LIKELY( i+LOOKAHEAD<i1 ) )
      mfc_get( &p_in[b].i,
               p_src+sizeof(particle_t)*(i+LOOKAHEAD)+3*sizeof(int32_t),
               sizeof(int32_t), b, 0, 0 );
    b = (b+1) & (LOOKAHEAD-1);
  }

  // Copy local coarse count to output
  mfc_put( count, args->coarse_partition + sizeof(int)*MP*pipeline_rank,
           sizeof(int)*MP, 31, 0, 0 );
}

void
_SPUEAR_coarse_sort_pipeline_spu( sort_p_pipeline_args_t * args,
                                  int pipeline_rank,
                                  int n_pipeline ) {
  MEM_PTR( const particle_t, 32 ) p_src = args->p;
  MEM_PTR(       particle_t, 32 ) p_dst = args->aux_p;
  int i, i1, n_subsort = args->n_subsort, vl = args->vl, vh = args->vh;
  int b, j;

  int        * RESTRICT ALIGNED(128) next;  SPU_MALLOC(next,  MP,        128);
  particle_t * RESTRICT ALIGNED(128) p_in;  SPU_MALLOC(p_in,  LOOKAHEAD, 128);
  particle_t * RESTRICT ALIGNED(128) p_out; SPU_MALLOC(p_out, LOOKAHEAD, 128);

  DISTRIBUTE( args->n, 1, pipeline_rank, n_pipeline, i, i1 ); i1 += i;

  // Begin getting particles for the initial block
  for( b=0; b<LOOKAHEAD; b++ )
    if( LIKELY( i+b<i1 ) ) 
      mfc_get( p_in+b, p_src+sizeof(particle_t)*(i+b), sizeof(particle_t),
               b, 0, 0 );
  
  // Load the local coarse partitioning into next
  mfc_get( next, args->coarse_partition + sizeof(int)*MP*pipeline_rank,
           sizeof(int)*MP, LOOKAHEAD, 0, 0 );
  mfc_write_tag_mask( 1<<LOOKAHEAD ); mfc_read_tag_status_all();
  
  // Copy particles into aux array in coarse sorted order 
  b = 0;
  for( ; i<i1; i++ ) {

    // Wait for the load of this particle and the store of the
    // corresponding particle in the previous block to complete,
    // determine where to store the this particle, copy it into a
    // local writeback buffer, begin the load of the corresponding
    // particle in the next block and the store of this particle

    mfc_write_tag_mask( (1<<b)|(1<<(b+LOOKAHEAD)) ); mfc_read_tag_status_all();     
    j = next[ V2P( p_in[b].i, n_subsort, vl, vh ) ]++;     
    ((vec_float4 *)(p_out+b))[0] = ((vec_float4 *)(p_in+b))[0];
    ((vec_float4 *)(p_out+b))[1] = ((vec_float4 *)(p_in+b))[1];
    mfc_put( p_out+b,
             p_dst+sizeof(particle_t)*j, sizeof(particle_t),
             b+LOOKAHEAD, 0, 0 );
    if( LIKELY( i+LOOKAHEAD<i1 ) )
      mfc_get( p_in+b,
               p_src+sizeof(particle_t)*(i+LOOKAHEAD), sizeof(particle_t),
               b, 0, 0 );
    b = (b+1) & (LOOKAHEAD-1);
  }
}

// MV limits local copies of segments of the next and partition
// to using at most 208KB of local store.  This restriction could be
// using more subsorts than SPE and having each SPE process multiple
// subsorts.  But this is not expected to be an issue in
// practice---using 8 SPU, this corresponds to ~459^2 local domains in
// 2d and ~57^3 local domains in 3d).

#define MV 26624

void
_SPUEAR_subsort_pipeline_spu( sort_p_pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline ) {
  MEM_PTR( const particle_t, 128 ) p_src = args->aux_p;
  MEM_PTR(       particle_t, 128 ) p_dst = args->p;
  int v0 = P2V( pipeline_rank,   n_pipeline, args->vl, args->vh );
  int v1 = P2V( pipeline_rank+1, n_pipeline, args->vl, args->vh );
  int i0, i1, i, j, v, sum, count;
  int b, nb;

  particle_t * RESTRICT ALIGNED(128) p_in;  SPU_MALLOC(p_in,  LOOKAHEAD, 128);
  particle_t * RESTRICT ALIGNED(128) p_out; SPU_MALLOC(p_out, LOOKAHEAD, 128);
  int * RESTRICT ALIGNED(128) next;         SPU_MALLOC(next,      MV,    128);
  int * RESTRICT              partition;    SPU_MALLOC(partition, MV+1+3,128);

  // +3 in partition to allow for EA alignment
  // FIXME: MV in next assumes MV>(n_pipeline+1)

  // This pipeline sorts particles in this voxel range in [i0,i1) in
  // the aux array.  These particles are in voxels [v0,v1).
  if( pipeline_rank==0            ) v0 = 0;
  if( pipeline_rank==n_pipeline-1 ) v1 = args->n_voxel;
  if( v1-v0 > MV ) return; // FIXME: DIAGNOSTIC ABORT
  mfc_get( next+pipeline_rank,
           args->coarse_partition+sizeof(int)*pipeline_rank,
           sizeof(int), 0, 0, 0 );
  mfc_get( next+pipeline_rank+1,
           args->coarse_partition+sizeof(int)*(pipeline_rank+1),
           sizeof(int), 1, 0, 0 );
  mfc_write_tag_mask( (1<<0)|(1<<1) ); mfc_read_tag_status_all();          
  i0 = next[pipeline_rank  ];
  i1 = next[pipeline_rank+1];

  next -= v0;             // Index next from v0.
  partition += (v0&3)-v0; // Index partition from v0 (and align ls with ea).

  // Begin getting the voxels for the initial particle block
  for( b=0; b<LOOKAHEAD; b++ )
    if( LIKELY(i0+b<i1) )
      mfc_get( &p_in[b].i, p_src+sizeof(particle_t)*(i0+b)+3*sizeof(int32_t),
               sizeof(int32_t), b, 0, 0 );

  // Clear fine grained count
  CLEAR( &next[v0], v1-v0 );

  // Fine grained count
  b = 0;
  for( i=i0; i<i1; i++ ) {
    mfc_write_tag_mask( 1<<b ); mfc_read_tag_status_all();
    next[ p_in[b].i ]++;
    if( LIKELY( i+LOOKAHEAD<i1 ) )
      mfc_get( &p_in[b].i,
               p_src+sizeof(particle_t)*(i+LOOKAHEAD)+3*sizeof(int32_t),
               sizeof(int32_t), b, 0, 0 );
    b = (b+1) & (LOOKAHEAD-1);
  }

  // Compute the partitioning
  sum = i0;
  for( v=v0; v<v1; v++ ) {
    count = next[v];
    next[v] = sum;
    partition[v] = sum;
    sum += count;
  }
  partition[v1] = sum; // All subsorts who write this agree

  // Begin getting particles for the initial particle block
  for( b=0; b<LOOKAHEAD; b++ )
    if( LIKELY(i0+b<i1) )
      mfc_get( p_in+b, p_src+sizeof(particle_t)*(i0+b), sizeof(particle_t),
               b, 0, 0 );

  // Write the ls partition to ea partition
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
    mfc_put( partition+v, args->partition+sizeof(int)*v, sizeof(int)*nb,
             b, 0, 0 ); b = (b+1) & 31;
  }

  // Local fine grained sort
  b = 0;
  for( i=i0; i<i1; i++ ) {

    // Wait for the load of this particle and the store of the 
    // corresponding particle in the previous block to complete,
    // determine where to store the arrived particle, copy it into a
    // local writeback buffer, begin the load of the corresponding
    // particle in the next block and the store of the particle in
    // this block.

    mfc_write_tag_mask( (1<<b)|(1<<(b+LOOKAHEAD)) ); mfc_read_tag_status_all(); 
    j = next[ p_in[b].i ]++;     
    ((vec_float4 *)(p_out+b))[0] = ((vec_float4 *)(p_in+b))[0];
    ((vec_float4 *)(p_out+b))[1] = ((vec_float4 *)(p_in+b))[1];
    mfc_put( p_out+b, p_dst+sizeof(particle_t)*j, sizeof(particle_t),
             b+LOOKAHEAD, 0, 0 );
    if( LIKELY(i+LOOKAHEAD<i1) )
      mfc_get( p_in+b,
               p_src+sizeof(particle_t)*(i+LOOKAHEAD), sizeof(particle_t),
               b, 0, 0 );
    b = (b+1) & (LOOKAHEAD-1);
  }
}

