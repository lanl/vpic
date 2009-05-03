#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include "../sf_interface_private.h"

#include <spu_mfcio.h>

#define NB accumulators_n_block
#define MAX_ARRAY 12
#define PUT_CH_MIN MAX_ARRAY
#define PUT_CH_MAX 31

// DMA tag usage:
// 0:MAX_ARRAY-1 -> Getting accumulator blocks
// MAX_ARRAY:31  -> Writing accumulator blocks

// Local store single precision accumulator data is organized
// 0 1 2 3 / 4 5 6 7 / 8 9 10 11 / x x x x / 12 13 14 15 / ...
// ---- A[0] -----------------------------   ---- A[1] -------
//
// Local store double precision accumulator data is organized
// 0 2 / 1 3 / 4 6 / 5 7 / 8 10 / 9 11 / 12 14 / 13 15 / ...
// ---- A[0] -------------------------   ----- A[1] --------

void
_SPUEAR_reduce_accumulators_pipeline_spu(
    MEM_PTR( accumulators_pipeline_args_t, 128 ) argp,
    int pipeline_rank,
    int n_pipeline ) {

  vec_uchar16 perm = {  4, 5, 6, 7,  0, 1, 2, 3, 12,13,14,15,   8, 9,10,11 };
  vec_float4 f0, f1, f2;
  MEM_PTR( accumulator_t, 128 ) a;
  int sj, nj, i, i1, j, n, k, c;

  DECLARE_ALIGNED_ARRAY( accumulators_pipeline_args_t, 128, args, 1 );
  DECLARE_ALIGNED_ARRAY( vec_float4,  128, f, 4*NB*MAX_ARRAY ); // 192KB
  DECLARE_ALIGNED_ARRAY( vec_double2, 128, d, 6*NB           ); //  24KB

  // Get pipeline args from the dispatcher

  mfc_get( args, argp, sizeof(*args), 31, 0, 0 );
  mfc_write_tag_mask( (1<<31) );
  mfc_read_tag_status_all();

  // Determine which accumulators are reduced by this pipeline

  a  = args->a;
  nj = args->n_array; if( nj>MAX_ARRAY ) return; // FIXME: ABORT W/ DIAGS! 
  sj = args->stride;
  DISTRIBUTE( sj, NB, pipeline_rank, n_pipeline, i, i1 ); i1 += i;

  // For all blocks of voxels reduced by this pipeline

  c = PUT_CH_MIN;
  for( ; i<i1; i+=NB ) {

    // Begin loading all the single precision accumulators for this block

    for( j=0; j<nj; j++ )
      mfc_get( f+4*NB*j, a+(i+j*sj)*sizeof(accumulator_t),
               NB*sizeof(accumulator_t), j, 0, 0 );

    // Wait for block 0 to arrive and convert it to double precision
    
    mfc_write_tag_mask( 1<<0 );
    mfc_read_tag_status_all();
    
    for( n=0; n<NB; n++ ) {
      f0 = f[4*n+0];
      f1 = f[4*n+1];
      f2 = f[4*n+2];
      d[6*n+0] = spu_extend(             f0          );
      d[6*n+1] = spu_extend( spu_shuffle(f0,f0,perm) );
      d[6*n+2] = spu_extend(             f1          );
      d[6*n+3] = spu_extend( spu_shuffle(f1,f1,perm) );
      d[6*n+4] = spu_extend(             f2          );
      d[6*n+5] = spu_extend( spu_shuffle(f2,f2,perm) );
    }

    // For each remaining block, wait for it to arrive and then sum it
    // in double precision with block 0
 
    for( j=1; j<nj; j++ ) {

      mfc_write_tag_mask( 1<<j );
      mfc_read_tag_status_all();

      k = j*NB;
      for( n=0; n<NB; n++ ) {       
        f0       = f[4*(k+n)+0];
        f1       = f[4*(k+n)+1];
        f2       = f[4*(k+n)+2];
        d[6*n+0] = spu_add( d[6*n+0], spu_extend(             f0          ) );
        d[6*n+1] = spu_add( d[6*n+1], spu_extend( spu_shuffle(f0,f0,perm) ) );
        d[6*n+2] = spu_add( d[6*n+2], spu_extend(             f1          ) );
        d[6*n+3] = spu_add( d[6*n+3], spu_extend( spu_shuffle(f1,f1,perm) ) );
        d[6*n+4] = spu_add( d[6*n+4], spu_extend(             f2          ) );
        d[6*n+5] = spu_add( d[6*n+5], spu_extend( spu_shuffle(f2,f2,perm) ) );
      }

    }

    // Convert the now reduced block 0 back to single precision and store it

    for( n=0; n<NB; n++ ) {
      f0       = spu_roundtf( d[6*n+1] );
      f1       = spu_roundtf( d[6*n+3] );
      f2       = spu_roundtf( d[6*n+5] );
      f[4*n+0] = spu_or( spu_roundtf( d[6*n+0] ), spu_shuffle(f0,f0,perm) );
      f[4*n+1] = spu_or( spu_roundtf( d[6*n+2] ), spu_shuffle(f1,f1,perm) );
      f[4*n+2] = spu_or( spu_roundtf( d[6*n+4] ), spu_shuffle(f2,f2,perm) );
    }

    mfc_put( f, a+i*sizeof(accumulator_t), NB*sizeof(accumulator_t), c, 0, 0 );
    c++; if( c>PUT_CH_MAX ) c = PUT_CH_MIN;
  }
  
  // Wait for all stores to complete

  mfc_write_tag_mask( 0xffffffff );
  mfc_read_tag_status_all();
}
