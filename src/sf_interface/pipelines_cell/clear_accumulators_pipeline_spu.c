#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include <sf_interface_private.h>

#include <spu_mfcio.h>

void
_SPUEAR_clear_accumulators_pipeline_spu( accumulators_pipeline_args_t * args,
                                         int pipeline_rank,
                                         int n_pipeline ) {
  MEM_PTR( accumulator_t, 128 ) a = args->a;
  int n = args->n, n_array = args->n_array, s_array = args->s_array*sizeof(accumulator_t), i;
  DISTRIBUTE(n, accumulators_n_block, pipeline_rank, n_pipeline, i, n); a += i*sizeof(accumulator_t);

  accumulator_t * ALIGNED(128) zeros; int c = 0;
  SPU_MALLOC( zeros, accumulators_n_block, 128 );
  CLEAR( zeros, accumulators_n_block );

  for( ; n_array; n_array--, a+=s_array ) /* CLEAR(a,n); */
    for( i=0; i<n; i+=accumulators_n_block ) {
      mfc_put( zeros, a+i*sizeof(accumulator_t),
               accumulators_n_block*sizeof(accumulator_t), c, 0, 0 );
      c = (c+1)&0x1f;
    }
}

