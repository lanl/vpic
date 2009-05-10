#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include "../sf_interface_private.h"

#include <spu_mfcio.h>

#define NB accumulators_n_block

void
_SPUEAR_clear_accumulators_pipeline_spu(
    MEM_PTR( accumulators_pipeline_args_t, 128 ) argp,
    int pipeline_rank,
    int n_pipeline ) {
  DECLARE_ALIGNED_ARRAY( accumulators_pipeline_args_t, 128, args,  1  );
  DECLARE_ALIGNED_ARRAY( accumulator_t,                128, zeros, NB );
  MEM_PTR( accumulator_t, 128 ) a;
  int n, n_array, s_array;
  int i, c;

  // Start getting the pipeline args from the dispatcher
  mfc_get( args, argp, sizeof(*args), 31, 0, 0 );

  // While waiting for the args, clear the local accumulator block
  CLEAR( zeros, NB );

  // Finish getting the pipeline args
  mfc_write_tag_mask( (1<<31) );
  mfc_read_tag_status_all();
  a       = args->a;
  n       = args->n;
  n_array = args->n_array;
  s_array = args->s_array * sizeof(accumulator_t);

  // Zero out the accumulators assigned to this pipeline
  DISTRIBUTE( n, NB, pipeline_rank, n_pipeline, i, n );
  a += i*sizeof(accumulator_t);
  c = 0;

  for( ; n_array; n_array--, a += s_array )
    for( i=0; i<n; i+=NB )
      mfc_put( zeros, a+i*sizeof(accumulator_t),
               NB*sizeof(accumulator_t), c, 0, 0 ), c = (c+1)&0x1f;

  // Wait for the zeroing to complete
  mfc_write_tag_mask( 0xffffffff );
  mfc_read_tag_status_all();
}

