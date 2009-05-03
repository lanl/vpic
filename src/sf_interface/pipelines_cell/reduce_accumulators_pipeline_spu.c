#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include "../sf_interface_private.h"

#include <spu_mfcio.h>

#define NB clear_accumulators_n_block

void
_SPUEAR_clear_accumulators_pipeline_spu(
    MEM_PTR( clear_accumulators_pipeline_args_t, 128 ) argp,
    int pipeline_rank,
    int n_pipeline ) {
  DECLARE_ALIGNED_ARRAY( clear_accumulators_pipeline_args_t, 128, args,  1  );
  DECLARE_ALIGNED_ARRAY( accumulator_t,                      128, zeros, NB );
  MEM_PTR( accumulator_t, 128 ) a;
  int i, n, c;

  // Start getting the pipeline args from the dispatcher
  mfc_get( args, argp, sizeof(*args), 31, 0, 0 );

  // While waiting for the args, clear the local accumulator block
  CLEAR( zeros, NB );

  // Finish getting the pipeline args
  mfc_write_tag_mask( (1<<31) );
  mfc_read_tag_status_all();

  // Zero out the accumulators assigned to this pipeline
  DISTRIBUTE( args->n, NB, pipeline_rank, n_pipeline, i, n );
  a = args->a + i*sizeof(accumulator_t);
  for( c=0; n>0; n-=NB, a+=NB*sizeof(accumulator_t), c=(c+1)&0x1f )
    mfc_put( zeros, a, NB*sizeof(accumulator_t), c, 0, 0 ); 

  // Wait for the zeroing to complete
  mfc_write_tag_mask( 0xffffffff );
  mfc_read_tag_status_all();
}

