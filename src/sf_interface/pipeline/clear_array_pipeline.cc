#define IN_sf_interface

#include "sf_interface_pipeline.h"

#include "../sf_interface_private.h"

#include "../../util/pipelines/pipelines_exec.h"

void
clear_array_pipeline_scalar( reduce_pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline )
{
  float * ALIGNED(16) a = args->a;

  int n       = args->n;
  int n_array = args->n_array;
  int s_array = args->s_array;
  int i;

  DISTRIBUTE( n, args->n_block, pipeline_rank, n_pipeline, i, n );

  a += i;

  for( ; n_array; n_array--, a += s_array )
  {
    CLEAR( a, n );
  }
}

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "V4 version not hooked up yet."

#endif

#define VOX( x, y, z ) VOXEL( x, y, z, aa->g->nx, aa->g->ny, aa->g->nz )

void
clear_accumulator_array_pipeline( accumulator_array_t * RESTRICT aa )
{
  DECLARE_ALIGNED_ARRAY( reduce_pipeline_args_t, 128, args, 1 );

  int i0;
  int na, nfloats;

  if ( ! aa )
  {
    ERROR( ( "Bad args." ) );
  }

  i0 = ( VOX( 1, 1, 1 ) / 2 ) * 2; // Round i0 down to even for 128 byte align.

  na = ( ( ( VOX( aa->g->nx, aa->g->ny, aa->g->nz ) - i0 + 1 ) + 1 ) / 2 ) * 2;

  nfloats       = sizeof(accumulator_t) / sizeof(float);

  args->a       = (float *) ( aa->a + i0 );
  args->n       = na * nfloats;
  args->n_array = aa->n_pipeline + 1;
  args->s_array = aa->stride * nfloats;
  args->n_block = accumulators_n_block;

  EXEC_PIPELINES( clear_array, args, 0 );

  WAIT_PIPELINES();
}

#undef VOX

void
clear_hydro_array_pipeline( hydro_array_t * RESTRICT ha )
{
  DECLARE_ALIGNED_ARRAY( reduce_pipeline_args_t, 128, args, 1 );

  int nfloats;

  if ( ! ha )
  {
    ERROR( ( "Bad args." ) );
  }

  nfloats       = sizeof(hydro_t) / sizeof(float);

  args->a       = (float *) ( ha->h );
  args->n       = ha->g->nv * nfloats;
  args->n_array = ha->n_pipeline + 1;
  args->s_array = ha->stride * nfloats;
  args->n_block = hydro_n_block;

  EXEC_PIPELINES( clear_array, args, 0 );

  WAIT_PIPELINES();
}
