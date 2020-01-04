#ifndef _sf_interface_pipeline_h_
#define _sf_interface_pipeline_h_

#ifndef IN_sf_interface
#error "Do not include sf_interface_pipeline.h; include sf_interface.h"
#endif

#include "../sf_interface.h"

///////////////////////////////////////////////////////////////////////////////
// load_interpolator_pipeline interface

typedef struct load_interpolator_pipeline_args
{
  MEM_PTR( interpolator_t, 128 ) fi;
  MEM_PTR( const field_t,  128 ) f;
  MEM_PTR( const int64_t,  128 ) nb;
  int nx;
  int ny;
  int nz;

  PAD_STRUCT( 3*SIZEOF_MEM_PTR + 3*sizeof(int) )

} load_interpolator_pipeline_args_t;

void
load_interpolator_pipeline_scalar( load_interpolator_pipeline_args_t * args,
				   int pipeline_rank,
				   int n_pipeline );

void
load_interpolator_pipeline_v4( load_interpolator_pipeline_args_t * args,
                               int pipeline_rank,
                               int n_pipeline );

///////////////////////////////////////////////////////////////////////////////

typedef struct unload_accumulator_pipeline_args
{
  MEM_PTR( field_t, 128 ) f;             // Reduce accumulators to this
  MEM_PTR( const accumulator_t, 128 ) a; // Accumulator array to reduce
  int nx;                                // Local domain x-resolution
  int ny;                                // Local domain y-resolution
  int nz;                                // Local domain z-resolution
  float cx;                              // x-axis coupling constant
  float cy;                              // y-axis coupling constant
  float cz;                              // z-axis coupling constant

  PAD_STRUCT( 2*SIZEOF_MEM_PTR + 3*sizeof(int) + 3*sizeof(float) )

} unload_accumulator_pipeline_args_t;

void
unload_accumulator_pipeline_scalar( unload_accumulator_pipeline_args_t * args,
                                    int pipeline_rank,
                                    int n_pipeline );

///////////////////////////////////////////////////////////////////////////////
// clear_array_pipeline interface

// Pipelines are be assigned accumulator blocks in multiples of 256
// (16KB) which is particularly convenient on Cell.  The pipeline
// dispatcher will handle any stragglers.

enum { accumulators_n_block = 256, hydro_n_block = 1 };

typedef struct reduce_pipeline_args
{
  MEM_PTR(float, 128) a;          // First array element to reduce
  int n;                          // Number of array elements to reduce
  int n_array;                    // Number of pipeline arrays
  int s_array;                    // Stride between each array
  int n_block;                    // Number of floats/block.

  PAD_STRUCT( SIZEOF_MEM_PTR + 4*sizeof(int) )

} reduce_pipeline_args_t;

void
clear_array_pipeline_scalar( reduce_pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline );

///////////////////////////////////////////////////////////////////////////////
// reduce_array_pipeline interface

void
reduce_array_pipeline_scalar( reduce_pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline );

void
reduce_array_pipeline_v4( reduce_pipeline_args_t * args,
                          int pipeline_rank,
                          int n_pipeline );

#endif // _sf_interface_pipeline_h_
