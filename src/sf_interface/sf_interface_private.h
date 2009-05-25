#ifndef _sf_interface_private_h_
#define _sf_interface_private_h_

#ifndef IN_sf_interface
#error "Do not include sf_interface_private.h; include sf_interface.h"
#endif

#include "sf_interface.h"

///////////////////////////////////////////////////////////////////////////////
// load_interpolator_pipeline interface

typedef struct load_interpolator_pipeline_args {

  MEM_PTR( interpolator_t, 128 ) fi;
  MEM_PTR( const field_t,  128 ) f;
  MEM_PTR( const int64_t,  128 ) nb;
  int nx;
  int ny;
  int nz;

  PAD_STRUCT( 3*SIZEOF_MEM_PTR + 3*sizeof(int) )

} load_interpolator_pipeline_args_t;

PROTOTYPE_PIPELINE( load_interpolator, load_interpolator_pipeline_args_t );

///////////////////////////////////////////////////////////////////////////////
// clear_accumulators_pipeline interface

// Pipelines are be assigned accumulator blocks in multiples of 256
// (16KB) which is particularly convenient on Cell.  The pipeline
// dispatcher will handle any stragglers.

enum { accumulators_n_block = 256 };

typedef struct accumulators_pipeline_args {

  MEM_PTR( accumulator_t, 128) a; // First accumulator to reduce
  int n;                          // Number of accumulators to reduce
  int n_array;                    // Number of accumulator arrays
  int s_array;                    // Stride between each array

  PAD_STRUCT( SIZEOF_MEM_PTR + 3*sizeof(int) )

} accumulators_pipeline_args_t;

PROTOTYPE_PIPELINE( clear_accumulators,  accumulators_pipeline_args_t );
PROTOTYPE_PIPELINE( reduce_accumulators, accumulators_pipeline_args_t );

///////////////////////////////////////////////////////////////////////////////

typedef struct unload_accumulator_pipeline_args {

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

PROTOTYPE_PIPELINE( unload_accumulator, unload_accumulator_pipeline_args_t );

#undef FOR_SPU

#endif // _sf_interface_private_h_
