/*
 * Written by Ben Bergen
 * Computational Physics Group (CCS-2)
 * Computer, Computational, and Statistical Sciences Division
 * Los Alamos National Laboratory
 * Original - February 2008
*/

#include <spu_mfcio.h>
#include "../pipelines/pipelines.h"

////////////////////////////////////////////////////////////////////////////////
// Main SPU program
//
// This will eventually replace the main in the particle directory, and will
// use overlays to call various operations, e.g., particle advance, sort, and
// field solve.  We should be able to do this with a single region under the
// current design.  Some of the structure of VPIC will need to be changed so
// that this event loop can be started outside of the particle advance.  This
// is the current structure.
////////////////////////////////////////////////////////////////////////////////

// Must match definition in pipelines_spu.c; Geddy Lee rules 
#define SPU_COMPLETE ((uint32_t)2112)

// Include the frandn Ziggurat table.  This is done here so that each
// SPE pipeline who uses frandn doesn't blow away local store for a
// redundant copy of frandn table.  This single instance of the frandn
// used by all blows away 0.5K of local store; in the long run, it
// would be better to use the SPU heap below for this table and load it
// in pipelines who need it (the current drandn table, which would blow
// away 8K of local store, virtually requires this).

#define IN_rng
#include "../rng/frandn_table.c"

// Declare the spe_heap

char _spu_heap[ 224*1024 ];
size_t _spu_heap_free;

int
main( uint64_t id,
      uint64_t main_argp,
      uint64_t main_envp ) {
  pipeline_func_t pipeline;
  MEM_PTR( void, 128 ) argp;
  char * args;
  int sz;
  int pipeline_rank, n_pipeline;

  for(;;) {

    // Read mailbox and interpret message as a function pointer
    // in the SPE symbol space and arguments for that function

    pipeline = (pipeline_func_t)spu_read_in_mbox();
    if( !pipeline ) break; // If we get a NULL pipeline, we are done
    argp = (((uint64_t)spu_read_in_mbox())<<32) |
            ((uint64_t)spu_read_in_mbox());
    sz            = spu_read_in_mbox();
    pipeline_rank = spu_read_in_mbox();
    n_pipeline    = spu_read_in_mbox();

    // Read the pipeline args from memory

    SPU_FREE_ALL();
    SPU_MALLOC( args, 1, sz );
    mfc_get( args, argp, sz, 31, 0, 0 );
    mfc_write_tag_mask( 1<<31 );
    mfc_read_tag_status_all();

    // Execute the pipeline

    pipeline( args, pipeline_rank, n_pipeline );

    // Wait for all actions done by the pipeline to complete

    mfc_write_tag_mask( 0xffffffff );
    mfc_read_tag_status_all();

    // Notify caller of completion

    spu_write_out_mbox(SPU_COMPLETE);

  } // for

  return 0;
} // main

