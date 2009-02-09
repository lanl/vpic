/*
 * Written by Ben Bergen
 * Computational Physics Group (CCS-2)
 * Computer, Computational, and Statistical Sciences Division
 * Los Alamos National Laboratory
 * Original - February 2008
*/

#include <spu_mfcio.h>

////////////////////////////////////////////////////////////////////////////////
// Main SPU program
//
// This will eventually replace the main in the particle directory, and will
// use overlays to call various operations, e.g., particle advance, sort, and
// field solve.  We should be able to do this with a single region under the
// current design.  Some of the structure of VPIC will need to be changed so
// that this event loop can be started outside of the particle advance.  This
// is the current structure.
//
////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>

typedef void (*fptr)(uint64_t id, uint64_t argp, uint64_t envp);

int main (uint64_t spu_id, uint64_t argp, uint64_t envp) {

	// Read mailbox and interpret message as a function pointer
	// in the SPE symbol space
	fptr f = (fptr)spu_read_in_mbox();

	while(f != NULL) {
		// Execute function
		f(id, argp, envp);

		// Write completion message to host
		spu_write_out_mbox(COMPLETE);

		// Read next function
		f = (fptr)spu_read_in_mbox();
	} // while

	return 0;
} // main
