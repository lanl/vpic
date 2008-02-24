/*
 * Written by Ben Bergen
 * Computational Physics Group (CCS-2)
 * Computer, Computational, and Statistical Sciences Division
 * Los Alamos National Laboratory
 * Original - February 2008
*/

#include <spu_mfcio.h>
#include <spe_events.h>

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

/* overlay functions */
extern int o1_test1(int);
extern int o1_test2(int);
extern int o2_test1(int);
extern int o2_test2(int);

#define TEST(x, y)                     \
	if ((x) != (y)) {                    \
		printf("olay func call failed! expecting %d, got %d\n", x, y);   \
		return (1);                        \
	}

int main (uint64_t spu_id, uint64_t argp, uint64_t envp) {
	uint32_t msg;

	int rc;

	/* bring in overlay 1 */
	rc = o1_test1(101);
	TEST(1,rc);

	rc = o1_test2(102);
	TEST(2,rc);

	/* bring in overlay 2 */
	rc = o2_test1(201);
	TEST(11, rc);

	rc = o2_test2(202);
	TEST(12, rc);

	/* back to overlay 1 */
	rc = o1_test1(301);
	TEST(1, rc);

	rc = o1_test2(302);
	TEST(2, rc);

	do {
		msg = spu_read_in_mbox();

		switch(msg) {
			case READ_ARGS:
				//mfc_get(args, argp, sizeof(*args), 31, 0, 0);
				//mfc_write_tag_mask((1<<31));
				//mfc_read_tag_status_all();
				break;
		} // switch

	} while(msg != END_EVENT_LOOP);

	return 0;
} // main
