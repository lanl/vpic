/*~-------------------------------------------------------------------------~~*
 *~-------------------------------------------------------------------------~~*/

#include <gtest/gtest.h>
#include "../src/vpic/vpic.h"

int main(int argc, char ** argv) {
	
	// Initialize the vpic runtime
	boot_services(&argc, &argv);

	// Initialize the GTest runtime
	::testing::InitGoogleTest(&argc, argv);
	int result = RUN_ALL_TESTS();

	// Shutdown the vpic runtime
	halt_services();

	return result;

} // main
