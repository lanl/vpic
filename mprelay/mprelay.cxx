/*------------------------------------------------------------------------------
------------------------------------------------------------------------------*/

#include <ConnectionManager.hxx>
#include <MPRelay.hxx>

int main(int argc, char *argv[]) {

	// initialize everything
	ConnectionManager::instance().init(argc, argv);

	// Relay object
	MPRelay relay;

	// start relay
	relay.start();

	// finalize
	ConnectionManager::instance().finalize();

	return 0;
} // main
