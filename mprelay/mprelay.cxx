/*------------------------------------------------------------------------------
------------------------------------------------------------------------------*/

#include <MPRelay.hxx>

int main(int argc, char *argv[]) {

	// initialize everything
	P2PConnection::instance().init(argc, argv);

	// Relay object
	MPRelay relay;

	// start relay
	relay.start();

	// finalize
	P2PConnection::instance().finalize();

	return 0;
} // main
