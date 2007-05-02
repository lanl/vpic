/*------------------------------------------------------------------------------
------------------------------------------------------------------------------*/

#include <iostream>
#include <P2PConnection.hxx>
#include <DMPRelay.hxx>

int main(int argc, char *argv[]) {

	// Relay object
	DMPRelay relay;

	// get an instance of the manager
	P2PConnection & p2p = P2PConnection::instance();		

	// initialize everything
	p2p.init(argc, argv);

	// start relay
	relay.start();

	// finalize
	p2p.finalize();

	return 0;
} // main

