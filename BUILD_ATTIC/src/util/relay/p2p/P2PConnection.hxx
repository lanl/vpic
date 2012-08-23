/*
	Definition of P2PConnection class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef P2PConnection_hxx
#define P2PConnection_hxx

/*!
	\class P2PConnection P2PConnection.h
	\brief  provides...
*/
template<class CommunicationPolicy> class P2PConnection_T
	: public CommunicationPolicy
	{
	public:

		static P2PConnection_T & instance() {
			static P2PConnection_T p2p;
			return p2p;
		} // instance

		// inherit public interface

	private:

		//! Constructor
		P2PConnection_T()
			: CommunicationPolicy()
			{}

		//! Constructor
		P2PConnection_T(const P2PConnection_T & p2p)
			: CommunicationPolicy()
			{}

		//! Destructor
		~P2PConnection_T() {}

	}; // class P2PConnection_T


#if defined USE_DACS_P2P
	#include <P2PPolicyDaCS.hxx>
	#if defined HOST_BUILD
		typedef P2PConnection_T<P2PPolicyDaCS<0> > P2PConnection;
	#else
		typedef P2PConnection_T<P2PPolicyDaCS<1> > P2PConnection;
	#endif // BUILD
#else
	#include <P2PPolicyMPI.hxx>
	#if defined HOST_BUILD
		typedef P2PConnection_T<P2PPolicyMPI<0> > P2PConnection;
	#else
		typedef P2PConnection_T<P2PPolicyMPI<1> > P2PConnection;
	#endif // BUILD
#endif // USE_DACS_P2P

#endif // P2PConnection_hxx
