/*
	Definition of DMPConnection class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$
	vim: set ts=3 :
*/

#ifndef DMPConnection_hxx
#define DMPConnection_hxx

/*!
	\class DMPConnection DMPConnection.h
	\brief  provides...
*/
template<class CommunicationPolicy> class DMPConnection_T
	: public CommunicationPolicy
	{
	public:

		static DMPConnection_T & instance() {
			static DMPConnection_T dmp;
			return dmp;
		} // instance

		// inherit public interface

	private:

		//! Constructor
		DMPConnection_T() {}

		//! Constructor
		DMPConnection_T(const DMPConnection_T & dmp) {}

		//! Destructor
		~DMPConnection_T() {}

	}; // class DMPConnection_T

#include "DMPPolicyMPI.hxx"

typedef DMPConnection_T<DMPPolicyMPI> DMPConnection;

#endif // DMPConnection_hxx
