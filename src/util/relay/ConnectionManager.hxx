/*
	Definition of ConnectionManager class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef ConnectionManager_hxx
#define ConnectionManager_hxx

#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <MPData.hxx>

/*!
	\struct ConnectionManager_T ConnectionManager_T.h
	\brief  provides...
*/
template<class ConnectionPolicy> class ConnectionManager_T
	: public ConnectionPolicy
	{
	public:

		static ConnectionManager_T & instance() {
			static ConnectionManager_T mpi;
			return mpi;
		} // instance

		// inherit public interface

	private:

		//! Constructor
		ConnectionManager_T()
			: ConnectionPolicy()
			{}

		ConnectionManager_T(const ConnectionManager_T & cm)
			: ConnectionPolicy()
			{}

		~ConnectionManager_T() {}

	}; // class ConnectionManager_T

#include <CMPolicyMPI.hxx>

typedef ConnectionManager_T<CMPolicyMPI> ConnectionManager;

#endif // ConnectionManager_hxx
