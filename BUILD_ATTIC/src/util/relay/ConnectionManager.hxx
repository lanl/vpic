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

/*!
	\struct ConnectionManager ConnectionManager.h
	\brief  provides...
*/
template<class ConnectionPolicy> class ConnectionManager_T
	: public ConnectionPolicy
	{
	public:

		static ConnectionManager_T & instance() {
			static ConnectionManager_T mgr;
			return mgr;
		} // instance

		// public interface inherited from policy

	private:

		ConnectionManager_T() {}
		ConnectionManager_T(const ConnectionManager_T & mgr) {}
		~ConnectionManager_T() {}

	}; // class ConnectionManager_T

#if defined USE_DACS_P2P
	#if defined HOST_BUILD
		#include "CMPolicyMPIDaCS.hxx"
		typedef ConnectionManager_T<CMPolicyMPIDaCS> ConnectionManager;
	#else
		#include "CMPolicyDaCS.hxx"
		typedef ConnectionManager_T<CMPolicyDaCS> ConnectionManager;
	#endif // BUILD TYPE
#else
	#include "CMPolicyMultipleContextMPI.hxx"
	typedef ConnectionManager_T<CMPolicyMultipleContextMPI>
		ConnectionManager;
#endif // P2P TYPE

#endif // ConnectionManager_hxx
