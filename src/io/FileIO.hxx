/*
	Definition of FileIO class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$
	vim: set ts=3 :
*/

#ifndef FileIO_hxx
#define FileIO_hxx

#include <FileIOData.hxx>

/*!
	\class FileIO FileIO.h
	\brief  provides...
*/
template<class ReadWritePolicy> class FileIO_T
	: public ReadWritePolicy
	{
	public:

		//! Constructor
		FileIO_T() {}

		//! Destructor
		~FileIO_T() {}

		// public interface inherited from policy

	private:

	}; // class FileIO_T

#include <StandardIOPolicy.hxx>
#include <P2PIOPolicy.hxx>

#if defined USE_AAIS_MP

#if defined HOST_BUILD
typedef FileIO_T<StandardIOPolicy> FileIO;
#else
typedef FileIO_T<P2PIOPolicy<true> > FileIOSwapped;
typedef FileIO_T<P2PIOPolicy<false> > FileIO;
#endif // BUILD

#else
typedef FileIO_T<StandardIOPolicy> FileIO;
#endif // MP Implementation

#endif // FileIO_hxx
