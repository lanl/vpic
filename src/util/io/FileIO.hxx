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

#include "FileIOData.hxx"

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

		FileIOStatus open(const char * filename, FileIOMode mode)
			{ return ReadWritePolicy::open(filename, mode); }
		void close()
			{ return ReadWritePolicy::close(); }

		bool isOpen()
			{ return ReadWritePolicy::isOpen(); }

		int size()
			{ return ReadWritePolicy::size(); }

		template<typename T> void read(T * data, size_t elements)
			{ return ReadWritePolicy::read(data, elements); }
		template<typename T> void write(const T * data, size_t elements)
			{ return ReadWritePolicy::write(data, elements); }

	private:

	}; // class FileIO_T


#if defined USE_MPRELAY

#if defined HOST_BUILD
#include <StandardIOPolicy.hxx>

typedef FileIO_T<StandardIOPolicy> FileIO;
#else
#include <P2PIOPolicy.hxx>

//typedef FileIO_T<P2PIOPolicy<true> > FileIOSwapped;
//typedef FileIO_T<P2PIOPolicy<true> > FileIO;
typedef FileIO_T<P2PIOPolicy<true> > FileIO;
typedef FileIO_T<P2PIOPolicy<false> > FileIOUnswapped;
#endif // BUILD

#else
#include <StandardIOPolicy.hxx>
typedef FileIO_T<StandardIOPolicy> FileIO;
typedef FileIO_T<StandardIOPolicy> FileIOUnswapped;
#endif // MP Implementation

#endif // FileIO_hxx
