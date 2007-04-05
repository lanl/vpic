/*
	Definition of Boundary_T class

	C++ version of boundary data structure originally
	invented by Kevin J. Bowers, Ph.D. March/April 2004.

	$Revision$
	$LastChagnedDate: 2007-04-05 08:54:14 -0600 (Thu, 05 Apr 2007) $
	$LastChangedBy$

	vim: set ts=3 :
*/
#ifndef Boundary_h
#define Boundary_h

/*!
	\class Boundary Boundary.h
	\brief Boundary handler class.
*/
class Boundary_T
	{
	public:

		//! Constructor
		Boundary_T();

		//! Destructor
		~Boundary_T();

	private:

	}; // class Boundary_T

typedef Boundary_T Boundary;

Boundary_T::Boundary_T()
	{
	} // Boundary_T::Boundary_T

Boundary_T::~Boundary_T()
	{
	} // Boundary_T::~Boundary_T

#endif // Boundary_h
