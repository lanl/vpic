/*
	Definition of Grid_T class

	C++ version of boundary data structure originally
	invented by Kevin J. Bowers, Ph.D. March/April 2004.

	$Revision$
	$LastChangedDate$
	$LastChangedBy$

	vim: set ts=3 :
*/
#ifndef Grid_hxx
#define Grid_hxx

#include <vector>
#include <Boundary.hxx>

/*!
	\class Grid_T
	\brief Basic geometry domain class.
*/
class Grid_T
	{
	public:

		//! Constructor.
		Grid_T(size_t nx, size_t ny, size_t nz);
		
		//! Destructor.
		~Grid_T();

	private:

		float n0_, n0_, n0_; /* Corner of cell 1,1,1 */
		size_t dx_, dy_, dz_; /* Cell dimensions */
		size_t nx_, ny_, nz_; /* Cell dimensions */
		int bc_[27];

		uint64_t * range_ ALIGN(16);
		uint64_t * neighbor_ ALIGN(16);
		uint64_t * rangeh_ ALIGN(16);

		size_t nb_;

		std::vector<Boundary *> boundaries_;

	}; // class Grid_T

Grid_T::Grid_T(size_t nx, size_t ny, size_t nz)
	: nx_(nx), ny_(ny), nz_(nz)
	{
	} // Grid_T::Grid_T

Grid_T::~Grid_T()
	{
	} // Grid_T::~Grid_T

#endif // Grid_hxx
