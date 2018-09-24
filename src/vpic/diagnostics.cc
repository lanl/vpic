/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "vpic.h"

# define RANK_TO_INDEX(rank,ix) BEGIN_PRIMITIVE {        \
	int _ix, _iy, _iz;                                    \
	_ix  = (rank);        /* ix = ix+gpx*( iy+gpy*iz ) */ \
	_iy  = _ix/int(px);   /* iy = iy+gpy*iz */            \
	_ix -= _iy*int(px);   /* ix = ix */                   \
	_iz  = _iy/int(py);   /* iz = iz */                   \
	_iy -= _iz*int(py);   /* iy = iy */                   \
	(ix) = _ix;                                           \
} END_PRIMITIVE 

/*------------------------------------------------------------------------------
 * Compute poynting flux at left boundary
 *
 * Note: This should maybe be made more general by adding which face to use.
 *
 * Note: This implementation is taken from Brian's GB deck
 * 
 * Inputs:
 *   e0 Peak instantaneous E field in "natural units"
 *----------------------------------------------------------------------------*/
// FIXME: Re-visit this implementation and replace with something more modern
// and concise.  People currently have versions in their input decks..
double vpic_simulation::poynting_flux(double e0) {
	double psum, gpsum;
	int stride = (grid->ny-1)*(grid->nz-1);

	float * pvec = new float[stride];

	if(!pvec) {
		ERROR(("Failed pvec allocation in poynting flux diagnostic"));
	} // if
	
	memset(pvec, 0, stride);

	int ix, k1, k2;
	RANK_TO_INDEX( int(rank()), ix);

	// Compute Poynting for domains on left of box
	if(ix==0) {
		for(int j(1); j<grid->ny; j++) {
			for(int k(1); k<grid->nz; k++) {
				k1 = INDEX_FORTRAN_3(1, j+1, k+1, 0, grid->nx+1,
					0, grid->ny+1, 0, grid->nz+1);
				k2 = INDEX_FORTRAN_3(2, j+1, k+1, 0, grid->nx+1,
					0, grid->ny+1, 0, grid->nz+1);

				pvec[(j-1)*(grid->nz-1)+k-1] =
					(field(k2).ey*0.5*(field(k1).cbz+field(k2).cbz) -
					field(k2).ez*0.5*(field(k1).cby+field(k2).cby)) /
					(grid->cvac*grid->cvac*e0*e0);
			} // for
		} // for
	} // if

	// Sum flux contributions
	for(int i(0); i<stride; i++) {
		psum += pvec[i];
	} // for
	
	// Collect sums over all ranks
	mp_allsum_d( &psum, &gpsum, 1 );

	// Normalize flux
	gpsum /= stride*py*pz;

	delete[] pvec;

	return gpsum;
} // poynting_flux

#undef RANK_TO_INDEX
