#include <boundary.h>
#include <stdio.h> // for debugging output

// Refluxing boundary condition on particles.  Calculate normalized
// momenta from the various perp and para temperatures.  Then, sample
// Maxwellian specra from Maxwellian distributsions that have the
// given normalized momenta.
//
// The perpendicular spectra are sampled from a bi-Maxwellian.  The
// parallel spectra are sampled from vpara = sqrt(2) vth sqrt(-log h)
// where h is uniformly distributed over range (0,1).  Strictly
// speaking, this method only works for non-relativistic distributions
// of edge particles.
//
// This routine requires that in the input deck one has defined arrays
// of id[], ut_para[], ut_perp[] for each particle species.  The
// routine creates a particle injector from the particle data passed
// to it.  Particle injection from custom boundary handlers is
// processed in boundary_p() in boundary_p.c after particles are
// exchanged between processors across domain boundaries.
//
// Note that the way maxwellian_reflux_t is defined we have a hard
// maximum of only 32 species handled by the boundary condition.  This
// is not dynamically sized at runtime since it makes it harder to get
// the restart working right if we do.  If more are needed, one needs
// to adjust the size of the maxwellian_reflux_t data space by
// changing the array sizes in boundary_handler.h
//
// Procedure to generate new particle displacements:
//
// 1) From face information, particle species type, and parameters for
//   species defined in the boundary handler data, obtain new
//   Maxwellian- refluxed particle momenta.
//
// 2) From the old particle displacements (stored in the mover data),
//   compute the new displacements according to:
//
// dx_new = dx_old * (ux_new/ux_old) * sqrt((1+|u_old|**2)/(1+|u_new|**2))
//
// Written by:  Brian J. Albright, X-1, LANL   April, 2005

#define maxwellian_rand(DEV)      ((DEV)*mt_normal_drand(rng))
#define SQRT_TWO                  1.41421356237
#define SET_VELOCITY(SIGN,X,Y,Z)  u##X=SIGN fabs(u0);u##Y=u1;u##Z=u2; break;
 
void maxwellian_reflux( void * _data, particle_t *r, 
                        particle_mover_t *pm, field_t *f, accumulator_t *a, 
			const grid_t *g, species_t *s, 
			particle_injector_t **ppi, mt_handle rng, int face ) {
  maxwellian_reflux_t * reflux_data = (maxwellian_reflux_t *)_data;
  float u0, u1, u2;       // u0 = para, u1 & u2 = perp
  float ux, uy, uz;       // x, y, z normalized momenta
  float h, gamma_ratio;  
  int ispec;
  particle_injector_t *pi; 

  // try to avoid having zero displacements
#define NO_ZERO_DISPLACEMENT 0

#if NO_ZERO_DISPLACEMENT 
  int pass=0; 
# define MAX_PASSES 5
# define EPS        1e-7
#endif 

  // sanity check: print handler data once at simulation start and restart
  {
    static int initted=0;
    if ( !initted && mp_rank(g->mp)==0 ) {
      initted=1;
      MESSAGE(("----------------boundary handler data-----------------")); 
      for ( ispec=0; ispec<reflux_data->nspec; ++ispec )
        MESSAGE(("Species: %i  ut_perp: %e  ut_para: %e",
	         reflux_data->id[ispec], reflux_data->ut_perp[ispec], 
                 reflux_data->ut_para[ispec])); 
      MESSAGE(("----------------boundary handler data-----------------")); 
    }
  }

  // obtain reflux data for particle according to particle species
  for ( ispec=0; ispec<reflux_data->nspec && s->id!=reflux_data->id[ispec]; ++ispec )
    ; 
  if ( ispec==reflux_data->nspec ) 
    ERROR(("Unknown species passed to boundary handler."));

#if NO_ZERO_DISPLACEMENT 
 sample_maxwellian_speeds: 
#endif 

  // compute velocity of injected particle
  do h=mt_drand(rng); while ( h==0 || h==1 ); 
  u0 = SQRT_TWO*reflux_data->ut_para[ispec]*sqrt(-log(h)); 
  u1 = maxwellian_rand(reflux_data->ut_perp[ispec]);
  u2 = maxwellian_rand(reflux_data->ut_perp[ispec]);

  // set velocity of particle according to which face it crosses
  switch( face ) {
  case 0:  SET_VELOCITY(+,x,y,z)  // negative x face
  case 1:  SET_VELOCITY(+,y,z,x)  // negative y face
  case 2:  SET_VELOCITY(+,z,x,y)  // negative z face
  case 3:  SET_VELOCITY(-,x,y,z)  // positive x face
  case 4:  SET_VELOCITY(-,y,z,x)  // positive y face
  case 5:  SET_VELOCITY(-,z,x,y)  // positive z face
      default: // should never happen
    ERROR(("Invalid face data passed to custom boundary routine."));
    ux = uy = uz = 0;
    break; 
  } 

  // insert values into particle injector
  pi = *ppi; 
  pi->dx = r->dx;
  pi->dy = r->dy;
  pi->dz = r->dz;
  pi->i  = r->i;
  pi->ux = ux;
  pi->uy = uy;
  pi->uz = uz;
  pi->q  = r->q;

  // compute particle displacements    
  gamma_ratio = sqrt((1+r->ux*r->ux+r->uy*r->uy+r->uz*r->uz)/(1+ux*ux+uy*uy+uz*uz)); 
  pi->dispx=pi->dispy=pi->dispz=0; 
# define ZWARN WARNING(("Zero velocity component encountered."))   
  if ( r->ux!=0 ) pi->dispx=pm->dispx*(ux/r->ux)*gamma_ratio; else ZWARN;
  if ( r->uy!=0 ) pi->dispy=pm->dispy*(uy/r->uy)*gamma_ratio; else ZWARN;
  if ( r->uz!=0 ) pi->dispz=pm->dispz*(uz/r->uz)*gamma_ratio; else ZWARN;
# undef ZWARN
  
#if NO_ZERO_DISPLACEMENT 
  if ( ++pass<MAX_PASSES ) {
    if ( fabs(pi->dispx)/g->dx<EPS || fabs(pi->dispy)/g->dy<EPS || fabs(pi->dispz)/g->dz<EPS )
      goto sample_maxwellian_speeds; 
  } else WARNING(("Infinitesmal particle displacement encountered.")); 
#endif 

  // increment pointer to injector
  (*ppi)++;  

  return;
}

