#include "boundary.h"

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
// Revamped by KJB, May 2008

#ifndef M_SQRT2
#define M_SQRT2 (1.4142135623730950488016887242096981)
#endif
 
// FIXME: MIGHT WANT TO CHECK THAT SP_ID DOESN'T OVERFLOW

void
maxwellian_reflux( void * _mr,
                   particle_t           * RESTRICT r, 
                   particle_mover_t     * RESTRICT pm,
                   field_t              * RESTRICT f,
                   accumulator_t        * RESTRICT a, 
                   const grid_t         * RESTRICT g,
                   species_t            * RESTRICT s, 
                   particle_injector_t ** RESTRICT ppi,
                   mt_rng_t             * RESTRICT rng,
                   int face ) {
  maxwellian_reflux_t * RESTRICT mr = (maxwellian_reflux_t *)_mr;
  int32_t sp_id = s->id;
  float ut_para = mr->ut_para[sp_id]; 
  float ut_perp = mr->ut_perp[sp_id];
  float u[3];                // u0 = para, u1 & u2 = perp
  float ux, uy, uz;          // x, y, z normalized momenta
  float dispx, dispy, dispz; // Particle displacement
  float ratio;  
  particle_injector_t * RESTRICT pi; 

  /**/                      // axis x  y  z 
  static const int perm[6][3] = { { 0, 1, 2 },   // -x face
                                  { 2, 0, 1 },   // -y face
                                  { 1, 2, 0 },   // -z face 
                                  { 0, 1, 2 },   // +x face
                                  { 2, 0, 1 },   // +y face
                                  { 1, 2, 0 } }; // +z face
  static const float scale[6] = {  M_SQRT2,  M_SQRT2,  M_SQRT2,
                                  -M_SQRT2, -M_SQRT2, -M_SQRT2 };

  // compute velocity of injected particle
  //
  // Suppose you have a Maxwellian at a boundary: p(u) ~ exp(-u^2/(2
  // ub^2)) where u is the || speed and ub is the thermal speed.  In a
  // time delta t, if the boundary has surface area delta A, there
  // will be
  //   
  //   p_inj(u) du ~ u exp(-u^2/(2 ub^2)) (delta t)(delta A) du
  //   
  // particles injected from the boundary between speeds u and
  // u+du. p_inj(u) is the distribution function we wish to sample.
  // It has a cumulative i distribution function
  //   
  //   cdf(u) = \int_0^u du p_inj(u) = 1 - exp(-u^2/(2 ub^2))
  //   
  // (I've adjusted the constants out front to give a proper cdf
  // ranging from 0 to 1, the range of h).
  //   
  // Let mu be a uniformly distributed random number from 0 to 1.
  // Setting cdf(u)=mu and solving for u gives the means for sampling
  // u:
  //   
  //   exp(-u^2/(2 ub^2)) = mu - 1 = mu
  //
  // (Note that 1-mu has same dist as mu.)  This implies that
  //
  //   u = sqrt(2) ub sqrt( -log(mu) ).
  //
  // Note that -log(mu) is an _exponentially_ distributed random
  // number.  In the long haul then, we should probably call
  // something like sqrt( mt_frande(rng) ) and mtrand is in a
  // better position to decide the best method to compute such
  // a rand (it too can be ziggurat accelerated) or even add a
  // maxwellian flux method to mtrand.

  // Note: This assumes ut_para > 0
  
  // Note: mt_frand_c1 may be more appropriate theoretically but
  // mt_frand is closer to original behavior (and guarantees a
  // non-zero normal displacement).

  u[0] = ut_para*scale[face]*sqrtf(-logf(mt_frand(rng)));
  u[1] = ut_perp*mt_frandn(rng);
  u[2] = ut_perp*mt_frandn(rng);
  ux   = u[perm[face][0]];
  uy   = u[perm[face][1]];
  uz   = u[perm[face][2]];

  // Compute the amount of aging to due of the refluxed particle.
  //
  // The displacement of the refluxed particle should be:
  //
  //   dr' = c dt u' (1-a) / gamma'
  //
  // where u' and gamma' refer to the refluxed 4-momentum and
  // a is when the particle's time step "age" when it hit the
  // boundary.
  //
  //   1-a = |remaining_dr| / ( c dt |u| / gamma )
  //
  // Thus, we have:
  //
  //   dr' = u' gamma |remaining_dr| / ( gamma' |u| )
  //
  // or:
  //
  //   dr' = u' sqrt(( (1+|u|^2) |remaining_dr|^2 ) / ( (1+|u'|^2) |u|^2 ))

  dispx = g->dx * pm->dispx;
  dispy = g->dy * pm->dispy;
  dispz = g->dz * pm->dispz;
  ratio = r->ux*r->ux + r->uy*r->uy + r->uz*r->uz;
  ratio = sqrtf( ( ( 1+ratio )*( dispx*dispx + dispy*dispy + dispz*dispz ) ) /
                 ( ( 1+(ux*ux+uy*uy+uz*uz) )*( FLT_MIN+ratio ) ) );
  dispx = ux * ratio * g->rdx;
  dispy = uy * ratio * g->rdy;
  dispz = uz * ratio * g->rdz;

  // If disp and u passed to this are consistent, ratio is sane in and
  // the displacment is non-zero in exact arithmetic.  However,
  // paranoid checking like the below can be done here if desired.
  //
  // if( ratio<=0 || ratio>=g->dt*g->cvac )
  //   WARNING(( "Bizarre behavior detected in maxwellian_reflux" ));

  pi        = (*ppi)++;
  pi->dx    = r->dx;
  pi->dy    = r->dy;
  pi->dz    = r->dz;
  pi->i     = r->i;
  pi->ux    = ux;
  pi->uy    = uy;
  pi->uz    = uz;
  pi->q     = r->q;
  pi->dispx = dispx;
  pi->dispy = dispy;
  pi->dispz = dispz;
  pi->sp_id = sp_id;
}

