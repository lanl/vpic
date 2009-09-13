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
// processed in boundary_p() in boundary_p.cxx after particles are
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
// Revamped by KJB, May 2008, Sep 2009

#define IN_boundary
#include "boundary_private.h"
 
/* Private interface ********************************************************/

typedef struct maxwellian_reflux {
  species_t * sp_list;
  mt_rng_t  * rng;
  float     * ut_para;
  float     * ut_perp;
} maxwellian_reflux_t;

#ifndef M_SQRT2
#define M_SQRT2 (1.4142135623730950488016887242096981)
#endif

/* FIXME: DON'T IGNORE MAX_PI */
int
interact_maxwellian_reflux( maxwellian_reflux_t * RESTRICT mr,
                            species_t           * RESTRICT sp,
                            particle_t          * RESTRICT p, 
                            particle_mover_t    * RESTRICT pm,
                            particle_injector_t * RESTRICT pi,
                            int                            max_pi,
                            int                            face ) {
  const grid_t   * RESTRICT g   = sp->g;
  /**/  mt_rng_t * RESTRICT rng = mr->rng;

  const int32_t sp_id   = sp->id;
  const float   ut_para = mr->ut_para[sp_id]; 
  const float   ut_perp = mr->ut_perp[sp_id];

  float u[3];                // u0 = para, u1 & u2 = perp
  float ux, uy, uz;          // x, y, z normalized momenta
  float dispx, dispy, dispz; // Particle displacement
  float ratio;  

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
  ratio = p->ux*p->ux + p->uy*p->uy + p->uz*p->uz;
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

  pi->dx    = p->dx;
  pi->dy    = p->dy;
  pi->dz    = p->dz;
  pi->i     = p->i;
  pi->ux    = ux;
  pi->uy    = uy;
  pi->uz    = uz;
  pi->w     = p->w;
  pi->dispx = dispx;
  pi->dispy = dispy;
  pi->dispz = dispz;
  pi->sp_id = sp_id;
  return 1;
}

void
checkpt_maxwellian_reflux( const particle_bc_t * RESTRICT pbc ) {
  const maxwellian_reflux_t * RESTRICT mr =
    (const maxwellian_reflux_t *)pbc->params;
  CHECKPT( mr, 1 );
  CHECKPT_PTR( mr->sp_list );
  CHECKPT_PTR( mr->rng     );
  CHECKPT( mr->ut_para, num_species( mr->sp_list ) );
  CHECKPT( mr->ut_perp, num_species( mr->sp_list ) );
  checkpt_particle_bc_internal( pbc );
}

particle_bc_t *
restore_maxwellian_reflux( void ) {
  maxwellian_reflux_t * mr;
  RESTORE( mr );
  RESTORE_PTR( mr->sp_list );
  RESTORE_PTR( mr->rng     );
  RESTORE( mr->ut_para );
  RESTORE( mr->ut_perp );
  return restore_particle_bc_internal( mr );
}

void
delete_maxwellian_reflux( particle_bc_t * RESTRICT pbc ) {
  FREE( pbc->params );
  delete_particle_bc_internal( pbc );
}

/* Public interface *********************************************************/

particle_bc_t *
maxwellian_reflux( species_t * RESTRICT sp_list,
                   mt_rng_t  ** rng ) {
  if( !sp_list || !rng ) ERROR(( "Bad args" ));
  maxwellian_reflux_t * mr;
  MALLOC( mr, 1 );
  mr->sp_list = sp_list;
  mr->rng     = rng[0];
  MALLOC( mr->ut_para, num_species( mr->sp_list ) );
  MALLOC( mr->ut_perp, num_species( mr->sp_list ) );
  CLEAR( mr->ut_para, num_species( mr->sp_list ) );
  CLEAR( mr->ut_perp, num_species( mr->sp_list ) );
  return new_particle_bc_internal( mr,
                                   (particle_bc_func_t)interact_maxwellian_reflux,
                                   delete_maxwellian_reflux,
                                   (checkpt_func_t)checkpt_maxwellian_reflux,
                                   (restore_func_t)restore_maxwellian_reflux,
                                   NULL );
}

/* FIXME: NOMINALLY, THIS INTERFACE SHOULD TAKE kT */
void
set_reflux_temp( /**/  particle_bc_t * RESTRICT pbc,
                 const species_t     * RESTRICT sp,
                 float ut_para,
                 float ut_perp ) {
  if( !pbc || !sp || ut_para<0 || ut_perp<0 ) ERROR(( "Bad args" ));
  maxwellian_reflux_t * RESTRICT mr = (maxwellian_reflux_t *)pbc->params;
  mr->ut_para[sp->id] = ut_para;
  mr->ut_perp[sp->id] = ut_perp;
}

