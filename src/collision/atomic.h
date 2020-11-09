#ifndef _atomic_h_
#define _atomic_h_
#include "collision.h"

typedef struct atomic_properties {
  float alpha_m2;    // alpha^-2 = 2*Ry/me c^2
  float Nca02;       // n0 c a_0^2 has units of 1/time

  // Atomic data routines. They are denoted by _ab_xx where a and b are the
  // species [ e: electrons, z: primary ions, z1: parent ions ] and xx are the processes:
  //  el : elastic scattering
  //  cx : charge exchange (implemented as backward-dominated scattering)
  //  ex : excitation (inelastic scattering)
  //  in : ionization
  //  rx : recombination
  // Note that elastic scattering cross-sections only make sense when z=0, such
  // that the primary ion is actually a neutral atom. Otherwise, elastic
  // collisions will likely be dominated by Coulomb collisions and a different
  // collisional model should be applied.

  // Cross sections. Each should take on input the collision energy in units of
  // Ry and return the cross-section in units of a0^2.
  float (*sigma_pt_el)(float E);
  float (*sigma_pt_cx)(float E);
  float (*sigma_pt_ex)(float E);
  float (*sigma_et_el)(float E);
  float (*sigma_et_ex)(float E);
  float (*sigma_et_in)(float E);
  float (*sigma_ep_rx)(float E);
  float (*sigma_tt_el)(float E);
  float (*sigma_tt_ex)(float E);

  // Delta = tan(theta/2). Each should take on input the collision energy.
  float (*delta_pt_el)(rng_t* rng, float E);
  float (*delta_pt_cx)(rng_t* rng, float E);
  float (*delta_et_el)(rng_t* rng, float E);
  float (*delta_tt_el)(rng_t* rng, float E);

  // Gamma is the coefficient of restitution for inelastic processes, defined
  // by gamma = sqrt(KE'/KE) where KE' is the post-collision energy and KE
  // is the collision energy. Each should take on input the collision energy.
  float (*gamma_pt_ex)(rng_t* rng, float E);
  float (*gamma_et_ex)(rng_t* rng, float E);
  float (*gamma_tt_ex)(rng_t* rng, float E);

  // For ionization, we also set the secondary energy (in Ry)
  float (*gamma_et_in)(rng_t* rng, float Eprimary, float *Esecondary);

} atomic_properties_t;

enum atomic_processes {
  parent_target_elastic          = 1 << 0,
  parent_target_charge_exchange  = 1 << 1,
  parent_target_inelastic        = 1 << 2,
  electron_target_elastic        = 1 << 3,
  electron_target_inelastic      = 1 << 4,
  electron_target_ionization     = 1 << 5,
  electron_parent_recombination  = 1 << 6,
  target_target_elastic          = 1 << 7,
  target_target_inelastic        = 1 << 8,
  three_body_recombination       = 1 << 9,
  all_atomic_processes           = 0xFFFF
};

BEGIN_C_DECLS

collision_op_t *
append_atomic_collision_operators( /**/  unsigned int                   processes,
                                   /**/  species_t           * RESTRICT electron,
                                   /**/  species_t           * RESTRICT target,
                                   /**/  species_t           * RESTRICT parent,
                                   const atomic_properties_t * RESTRICT base,
                                   /**/  field_array_t       * RESTRICT fa,
                                   /**/  rng_pool_t          * RESTRICT rp,
                                   /**/  collision_op_t     **          cop_list,
                                   const double                         Na02,
                                   const double                         alpha,
                                   const int                            interval,
                                   const double                         sample,
                                   const double                         Emax
                                 );

END_C_DECLS

#endif /* _atomic_h_ */
