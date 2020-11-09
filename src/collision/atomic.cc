#include "atomic.h"

// We are splitting the colllsion process up into discrete sections
// because some processes only use part, and others need to do work
// in-between.
#define CMOV(a,b) if(t0<t1) a=b

#define BEGIN_COLLISION(AP,SPI,SPJ,PI,PJ) do{                                  \
  const species_t  * _spi = (SPI);                                             \
  const species_t  * _spj = (SPJ);                                             \
  /**/  particle_t * _pi  = (PI);                                              \
  /**/  particle_t * _pj  = (PJ);                                              \
                                                                               \
  const float mu   = (_spi->m * _spj->m)/(_spi->m + _spj->m);                  \
  const float mu_i = _spj->m/(_spi->m + _spj->m);                              \
  const float mu_j = _spi->m/(_spi->m + _spj->m);                              \
                                                                               \
  float E, dd, R, R1, ur, urx, ury, urz, tx, ty, tz, t0, t1, t2, stack[3];     \
  float cmx, cmy, cmz;                                                         \
  int d0, d1, d2;                                                              \
                                                                               \
  R  = 1;                                                                      \
  R1 = 0;                                                                      \
                                                                               \
  /* Relative velocity. Needed for all collisions. */                          \
  urx = _pi->ux - _pj->ux;                                                     \
  ury = _pi->uy - _pj->uy;                                                     \
  urz = _pi->uz - _pj->uz;                                                     \
                                                                               \
  /* Center of mass velocity. Needed for inelastic collisions. */              \
  cmx = mu_j*_pi->ux + mu_i*_pj->ux;                                           \
  cmy = mu_j*_pi->uy + mu_i*_pj->uy;                                           \
  cmz = mu_j*_pi->uz + mu_i*_pj->uz;                                           \
                                                                               \
  /* There are lots of ways to formulate T vector formation    */              \
  /* This has no branches (but uses L1 heavily)                */              \
                                                                               \
  t0 = urx*urx;      d0=0;       d1=1;       d2=2;       t1=t0;  ur  = t0;     \
  t0 = ury*ury; CMOV(d0,1); CMOV(d1,2); CMOV(d2,0); CMOV(t1,t0); ur += t0;     \
  t0 = urz*urz; CMOV(d0,2); CMOV(d1,0); CMOV(d2,1);              ur += t0;     \
  E  = (AP)->alpha_m2*mu*ur;                                                   \
  ur = sqrtf( ur );                                                            \
                                                                               \
  stack[0] = urx;                                                              \
  stack[1] = ury;                                                              \
  stack[2] = urz;                                                              \
  t1  = stack[d1];                                                             \
  t2  = stack[d2];                                                             \
  t0  = 1 / sqrtf( t1*t1 + t2*t2 + FLT_MIN );                                  \
  stack[d0] =  0;                                                              \
  stack[d1] =  t0*t2;                                                          \
  stack[d2] = -t0*t1;                                                          \
  tx = stack[0];                                                               \
  ty = stack[1];                                                               \
  tz = stack[2]

#define SCATTER(DELTA,PHI)                                                     \
  dd = (DELTA);                                                                \
  if(!isfinite(dd) || dd > 1.30e19f) dd = 1.30e19f;                            \
  t0 = 2*dd/(1+dd*dd);                                                         \
  t1 = (PHI);                                                                  \
  t2 = t0*sin(t1);                                                             \
  t1 = t0*ur*cos(t1);                                                          \
  t0 *= -dd;                                                                   \
                                                                               \
  stack[0] = (t0*urx + t1*tx) + t2*( ury*tz - urz*ty );                        \
  stack[1] = (t0*ury + t1*ty) + t2*( urz*tx - urx*tz );                        \
  stack[2] = (t0*urz + t1*tz) + t2*( urx*ty - ury*tx )

#define REDUCE_ENERGY(GAMMA)                                                   \
  R   = (GAMMA);                                                               \
  R1  = 1-R                                                                    \

#define FINALIZE_COLLISION(TYPE)                                               \
  d0 = (TYPE);                                                                 \
  if( d0 & 1 ){                                                                \
    _pi->ux = (_pi->ux + mu_i*stack[0])*R + R1*cmx;                            \
    _pi->uy = (_pi->uy + mu_i*stack[1])*R + R1*cmy;                            \
    _pi->uz = (_pi->uz + mu_i*stack[2])*R + R1*cmz;                            \
  }                                                                            \
  if( d0 & 2 ){                                                                \
    _pj->ux = (_pj->ux - mu_j*stack[0])*R + R1*cmx;                            \
    _pj->uy = (_pj->uy - mu_j*stack[1])*R + R1*cmy;                            \
    _pj->uz = (_pj->uz - mu_j*stack[2])*R + R1*cmz;                            \
  }                                                                            \
} while(0)

#define ELASTIC_COLLISION(AP, SPI, SPJ, PI, PJ, RNG, TYPE, DELTA)              \
  BEGIN_COLLISION(AP,SPI,SPJ,PI,PJ);                                           \
    SCATTER( DELTA(RNG, E), 2*M_PI*frand_c0(RNG) );                            \
  FINALIZE_COLLISION(TYPE)

#define INELASTIC_COLLISION(AP, SPI, SPJ, PI, PJ, RNG, TYPE, DELTA, GAMMA)     \
  BEGIN_COLLISION(AP,SPI,SPJ,PI,PJ);                                           \
    SCATTER( DELTA(RNG, E), 2*M_PI*frand_c0(RNG) );                            \
    REDUCE_ENERGY( GAMMA(RNG, E) );                                            \
  FINALIZE_COLLISION(TYPE)

// Convenience factory for binary rate constants.
#define atomic_rate_constant(sigma) do {                                       \
  float dux = pi->ux - pj->ux;                                                 \
  float duy = pi->uy - pj->uy;                                                 \
  float duz = pi->uz - pj->uz;                                                 \
  float du2 = dux*dux + duy*duy + duz*duz;                                     \
  float mu  = (spi->m*spj->m)/(spi->m+spj->m);                                 \
  return ap->Nca02 * (*(sigma))(ap->alpha_m2*mu*du2) * sqrtf(du2);             \
} while(0);

#define rate_factory(ABXX)                                                     \
float                                                                          \
ABXX##_rate_constant( const atomic_properties_t * RESTRICT ap,                 \
                      const species_t           * RESTRICT spi,                \
                      const species_t           * RESTRICT spj,                \
                      const particle_t          * RESTRICT pi,                 \
                      const particle_t          * RESTRICT pj ) {              \
  atomic_rate_constant(ap->sigma_##ABXX);                                      \
}

// Convenience factory for binary collisions. We should _not_ use
// restricted pointers because technically they are aliased in the
// macros.
#define inelastic_collision_factory(ABXX, DELTA)                               \
void                                                                           \
ABXX##_collision( const atomic_properties_t * ap,                              \
                  const species_t           * spi,                             \
                  const species_t           * spj,                             \
                  /**/  particle_t          * pi,                              \
                  /**/  particle_t          * pj,                              \
                  /**/  rng_t               * rng,                             \
                  const int                   type ) {                         \
  INELASTIC_COLLISION(ap, spi, spj, pi, pj, rng, type,                         \
                      ap->delta_##DELTA, ap->gamma_##ABXX );                   \
}

#define elastic_collision_factory(ABXX)                                        \
void                                                                           \
ABXX##_collision( const atomic_properties_t * ap,                              \
                  const species_t           * spi,                             \
                  const species_t           * spj,                             \
                  /**/  particle_t          * pi,                              \
                  /**/  particle_t          * pj,                              \
                  /**/  rng_t               * rng,                             \
                  const int                   type ) {                         \
  ELASTIC_COLLISION(ap, spi, spj, pi, pj, rng, type, ap->delta_##ABXX);        \
}

// Now actually create the wrappers.
rate_factory(pt_el);
rate_factory(pt_cx);
rate_factory(pt_ex);
rate_factory(et_el);
rate_factory(et_ex);
rate_factory(tt_el);
rate_factory(tt_ex);

elastic_collision_factory(pt_el);
elastic_collision_factory(pt_cx);
elastic_collision_factory(et_el);
elastic_collision_factory(tt_el);

// Inelastic collisions use the scattering angle
// from their elastic counterpart
inelastic_collision_factory(pt_ex, pt_el);
inelastic_collision_factory(et_ex, et_el);
inelastic_collision_factory(tt_ex, tt_el);

#undef rate_factory
#undef elastic_collision_factory
#undef inelastic_collision_factory

// Ionization and recombination collisions are special
// and use the chemical collision operator
float
ep_rx_rate_constant( atomic_properties_t  * RESTRICT ap,
                     species_t           ** RESTRICT reactants,
                     particle_t          ** RESTRICT particles ) {
   particle_t * RESTRICT pi = particles[0];
   particle_t * RESTRICT pj = particles[1];
   species_t  * RESTRICT spi = reactants[0];
   species_t  * RESTRICT spj = reactants[1];
   atomic_rate_constant( ap->sigma_ep_rx );
}

void
ep_rx_collision( atomic_properties_t  * RESTRICT ap,
                 species_t           ** RESTRICT reactants,
                 particle_t          ** RESTRICT reactant_particles,
                 species_t           ** RESTRICT products,
                 particle_t          ** RESTRICT product_particles,
                 rng_t                * RESTRICT rng,
                 int                             update ) {

    // In free recombination, momentum conservation tells us exactly what the
    // product velocity should be (assuming photon momentum is negligible).
    particle_t * RESTRICT pi = reactant_particles[0];
    particle_t * RESTRICT pj = reactant_particles[1];
    particle_t * RESTRICT pk = product_particles[0];
    species_t  * RESTRICT spi = reactants[0];
    species_t  * RESTRICT spj = reactants[1];
    species_t  * RESTRICT spk = products[0];

    float mi = spi->m / spk->m;
    float mj = spj->m / spk->m;

    pk->ux = mi * pi->ux + mj * pj->ux ;
    pk->uy = mi * pi->uy + mj * pj->uy ;
    pk->uz = mi * pi->uz + mj * pj->uz ;

}

float
et_in_rate_constant( atomic_properties_t  * RESTRICT ap,
                     species_t           ** RESTRICT reactants,
                     particle_t          ** RESTRICT particles ) {
   particle_t * RESTRICT pi = particles[0];
   particle_t * RESTRICT pj = particles[1];
   species_t  * RESTRICT spi = reactants[0];
   species_t  * RESTRICT spj = reactants[1];
   atomic_rate_constant( ap->sigma_et_in );
}

void
et_in_collision( atomic_properties_t  * RESTRICT ap,
                 species_t           ** RESTRICT reactants,
                 particle_t          ** RESTRICT reactant_particles,
                 species_t           ** RESTRICT products,
                 particle_t          ** RESTRICT product_particles,
                 rng_t                * RESTRICT rng,
                 int                             update ) {

    // FIXME: Designation of electron, target, parent are hardcoded
    //        and assuming the order specified in the public interface
    //        below.

    particle_t * pe1 = reactant_particles[0];
    particle_t * pt  = reactant_particles[1];
    particle_t * pe2 = product_particles[0];
    particle_t * pp  = product_particles[1];

    const species_t * spe1 = reactants[0];
    const species_t * spt  = reactants[1];
    const species_t * spe2 = products[0];
    const species_t * spp  = products[1];

    float me2, mt, dE, gamma;
    int type;

    // Scatter the primary electron off of the target and reduce the energy,
    // but only update the electron if we are told to. In the current
    // implementation, we are assuming that we can write the triple differential
    // cross-section (d^3 sigma(E1)/ dE2 dOmega1 dOmega2) as the seperable product
    // of the singly differential cross section d sigma(E1)/ dE2 and the
    // elastic scattering cross-sections. This is not a good approximation in
    // general. It would be better to fully implement the triply differential
    // cross-section, but data is sparse. Testing this for now.

    if( update & 1 ) type = 1; type += 2;
    BEGIN_COLLISION(ap, spe1, spt, pe1, pt);
      REDUCE_ENERGY(ap->gamma_et_in(rng, E, &dE));
      SCATTER( ap->delta_et_el(rng, dE), 2*M_PI*frand_c0(rng) ); // Scatter pe2

      dE = sqrtf(dE*mu/(E*spe2->m)); // Secondary v/ur
      pe2->ux = cmx - dE*stack[0];   // Initialize secondary ux moving with pt
      pe2->uy = cmy - dE*stack[1];   // Initialize secondary uy moving with pt
      pe2->uz = cmz - dE*stack[2];   // Initialize secondary uz moving with pt

      SCATTER( ap->delta_et_el(rng, E), 2*M_PI*frand_c0(rng) ); // Scatter pe1
    FINALIZE_COLLISION(type);

    // Finally, copy target momentum to parent and reduce it by the secondary
    // momentum. Note that if implemented correctly at the user level,
    // m_p + m_e = m_t in order to ensure mass conservation, but even if this
    // is not true, we always conserve momentum. How this is specifically
    // implemented ie debatble.

    me2 = spe2->m / spp->m;
    mt  = spt->m / spp->m;
    pp->ux = mt * pt->ux - me2 * pe2->ux;
    pp->uy = mt * pt->uy - me2 * pe2->uy;
    pp->uz = mt * pt->uz - me2 * pe2->uz;

}

float
eep_rx_rate_constant( atomic_properties_t  * RESTRICT ap,
                      species_t           ** RESTRICT reactants,
                      particle_t          ** RESTRICT particles ){

  ERROR(("Not implemented."));
  return 0;

}

void
eep_rx_collision( atomic_properties_t  * RESTRICT ap,
                 species_t           ** RESTRICT reactants,
                 particle_t          ** RESTRICT reactant_particles,
                 species_t           ** RESTRICT products,
                 particle_t          ** RESTRICT product_particles,
                 rng_t                * RESTRICT rng,
                 int                             update ) {

  ERROR(("Not implemented."));

}

#undef CMOV
#undef BEGIN_COLLISION
#undef SCATTER
#undef REDUCE_ENERGY
#undef FINALIZE_COLLISION
#undef ELASTIC_COLLISION
#undef INELASTIC_COLLISION
#undef atomic_rate_constant


// Helper routine to decide how to pick the sampling fraction for each process.
double
sampling_fraction( const atomic_properties_t * RESTRICT ap,
                   const species_t           * RESTRICT spi,
                   const species_t           * RESTRICT spj,
                   /**/  float             (*const sigma)(float),
                   const int                            interval,
                   const double                         Emax,
                   const double                         Emin,
                   const int                            nsample){

  double E, K, dE, Kmax = 0;
  const float mu  = (spi->m*spj->m)/(spi->m+spj->m);

  // Find the maxmium sigma v over the energy range using brute force.
  E  = Emin;
  dE = (Emax-Emin)/nsample;
  for(int i=0 ; i<=nsample ; ++i, E += dE){
    K = (*sigma)(ap->alpha_m2*E) * sqrt(E);
    Kmax = K > Kmax ? K : Kmax;
  }
  Kmax *= ap->Nca02 * sqrt(2/mu) ;

  // Now figure out sampling fraction. The number of collisions between species
  // i and j is equal to sample np = sample ni nj, which should be larger than
  // max(wi,wj) ni nj Kmax interval dt / dV. So,
  // sample > Kmax dtinterval_dV * max(wi,wj)

  if( spi->np && spj->np ){
    // Implicit assumption that all particles in each species have the same weight.
    double wi = spi->p[0].w;
    double wj = spj->p[0].w;

    // Add a factor of 2 for safety.
    Kmax *= (wi > wj ? wi : wj) * spi->g->dt * interval / spi->g->dV;
    return Kmax;
  }
  else {
    WARNING(("Cannot determine an appropriate sampling fraction for collisions "
             "between %s and %s because one (or more) species are empty. "
             "Setting to 1\n", spi->name, spj->name));
    return 1;
  }
}

void
checkpt_atomic_properties( const atomic_properties_t * ap ){
  CHECKPT(ap, 1);
  CHECKPT_SYM( ap->sigma_pt_el ); CHECKPT_SYM( ap->delta_pt_el );
  CHECKPT_SYM( ap->sigma_pt_cx ); CHECKPT_SYM( ap->delta_pt_cx );
  CHECKPT_SYM( ap->sigma_et_el ); CHECKPT_SYM( ap->delta_et_el );
  CHECKPT_SYM( ap->sigma_tt_el ); CHECKPT_SYM( ap->delta_tt_el );
  CHECKPT_SYM( ap->sigma_pt_ex ); CHECKPT_SYM( ap->gamma_pt_ex );
  CHECKPT_SYM( ap->sigma_et_ex ); CHECKPT_SYM( ap->gamma_et_ex );
  CHECKPT_SYM( ap->sigma_et_in ); CHECKPT_SYM( ap->gamma_et_in );
  CHECKPT_SYM( ap->sigma_tt_ex ); CHECKPT_SYM( ap->gamma_tt_ex );
  CHECKPT_SYM( ap->sigma_ep_rx );
}

atomic_properties_t *
restore_atomic_properties( void ){
  atomic_properties_t * ap;
  RESTORE( ap );
  RESTORE_SYM( ap->sigma_pt_el ); RESTORE_SYM( ap->delta_pt_el );
  RESTORE_SYM( ap->sigma_pt_cx ); RESTORE_SYM( ap->delta_pt_cx );
  RESTORE_SYM( ap->sigma_et_el ); RESTORE_SYM( ap->delta_et_el );
  RESTORE_SYM( ap->sigma_tt_el ); RESTORE_SYM( ap->delta_tt_el );
  RESTORE_SYM( ap->sigma_pt_ex ); RESTORE_SYM( ap->gamma_pt_ex );
  RESTORE_SYM( ap->sigma_et_ex ); RESTORE_SYM( ap->gamma_et_ex );
  RESTORE_SYM( ap->sigma_et_in ); RESTORE_SYM( ap->gamma_et_in );
  RESTORE_SYM( ap->sigma_tt_ex ); RESTORE_SYM( ap->gamma_tt_ex );
  RESTORE_SYM( ap->sigma_ep_rx );
  return ap;
}

/* Public interface **********************************************************/

#define add_binary_op(SPI,SPJ,ABXX) if(SPI&&SPJ&&ap->sigma_##ABXX)             \
append_collision_op(                                                           \
    binary_collision_model(                                                    \
      #ABXX,                                                                   \
      (binary_rate_constant_func_t) ABXX##_rate_constant,                      \
      (binary_collision_func_t)     ABXX##_collision,                          \
      ap, SPI, SPJ, rp,                                                        \
      sample > 0 ? sample : sampling_fraction(ap,SPI,SPJ,ap->sigma_##ABXX,     \
                                              interval,Emax,0,100),            \
      interval,                                                                \
      mass_action ),                                                           \
    cop_list                                                                   \
  )

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
                                 ) {

  if( !electron || !electron->g || !target || !target->g ||
      !parent || !parent->g || !base || !rp ) ERROR(( "Bad args" ));

  if( sample <= 0 && Emax <= 0 ) ERROR(("Either sample or Emax must be given."));

  atomic_properties_t *ap;
  MALLOC( ap, 1 );
  COPY( ap, base, 1 );

  // Set up properties.
  ap->alpha_m2   = 1/(alpha*alpha);
  ap->Nca02      = Na02 * electron->g->cvac;

  // Register the object.
  REGISTER_OBJECT( ap,
                   checkpt_atomic_properties,
                   restore_atomic_properties, NULL );

  // Create the collision operators.
  if( processes & parent_target_elastic )         add_binary_op(parent,   target, pt_el);
  if( processes & parent_target_inelastic )       add_binary_op(parent,   target, pt_ex);
  if( processes & parent_target_charge_exchange ) add_binary_op(parent,   target, pt_cx);
  if( processes & electron_target_elastic )       add_binary_op(electron, target, et_el);
  if( processes & electron_target_inelastic )     add_binary_op(electron, target, et_ex);
  if( processes & target_target_elastic )         add_binary_op(target,   target, tt_el);
  if( processes & target_target_inelastic )       add_binary_op(target,   target, tt_ex);

  // Non-number preserving operators. These use the chemical collision operator.

  // Ionization.
  if( processes & electron_target_ionization && ap->sigma_et_in ){
    species_t *reactants[2] = { electron, target };
    species_t  *products[2] = { electron, parent };
    int       consumable[2] = {0, 1};
    double s = sample;

    if( sample <= 0 )
      s = sampling_fraction(ap, electron, target, ap->sigma_et_in,
                            interval, Emax, 0, 100);

    collision_op_t * et_in = chemical_collision_model("et_in", ap,
                                                      reactants, consumable, 2,
                                                      products, 2,
                      (chemical_rate_constant_func_t) et_in_rate_constant,
                          (chemical_collision_func_t) et_in_collision,
                                                      rp, fa, s, 1, interval);
    append_collision_op(et_in, cop_list);
  }

  // Free recombination.
  if( processes & electron_parent_recombination && ap->sigma_ep_rx ){
    species_t *reactants[2] = { electron, parent };
    species_t  *products[1] = { target };
    int       consumable[2] = {1, 1};
    double s = sample;

    if( sample <= 0 )
      // Divide by 0 at E=0 !
      s = sampling_fraction(ap, electron, parent, ap->sigma_ep_rx,
                            interval, Emax, 1e-2*Emax, 99);

    collision_op_t * ep_rx = chemical_collision_model("ep_rx", ap,
                                                      reactants, consumable, 2,
                                                      products, 1,
                      (chemical_rate_constant_func_t) ep_rx_rate_constant,
                          (chemical_collision_func_t) ep_rx_collision,
                                                      rp, fa, s, 1, interval);
    append_collision_op(ep_rx, cop_list);
  }

  // Three-body recombination.
  if( processes & three_body_recombination ){
    ERROR(("Three body recombination is not implemented."));
  }

  return *cop_list;
}

#undef add_binary_op
