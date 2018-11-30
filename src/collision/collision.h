#ifndef _collision_h_
#define _collision_h_

/* Note that it is possible to generalize these interfaces to
   accomdate collisional processes involving an arbitrary number of
   bodies (e.g. 3-body recombination processes). */

#include "../species_advance/species_advance.h"

#include <math.h>
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

struct collision_op;
typedef struct collision_op collision_op_t;

BEGIN_C_DECLS

/* In collision.c */

int
num_collision_op( const collision_op_t * RESTRICT cop_list );

void
apply_collision_op_list( collision_op_t * RESTRICT cop_list );

void
delete_collision_op_list( collision_op_t * RESTRICT cop_list );

collision_op_t *
append_collision_op( collision_op_t * cop,
                     collision_op_t ** cop_list );

/* In takizuka_abe.c */

/* The Takizuka-Abe collision model is based on Takizuka and Abe, JCP 1977
   and efficiently models small-angle Coulomb scattering by randomly pairing
   particles. On average, each particle is scattered once per call. The model
   is fully defined by a single parameter, the base collision frequency nu0.
   In SI units, nu0 is defined by

   nu0 = log(Lambda) / 8 pi sqrt(2) c^3 eps0^2

   where log(Lambda) is the Coulomb logarithm. For a thermal species with
   temperature T, normalized mass m and charge q, the self-scattering momentum
   transfer rate (i.e., "the" collision rate) is related to nu0 via

   nu_s = 4 (mc^2 / T)^3 nu0 / 3 sqrt(pi)
*/

collision_op_t *
takizuka_abe( const char       * RESTRICT name,
              /**/  species_t  * RESTRICT spi,
              /**/  species_t  * RESTRICT spj,
              /**/  rng_pool_t * RESTRICT rp,
              const double                nu0,
              const int                   interval );

/* In langevin.c */

/* The most basic collision model (but implemented with numerical
   sophistication).  nu is the collision frequency of the particles
   with some unresolved stationary large thermal bath.  kT is the
   thermal bath temperature.  This method is stable (e.g. if you set
   nu very large, it is equivalent to resampling your particle
   normal momenta from a normal distribution with
   uth = sqrt(kT/mc)) every time this operator is applied (in MD,
   this is called an Anderson thermostat).  This method is only
   intended to be used when the temperature is non-relativistic.

   For the pedants, this operator applies exactly (in exact
   arithmetic), the stochastic operator:
     du = -nu dt + sqrt((2kT)/(mc)) dW
   for the finite duration:
     sp->g->dt * interval
   every interval time step to all the particle momenta in the
   species.  Above, dW is a basic Weiner process. */

collision_op_t *
langevin( float                 kT,
          float                 nu,
          species_t  * RESTRICT sp,
          rng_pool_t * RESTRICT rp,
          int                   interval );

/* In unary.c */

/* A unary_rate_constant_func_t returns the lab-frame rate constant for
   collisions between a monochromatic beam of particles (with rest mass
   sp->m) and momentum p->u{xyz} (normalized to sp->m sp->g->cvac where
   cvac is the speed of light in vacuum) against some background whose
   properties are determined by the specific collision model.

   The returned value has units of FREQUENCY.

   In the case of collisions with a static background of density
   n_background, the rate constant is:

     vi sigma(vi) n_background

   where vi is cvac |ui| / gamma_i is the physical velocity of particle
   vi in the lab frame and sigma is the model specific collision
   cross section for collisions of particles of that velocity with the
   background.

   The basic control flow for a unary_collision_func_t should be:

     void
     my_unary_rate_constant( my_collision_model_params_t * RESTRICT params,
                             species_t * RESTRICT sp,
                             particle_t * RESTRICT ALIGNED(32) p ) {
       return vi sigma(vi) n_background;
     } */

typedef float
(*unary_rate_constant_func_t)( /**/  void       * RESTRICT params,
                               const species_t  * RESTRICT sp,
                               const particle_t * RESTRICT ALIGNED(32) p );

/* A unary_collision_func_t implements the microscopic physics of a
   collision between a particle and some background whose properties
   are determined by the specific collision model.

   The basic control flow for a unary_collision_func_t should be:

     void
     my_unary_collide( my_collision_model_params_t * RESTRICT params,
                       const species_t * RESTRICT sp,
                       particle_t * RESTRICT ALIGNED(32) p,
                       rng_t * RESTRICT rng ) {
       pi->u{xyz} = final momentum, ui{xyz} of particle i
                    given fluid background and initial momentum
                    pi->u{xyz}
     } */

typedef void
(*unary_collision_func_t)( /**/  void       * RESTRICT params,
                           const species_t  * RESTRICT sp,
                           /**/  particle_t * RESTRICT ALIGNED(32) p,
                           /**/  rng_t      * RESTRICT rng );

/* Declare a unary collision model with the given microscopic physics.
   params must be a registered object or NULL.  Every particle is
   tested for collision on every "interval" timesteps.  */

collision_op_t *
unary_collision_model( const char       * RESTRICT name,
                       unary_rate_constant_func_t  rate_constant,
                       unary_collision_func_t      collision,
                       /**/  void       * RESTRICT params,
                       /**/  species_t  * RESTRICT sp,
                       /**/  rng_pool_t * RESTRICT rp,
                       int                         interval );

/* In binary.c */

/* A binary_rate_constant_func_t returns the lab-frame rate constant
   for the collisions between a monochromatic beam of species i
   _physical_ particles (with mass spi->m) and momentum pi->u{xyz}
   (normalized to spi->m spi->g->cvac where cvac is the speed of light
   in vacuum) and a monochromatic beam of species j _physical_
   particles (with mass spj->m) with momentum pj->u{xyz}
   (normalized to spj->m spj->g->cvac).

   The returned value has units of VOLUME / TIME.

   For simple non-relativistic collisions, the rate constant, K, is
   related to the total collision cross section by:

     K = vr sigma(vr)

   where vr = cvac | ui - uj | and sigma( vr ) is the collision cross
   section for particles of that relative velocity.

   For relativistic collisions, this needs to be modified to:

     K = vr sigma( vr ) [ 1 - vi.vj / c^2 ]

   where vr = sqrt( |vi-vj|^2 - |vi x vj|^2/c^2 ) is the relative
   particle velocity in a frame in which one particle is at rest, vi
   = c ui / gamma_i is particle i's lab frame velocity, gamma_i =
   sqrt(1+ui^2) is particle i's lab frame relativistic factor and
   similarly for vj and gamma_j.  [CITE: Peano et al, ARXIV
   pre-print, 2009].

   This is both inefficient and numerically unsafe method to compute
   K in this regime (it misbehaves badly in finite precision for
   relativistic particles ... which is the whole point of doing this
   correct relativistically).  Relativistically, K is better
   computed via:

     s  = gamma_i gamma_j - ui.uj - 1
     vr = cvac sqrt{ s / [ s + 1/(2+s) ] }
     K  = vr sigma( vr ) (1+s) / (gamma_i gamma_j)

   which, provided s is computed with care, has _no_ catastropic
   cancellations, no near singular divisions, behaves well for
   non-relativistic particles and behaves wells for
   ultra-relativistic particles.

   For relativity afficinados, note that s is related to the Lorentz
   boost factor and has the manifestly covariant expression s =
   Ui.Uj - 1 where Ui = (gamma_i,ui) and Uj = (gamma_j,uj) are the
   normalized 4-momenta of particle i and particle j respectively
   and the Minkowski 4-dot product has a +--- signature.

   The basic control flow for a unary_collision_func_t should be:

     void
     my_binary_rate_constant( my_collision_model_params_t * RESTRICT params,
                              const species_t * RESTRICT spi,
                              const species_t * RESTRICT spj,
                              const particle_t * RESTRICT ALIGNED(32) pi,
                              const particle_t * RESTRICT ALIGNED(32) pj ) {
       return vr sigma(vr);
     } */

typedef float
(*binary_rate_constant_func_t)( /**/  void       * RESTRICT params,
                                const species_t  * RESTRICT spi,
                                const species_t  * RESTRICT spj,
                                const particle_t * RESTRICT ALIGNED(32) pi,
                                const particle_t * RESTRICT ALIGNED(32) pj );

/* A binary_collision_func_t implements the microscopic physics of a
   collision between two particles, pi and pj.  The basic control
   flow for a binary_collision_func_t should be:

     void
     my_collide_func( my_collision_model_params_t * RESTRICT params,
                      const species_t * RESTRICT spi,
                      const species_t * RESTRICT spj,
                      particle_t * RESTRICT ALIGNED(32) pi,
                      particle_t * RESTRICT ALIGNED(32) pj,
                      rng_t * RESTRICT rng,
                      int type ) {

       ... compute the final normalized momenta, ui{xyz} and uj{xyz}
       ... between two colliding _physical_ particles, one from
       ... species spi with initial normalized momenta pi->u{xyz} and
       ... the other from spj with initial momentum pj->u{xyz}.

       if( type & 1 ) pi->u{xyz} = ui{xyz};
       if( type & 2 ) pj->u{xyz} = uj{xyz};
     } */

typedef void
(*binary_collision_func_t)( /**/  void       * RESTRICT params,
                            const species_t  * RESTRICT spi,
                            const species_t  * RESTRICT spj,
                            /**/  particle_t * RESTRICT ALIGNED(32) pi,
                            /**/  particle_t * RESTRICT ALIGNED(32) pj,
                            /**/  rng_t      * RESTRICT rng,
                            int type );

/* Declare a binary collision model with the given microscopic physics.
   params must be a registered object or NULL.  A particle in a species
   will be tested for collision on average at least "sample" times every
   "interval" timesteps.  */

collision_op_t *
binary_collision_model( const char       * RESTRICT name,
                        binary_rate_constant_func_t rate_constant,
                        binary_collision_func_t     collision,
                        /**/  void       * RESTRICT params,
                        /**/  species_t  * RESTRICT spi,
                        /**/  species_t  * RESTRICT spj,
                        /**/  rng_pool_t * RESTRICT rp,
                        double                      sample,
                        int                         interval );

/* In hard_sphere.c */

/* Based on unary_collision_model */

collision_op_t *
hard_sphere_fluid( const char * RESTRICT name, /* Model name */
                   const float n0,             /* Fluid density (#/VOLUME) */
                   const float v0x,            /* Fluid x-drift (VELOCITY) */
                   const float v0y,            /* Fluid y-drift (VELOCITY) */
                   const float v0z,            /* Fluid z-drift (VELOCITY) */
                   const float kT0,            /* Fluid temperature (ENERGY) */
                   const float m0,             /* Fluid p. mass (MASS) */
                   const float r0,             /* Fluid p. radius (LENGTH) */
                   species_t * RESTRICT sp,    /* Species */
                   const float rsp,            /* Species p. radius (LENGTH) */
                   rng_pool_t * RESTRICT rp,   /* Entropy pool */
                   const int interval );       /* How often to apply this */

/* Based on binary_collision_model */

collision_op_t *
hard_sphere( const char * RESTRICT name, /* Model name */
             species_t * RESTRICT spi,   /* Species-i */
             const float ri,             /* Species-i p. radius (LENGTH) */
             species_t * RESTRICT spj,   /* Species-j */
             const float rj,             /* Species-j p. radius (LENGTH) */
             rng_pool_t * RESTRICT rp,   /* Entropy pool */
             const double sample,        /* Sampling density */
             const int interval );       /* How often to apply this */

/* In large_angle_coulomb.c */

/* Based on unary_collision_model */

collision_op_t *
large_angle_coulomb_fluid(
    const char * RESTRICT name, /* Model name */
    const float n0,             /* Fluid density (#/VOLUME) */
    const float vdx,            /* Fluid x-drift (VELOCITY) */
    const float vdy,            /* Fluid y-drift (VELOCITY) */
    const float vdz,            /* Fluid z-drift (VELOCITY) */
    const float kT0,            /* Fluid temperature (ENERGY) */
    const float q0,             /* Fluid particle charge (CHARGE) */
    const float m0,             /* Fluid particle mass (MASS) */
    species_t * RESTRICT sp,    /* Species */
    const float bmax,           /* Impact parameter cutoff */
    rng_pool_t * RESTRICT rp,   /* Entropy pool */
    const int interval );       /* How often to apply this */

/* Based on binary_collision_model */

collision_op_t *
large_angle_coulomb( const char * RESTRICT name, /* Model name */
                     species_t * RESTRICT spi,   /* Species-i */
                     species_t * RESTRICT spj,   /* Species-j */
                     const float bmax,           /* Impact parameter cutoff */
                     rng_pool_t * RESTRICT rp,   /* Entropy pool */
                     const double sample,        /* Sampling density */
                     const int interval );       /* How often to apply this */

END_C_DECLS

#endif /* _collision_h_ */
