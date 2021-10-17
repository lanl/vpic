#include "collision.h"

/* Private interface *********************************************************/

typedef struct large_angle_coulomb
{
  float cc, twomu_mi, twomu_mj, Kc;
  float udx, udy, udz, ut;
  float ut2, alpha_Kt2ut4, beta_Kt2ut2, gamma_Kt2;
} large_angle_coulomb_t;

/* See hard_sphere.c for a derivation of these. */

float
large_angle_coulomb_fluid_rate_constant(
    const large_angle_coulomb_t * RESTRICT lac,
    const species_t             * RESTRICT spi,
    const particle_t            * RESTRICT pi )
{
  static const float gamma = (3.*M_PI-8.)/(24.-6*M_PI);
  float urx = pi->ux - lac->udx;
  float ury = pi->uy - lac->udy;
  float urz = pi->uz - lac->udz;
  float ur2 = urx*urx + ury*ury + urz*urz;
  return sqrtf((lac->alpha_Kt2ut4+ur2*(lac->beta_Kt2ut2+ur2*lac->gamma_Kt2))/
               (lac->ut2+ur2*gamma));
}

float
large_angle_coulomb_rate_constant(
    const large_angle_coulomb_t * RESTRICT lac,
    const species_t             * RESTRICT spi,
    const species_t             * RESTRICT spj,
    const particle_t            * RESTRICT pi,
    const particle_t            * RESTRICT pj )
{
  float urx = pi->ux - pj->ux;
  float ury = pi->uy - pj->uy;
  float urz = pi->uz - pj->uz;
  return lac->Kc*sqrtf( urx*urx + ury*ury + urz*urz );
}

/* See hard_sphere.c for a derivation of the basic structure of this.

   The only difference is that, in a large angle Coulomb collision,
   the relationship between the impact parameter and b is:
     b = ( qi qj cot theta/2 ) / ( 4 pi eps0 mu vr0^2 )
   or, given b:
     theta = 2 acot ( ( 4 pi eps0 mu vr0^2 b ) / ( qi qj ) )
           = 2 acot B

   Given theta, we know that, for an elastic collision, a
   computationally efficient expression for the change in relative
   velocity is:
     dvr = vr0 ( cos theta - 1 ) + sin theta cos phi |vr0| T1 +
                                   sin theta sin phi vr0 x T1 +
   where phi is picked randomly on [0,2pi).  Noting that via
   basic trig:
     cos acot B = B / sqrt( 1 + B^2 )
   and:
     sin acot B = 1 / sqrt( 1 + B^2 )
   we have:
     cos theta - 1 = cos 2 acot B - 1
                   = cos^2 acot B - sin^2 acot B - 1
                   = B^2/(1+B^2) - 1/(1+B^2) - 1
                   = ( B^2 - 1 - B^2 - 1 ) / ( B^2 + 1 )
                   = -2 / ( B^2 + 1 )
   and:
     sin theta     = sin 2 acot B
                   = 2 cos acot B sin acot B
                   = 2 B / ( B^2+1 )
   such that:
     dvr = -2 (1 / (B^2+1)) vr0
         +  2 (B / (B^2+1)) cos phi |vr0| T1
         +  2 (B / (B^2+1)) sin phi vr0 x T1

   Normalizing to -2 c:

     dur = (1 / (B^2+1)) ur0
         - (B / (B^2+1)) cos phi |ur0| T1
         - (B / (B^2+1)) sin phi  ur0 x T1 */

#define CMOV(a,b) if(t0<t1) a=b

#define COMPUTE_MOMENTUM_TRANSFER(urx,ury,urz,ax,ay,az,rng)                   \
  do {									      \
    float bcs_bmax, bsn_bmax, b2_bmax2, ur2, ur, tx, ty, tz;                  \
    float t0, t1, t2, stack[3];                                               \
    int d0, d1, d2;                                                           \
                                                                              \
    do {                                                                      \
      bcs_bmax = 2*frand_c0(rng) - 1;                                         \
      bsn_bmax = 2*frand_c0(rng) - 1;                                         \
      b2_bmax2 = bcs_bmax*bcs_bmax + bsn_bmax*bsn_bmax;                       \
    } while( b2_bmax2>=1 );                                                   \
                                                                              \
    /* There are lots of ways to formulate T vector formation    */           \
    /* This has no branches (but uses L1 heavily)                */           \
                                                                              \
    t0 = urx*urx;      d0=0;       d1=1;       d2=2;       t1=t0;  ur2  = t0; \
    t0 = ury*ury; CMOV(d0,1); CMOV(d1,2); CMOV(d2,0); CMOV(t1,t0); ur2 += t0; \
    t0 = urz*urz; CMOV(d0,2); CMOV(d1,0); CMOV(d2,1);              ur2 += t0; \
    ur = sqrtf( ur2 );                                                        \
                                                                              \
    stack[0] = urx;                                                           \
    stack[1] = ury;                                                           \
    stack[2] = urz;                                                           \
    t1  = stack[d1];                                                          \
    t2  = stack[d2];                                                          \
    t0  = 1 / sqrtf( t1*t1 + t2*t2 + FLT_MIN );                               \
    stack[d0] =  0;                                                           \
    stack[d1] =  t0*t2;                                                       \
    stack[d2] = -t0*t1;                                                       \
    tx = stack[0];                                                            \
    ty = stack[1];                                                            \
    tz = stack[2];                                                            \
                                                                              \
    t2 = lac->cc;                /* 4 pi eps0 mu c^2 bmax / (qi qj) */        \
    t1 = t2 * ur2;               /*  B (bmax / b)                   */        \
    t0 = 1/(1+(t1*t1)*b2_bmax2); /*  1 / ( B^2 + 1 )                */        \
    t2 = t0*t1;                  /* (B / ( B^2 + 1 ))(bmax / b)     */        \
    t1 = t2*bcs_bmax*ur;         /* (B / (B^2+1)) cos phi |ur0|     */        \
    t2 = t2*bsn_bmax;            /* (B / (B^2+1)) sin phi           */        \
                                                                              \
    ax = (t0*urx - t1*tx) - t2*( ury*tz - urz*ty );                           \
    ay = (t0*ury - t1*ty) - t2*( urz*tx - urx*tz );                           \
    az = (t0*urz - t1*tz) - t2*( urx*ty - ury*tx );                           \
  } while(0)

/* It would be nice to preserve redundant rate constant
   computations */

void
large_angle_coulomb_fluid_collision(
    const large_angle_coulomb_t * RESTRICT lac,
    const species_t             * RESTRICT spi,
    /**/  particle_t            * RESTRICT pi,
    /**/  rng_t                 * RESTRICT rng )
{
  float urx, ury, urz, ax, ay, az, w;

  urx = pi->ux - lac->udx;
  ury = pi->uy - lac->udy;
  urz = pi->uz - lac->udz;

  w = lac->ut;
  if( w ) {
    urx -= w*frandn(rng);
    ury -= w*frandn(rng);
    urz -= w*frandn(rng);
  }

  COMPUTE_MOMENTUM_TRANSFER(urx,ury,urz,ax,ay,az,rng);

  w = lac->twomu_mi;
  pi->ux -= w*ax;
  pi->uy -= w*ay;
  pi->uz -= w*az;
}

void
large_angle_coulomb_collision(
    const large_angle_coulomb_t * RESTRICT lac,
    const species_t             * RESTRICT spi,
    const species_t             * RESTRICT spj,
    /**/  particle_t            * RESTRICT pi,
    /**/  particle_t            * RESTRICT pj,
    /**/  rng_t                 * RESTRICT rng,
    const int                              type )
{
  float urx, ury, urz, ax, ay, az, w;

  urx = pi->ux - pj->ux;
  ury = pi->uy - pj->uy;
  urz = pi->uz - pj->uz;

  COMPUTE_MOMENTUM_TRANSFER(urx,ury,urz,ax,ay,az,rng);

  if( type & 1 ) {
    w = lac->twomu_mi;
    pi->ux -= w*ax;
    pi->uy -= w*ay;
    pi->uz -= w*az;
  }

  if( type & 2 ) {
    w = lac->twomu_mj;
    pj->ux += w*ax;
    pj->uy += w*ay;
    pj->uz += w*az;
  }
}

#undef CMOV

void
checkpt_large_angle_coulomb( const large_angle_coulomb_t * lac )
{
  CHECKPT( lac, 1 );
}

large_angle_coulomb_t *
restore_large_angle_coulomb( void )
{
  large_angle_coulomb_t * lac;
  RESTORE( lac );
  return lac;
}

/* Public interface **********************************************************/

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
    const int interval )        /* How often to apply this */
{
  large_angle_coulomb_t * lac;

  if( n0<0 || kT0<0 || !q0 || m0<=0 || !sp || !sp->q || sp->m<=0 || bmax<0 )
    ERROR(( "Bad args" ));

  MALLOC( lac, 1 );

  lac->cc       = (4.*M_PI*sp->g->eps0*sp->m*m0*sp->g->cvac*sp->g->cvac*bmax) /
                  ((sp->m + m0)*sp->q*q0);
  lac->twomu_mi = 2.*m0    / (sp->m + m0);
  lac->twomu_mj = 2.*sp->m / (sp->m + m0);
  lac->Kc       = M_PI*bmax*bmax*sp->g->cvac;

  lac->udx      = vdx / sp->g->cvac;
  lac->udy      = vdy / sp->g->cvac;
  lac->udz      = vdz / sp->g->cvac;
  lac->ut       = sqrtf(kT0/(m0*sp->g->cvac*sp->g->cvac));

  lac->ut2           = kT0/(m0*sp->g->cvac*sp->g->cvac);
  lac->alpha_Kt2ut4  = (8./M_PI)*(lac->Kc*n0*lac->Kc*n0)*(lac->ut2*lac->ut2);
  lac->beta_Kt2ut2   = (4./(12.-3.*M_PI))*(lac->Kc*n0*lac->Kc*n0)*lac->ut2;
  lac->gamma_Kt2     = ((3.*M_PI-8.)/(24.-6*M_PI))*(lac->Kc*n0*lac->Kc*n0);
  lac->ut2          += FLT_MIN;

  REGISTER_OBJECT( lac,
                   checkpt_large_angle_coulomb,
                   restore_large_angle_coulomb, NULL );
  return unary_collision_model( name,
           (unary_rate_constant_func_t)large_angle_coulomb_fluid_rate_constant,
           (unary_collision_func_t)    large_angle_coulomb_fluid_collision,
                                lac, sp, rp, interval );
}

collision_op_t *
large_angle_coulomb( const char * RESTRICT name, /* Model name */
                     species_t * RESTRICT spi,   /* Species-i */
                     species_t * RESTRICT spj,   /* Species-j */
                     const float bmax,           /* Impact parameter cutoff */
                     rng_pool_t * RESTRICT rp,   /* Entropy pool */
                     const double sample,        /* Sampling density */
                     const int interval,         /* How often to apply this */
                     const int strategy )        /* Sampling strategy */
{
  large_angle_coulomb_t * lac;

  if( !spi || !spi->q || spi->m<=0 ||
      !spj || !spj->q || spj->m<=0 || spi->g!=spj->g ) ERROR(( "Bad args" ));

  MALLOC( lac, 1 );
  CLEAR(  lac, 1 );
  lac->cc       = (4.*M_PI*spi->g->eps0*spi->m*spj->m*spi->g->cvac*spi->g->cvac*bmax) /
                  ((spi->m + spj->m)*spi->q*spj->q);
  lac->twomu_mi = 2*spj->m/(spi->m+spj->m);
  lac->twomu_mj = 2*spi->m/(spi->m+spj->m);
  lac->Kc       = spi->g->cvac*M_PI*bmax*bmax;

  REGISTER_OBJECT( lac,
                   checkpt_large_angle_coulomb,
                   restore_large_angle_coulomb, NULL );
  return binary_collision_model( name,
                (binary_rate_constant_func_t)large_angle_coulomb_rate_constant,
                (binary_collision_func_t)    large_angle_coulomb_collision,
                                 lac, spi, spj, rp, sample, interval, strategy );
}
