#include "collision.h"

/* Private interface *********************************************************/

typedef struct hard_sphere {
  float twomu_mi, twomu_mj, Kc;
  float udx, udy, udz, ut;
  float ut2, alpha_Kt2ut4, beta_Kt2ut2, gamma_Kt2;
} hard_sphere_t;

/* The rate constant for particle i to scatter off particle j is:
     K = pi R^2 |vi-vj|

   where R = ri+rj is the sum of the two particle's radii.  For the
   unary model, j is conceptually drawn from a drifting Maxwellian
   (with density nj, drift vdj and temperature kTj).  The rate
   constant for a particle i to collide with beam particles is:
     K = integral pi R^2 |vi-vj| { nj / [ vtj^3 (2 pi)^(3/2) ] }
                  exp( -(vj-vdj)^2 / (2 vtj^2) ) dvj
   where vtj = sqrt( kTj / mj ).  This integral can be rewritten in
   the frame of the drift as:
     K = { [ nj pi R^2 ] / [ vtj^3 (2 pi)^(3/2) ] }
       integral |vr-vj| exp( -vj^2 / (2 vtj^2) ) dvj
   where vr = vi - vdj.  Split vj into components perpendicular, vjP,
   and parallel, vjp, to vr.  This yields:
     K = { [ nj pi R^2 ] / [ vtj^3 (2 pi)^(3/2) ] }
       integral ( vr^2 - 2|vr|vjp + vjp^2 + vjP^2 )^(1/2)
                exp( -vjp^2 / ( 2 vtj^2 ) ) 
                exp( -vjP^2 / ( 2 vtj^2 ) )
                2 pi vjP dvjP dvjp
   where the vjP integration is from 0 to infinity and the vjp is
   -infinity to infinity.  Let z = vjP^2, noting that dz = 2 vjP dvjP,
   we have:
     K = { [ nj pi^2 R^2 ] / [ vtj^3 (2 pi)^(3/2) ] }
       integral ( vr^2 - 2|vr|vjp + vjp^2 + z )^(1/2)
                exp( -vjp^2 / ( 2 vtj^2 ) ) 
                exp( -z / ( 2 vtj^2 ) ) dz dvjp

   Performing the z integral, with the help of Mathematica, we get:
     K = { [ nj pi^2 R^2 ] / [ vtj^3 (2 pi)^(3/2) ] }
       integral [ 2 vtj^2 ||vr|-vjp| exp( -vjp^2 / ( 2 vtj^2 ) ) +
                 (2 pi)^(1/2) vtj^3 x
                 exp( [(|vr|-vjp)^2-vjp^2] / ( 2 vtj^2 ) ) x
                 erfc( ||vr|-vjp| / ( 2^(1/2) vtj ) ) ]
   With the help of Mathematica again, performing the second integral
   and simplifying we get:

     K = nj pi R^2 [ ( |vr| + vtj^2 / |vr| ) x
                     erf( |vr| / [ 2^(1/2) vtj ] ) +
                     vtj (2/pi)^(1/2) exp( - vr^2 / ( 2 vtj^2 ) ) ]

   Sanity checking, in the limit vtj -> 0, we have:
     K = nj pi R^2 |vr|
   For |vr|>>vt, dropping all terms 3rd order and higher in vtj/|vr|:

     K = nj pi R^2 |vr| [ 1 + vtj^2/|vr|^2 ]

   These make sense.  For |vr|<<vt, dropping all terms 3rd order and
   higher in |vr|/vtj:

     K = nj pi R^2 vtj [ (8/pi)^(1/2) + [2/(9 pi)]^(1/2) |vr|^2/vtj^2 ]

   Noting (8/pi)^(1/2) vtj is the mean particle velocity magnitude,
   this also makes sense.

   The raw expression above is a pain to compute numerically.  But we
   have 3rd order accurate expansions for both large and small |vr|
   which we can use to design a very efficient function.  We first
   note that K can be written as:

     K = ( ni pi^2 R^2 vtj ) f( |vr|/vtj )

   Noting that f( x ) = sqrt( 8/pi + x^2 ) matches the zeorth order
   behavior for x=0 and x->infinity (indeed, this simple function
   is accurate to about 3% over the whole domain, we consider a
   more general functional form:

     f(x) = sqrt( ( a + b x^2 + c x^4 ) / ( 1 + d x^2 ) )

   To reproduce the the zeroth order behavior at the order and
   asymptotically, obviously, we must have a = 8/pi and c=d.  This
   leaves two degrees of freedom (b and c).  Expanding f(x) about the
   origin and infinity and solving the system of equations that comes
   about when matching this approximation's expansions to the above
   expansions, we get:

     a = 8/pi
     b = 4/(12 - 3 pi)
     c = (3 pi - 8)/(24 - 6 pi)
     d = c

   This simple expression is better than 0.3% accurate over the whole
   domain.  (Note that the above is not strictly a Pade approximant
   for K or for some function of K as we matched our approximation at
   two distinct points.  If we wanted to be tedious, we could match
   higher order terms (we can basically get to floating point single
   precision with a two extra terms in the above.) */

float
hard_sphere_fluid_rate_constant( const hard_sphere_t * RESTRICT hs,
                                 const species_t     * RESTRICT spi,
                                 const particle_t    * RESTRICT pi ) {
  static const float gamma = (3.*M_PI-8.)/(24.-6*M_PI);
  float urx = pi->ux - hs->udx;
  float ury = pi->uy - hs->udy;
  float urz = pi->uz - hs->udz;
  float ur2 = urx*urx + ury*ury + urz*urz;
  return sqrtf((hs->alpha_Kt2ut4+ur2*(hs->beta_Kt2ut2+ur2*hs->gamma_Kt2))/
               (hs->ut2+ur2*gamma));
}

/* The particle-particle case is much easier theoretically. */

float
hard_sphere_rate_constant( const hard_sphere_t * RESTRICT hs,
                           const species_t     * RESTRICT spi,
                           const species_t     * RESTRICT spj,
                           const particle_t    * RESTRICT pi,
                           const particle_t    * RESTRICT pj ) {
  float urx = pi->ux - pj->ux;
  float ury = pi->uy - pj->uy;
  float urz = pi->uz - pj->uz;
  return hs->Kc*sqrtf( urx*urx + ury*ury + urz*urz );
}

/* Conservation of momentum, pi0+pj0 = pi1+pj1, implies that the
   center of mass velocity vcm = (pi+pj)/(mi+mj) is unchanged by the
   collision.  In an elastic collision, conservation of energy further
   implies that magnitude of the relative velocity is conserved.
   Namely, with dvi0 = vi0-vcm, dvj0 = vj0-vcm:
      0.5 mi vi0^2 + 0.5 mj vj0^2 = 0.5 mi vi1^2 + 0.5 mj vj1^2
   => 0.5 mi (dvi0+vcm)^2 + 0.5 mj (dvj0+vcm)^2 =
        0.5 mi (dvi0+vcm)^2 + 0.5 mj (dvj0+vcm)^2
   => 0.5 ( mi+mj ) vcm^2 + 0.5 mi dvi0^2 + 0.5 mj dvj0^2 +
        vcm.( mi dvi0 + mj dvj0 ) =
      0.5 ( mi+mj ) vcm^2 + 0.5 mi dvi1^2 + 0.5 mj dvj1^2 +
        vcm.( mi dvi1 + mj dvj1 )
   But mi dvi0 + mj dvj0 = mi vi0 + mj vj0 - mi vcm - mj vcm = 0
   by conservation of momentum (and similarly for
   mi dvi1 + mj dvj1), leaving:
      0.5 mi dvi0^2 + 0.5 mj dvj0^2 = 0.5 mi dvi1^2 + 0.5 mj dvj1^2
   Noting that:
      dvi0 = vi0-vcm
           = vi0 - ( mi vi0 + mj vj0 ) / (mi+mj)
           = mj (vi0 - vj0 ) / (mi+mj)
           = mu vr0 / mi
   where mu = mi mj / (mi+mj) is the reduced mass and vr0 = vi-vj is
   the relative velocity (and similarly for dvj0, dvi1, dvj1), we
   have:
      0.5 mu^2 (1/mi + 1/mj) vr0^2 = 0.5 mu^2 (1/mi + 1/mj) vr1^2
   or:
      0.5 mu vr0^2 = 0.5 mu vr1^2
   => |vr0| = |vr1|

   In a hard sphere collision, if two particles are known to have
   collided (but we know nothing about the details of the collisions),
   we know that in the reduced mass system, the reduced mass particle
   was incident transversely on the origin particle somewhere in the
   circle from [0,R).  By picking a point uniformly in the unit
   circle, the radius of this point has the same distribution as the
   impact parameter (normalized by R).  Further, the angle of this
   point is uniform on [0,2pi) and uncorrelated with the radius and
   can thus can be used to specify the transverse scattering angle.

   Given this normalized impact parameter, b/R, the collision angle
   theta in the reduced mass system is given by:
     sin( pi/2 - theta/2 ) = b/R
   or:
     cos( theta/2 ) = b/R
   or:
     theta = 2 acos b/R
   This, and the above, imply the final relative velocity is given by:
     vr1 = |vr0| ( cos theta || + sin theta T )
   where || is a unit vector parallel to vr and T is one of the unit
   vectors perpendicular to vr (picked at random).  The change in
   relative velocity is:
     dvr = vr0 (cos theta-1) + |vr0| sin theta T
   Apply trig identities to the above, we find that:
     cos theta - 1 = cos 2 acos b/R - 1
                   = cos^2 acos b/R - sin^2 acos b/R - 1
                   = (b/R)^2 - ( 1 - (b/R)^2 ) - 1
                   = -2 [ 1 - (b/R)^2 ]
     sin theta     = sin 2 acos b/R
                   = 2 cos acos b/R sin acos b/R
                   = 2 (b/R) [ 1 - (b/R)^2 ]^(1/2)
   Yielding:

     dvr = -2 [ 1 - (b/R)^2 ] vr0 + |vr0| 2 (b/R) [ 1 - (b/R)^2 ]^(1/2) T
     
   To convert the scattering angle above into T, we need to determine
   to two vectors perpendicular to the incident velocity, vr0.  Given
   the desire to not introduce any preferred directions in the
   calculation, this is tricky as the only directional cue we have is
   the vr0.  Picking a reference vector in the direction of the axis
   aligned unit vector in which vr0 is most deficient (this unit
   vector by construction is well separated in angle from vr0), we can
   construct a perpendicular vector T1 safely and without bias by
   forming the cross product of this and normalizing (mathematically,
   this reduces to zeroing out the smallest component of vr0, rotating
   the other two components 90 degrees and normalizing).  T2 can be
   formed as the cross product of VR and T1 (and normalizing).  This
   given:

     dvr = -2 [ 1 - (b/R)^2 ] vr0 + 
       |vr0| 2 (b/R) [ 1 - (b/R)^2 ]^(1/2) ( cs T1 + sn T2 )

   Noting that (b/R) cs and (b/R) sn are just the coordinates of the
   point picked in the unit circle originally and that
   T2 = vr0 x T1 / |vr0|, we can simplify the above:

     dvr = -2 [ 1 - (b/R)^2 ] vr0 + 
       2 [ 1 - (b/R)^2 ]^(1/2) [ |vr0| bcs_R T1 + bsn_R ( vr0 x T1 ) ]

   Requiring conservation of momentum (which implies
   mi vi0 + mj vi0 = mi (vi0 + wi dvr) + mj (vj0 - wj dvr) and
   that the relative velocities change as given above (which implies
   (vi0 + wi dvr) - (vj0 - wj dvr) = (vi0-vj0) + dvr, implies:

      mi wi - mj wj = 0
         wi +    wj = 1

   or wi = mj / (mi+mj) = mu / mi and similarly for wj.

   This gives, normalizing to -2 c:

     a = [ 1 - b2_R2 ] ur0 -
         [ 1 - b2_R2 ]^(1/2) bcs_R  |ur0|  T1   -
         [ 1 - b2_R2 ]^(1/2) bsn_R ( ur0 x T1 )

   and:

     ui1 = ui0 - (2 mu/mi) a
     uj1 = uj0 + (2 mu/mj) a

   of the below.  (In short, its fast and correct.) */

#define CMOV(a,b) if(t0<t1) a=b

#define COMPUTE_MOMENTUM_TRANSFER(urx,ury,urz,ax,ay,az,rng) do {        \
    float bcs_R, bsn_R, b2_R2, ur, tx, ty, tz, t0, t1, t2, stack[3];    \
    int d0, d1, d2;                                                     \
                                                                        \
    do {                                                                \
      bcs_R = 2*frand_c0(rng) - 1;                                      \
      bsn_R = 2*frand_c0(rng) - 1;                                      \
      b2_R2 = bcs_R*bcs_R + bsn_R*bsn_R;                                \
    } while( b2_R2>=1 );                                                \
                                                                        \
    /* There are lots of ways to formulate T vector formation    */     \
    /* This has no branches (but uses L1 heavily)                */     \
                                                                        \
    t0 = urx*urx;      d0=0;       d1=1;       d2=2;       t1=t0;  ur  = t0; \
    t0 = ury*ury; CMOV(d0,1); CMOV(d1,2); CMOV(d2,0); CMOV(t1,t0); ur += t0; \
    t0 = urz*urz; CMOV(d0,2); CMOV(d1,0); CMOV(d2,1);              ur += t0; \
    ur = sqrtf( ur );                                                   \
                                                                        \
    stack[0] = urx;                                                     \
    stack[1] = ury;                                                     \
    stack[2] = urz;                                                     \
    t1  = stack[d1];                                                    \
    t2  = stack[d2];                                                    \
    t0  = 1 / sqrtf( t1*t1 + t2*t2 + FLT_MIN );                         \
    stack[d0] =  0;                                                     \
    stack[d1] =  t0*t2;                                                 \
    stack[d2] = -t0*t1;                                                 \
    tx = stack[0];                                                      \
    ty = stack[1];                                                      \
    tz = stack[2];                                                      \
                                                                        \
    t0  = 1 - b2_R2;                                                    \
    t2  = sqrtf( t0 );                                                  \
    t1  = t2*bcs_R*ur;                                                  \
    t2 *= bsn_R;                                                        \
                                                                        \
    ax = (t0*urx - t1*tx) - t2*( ury*tz - urz*ty );                     \
    ay = (t0*ury - t1*ty) - t2*( urz*tx - urx*tz );                     \
    az = (t0*urz - t1*tz) - t2*( urx*ty - ury*tx );                     \
  } while(0)

/* It would be nice to preserve redundant rate constant
   computations */

void
hard_sphere_fluid_collision( const hard_sphere_t * RESTRICT hs,
                             const species_t     * RESTRICT spi,
                             /**/  particle_t    * RESTRICT pi,
                             /**/  rng_t         * RESTRICT rng ) {
  float urx, ury, urz, ax, ay, az, w;

  urx = pi->ux - hs->udx;
  ury = pi->uy - hs->udy;
  urz = pi->uz - hs->udz;

  w = hs->ut;
  if( w ) {
    urx -= w*frandn(rng);
    ury -= w*frandn(rng);
    urz -= w*frandn(rng);
  }
  
  COMPUTE_MOMENTUM_TRANSFER(urx,urz,urz,ax,ay,az,rng);

  w = hs->twomu_mi;
  pi->ux -= w*ax;
  pi->uy -= w*ay;
  pi->uz -= w*az;
}

void
hard_sphere_collision( const hard_sphere_t * RESTRICT hs,
                       const species_t     * RESTRICT spi,
                       const species_t     * RESTRICT spj,
                       /**/  particle_t    * RESTRICT pi,
                       /**/  particle_t    * RESTRICT pj,
                       /**/  rng_t         * RESTRICT rng,
                       const int                      type ) {
  float urx, ury, urz, ax, ay, az, w;

  urx = pi->ux - pj->ux;
  ury = pi->uy - pj->uy;
  urz = pi->uz - pj->uz;

  COMPUTE_MOMENTUM_TRANSFER(urx,ury,urz,ax,ay,az,rng);

  if( type & 1 ) {
    w = hs->twomu_mi; 
    pi->ux -= w*ax;
    pi->uy -= w*ay;
    pi->uz -= w*az;
  }

  if( type & 2 ) {
    w = hs->twomu_mj; 
    pj->ux += w*ax;
    pj->uy += w*ay;
    pj->uz += w*az;
  }
}

#undef CMOV

void
checkpt_hard_sphere( const hard_sphere_t * hs ) {
  CHECKPT( hs, 1 );
}

hard_sphere_t *
restore_hard_sphere( void ) {
  hard_sphere_t * hs;
  RESTORE( hs );
  return hs;
}

/* Public interface **********************************************************/

collision_op_t *
hard_sphere_fluid( const char * RESTRICT name, /* Model name */
                   const float n0,             /* Fluid density (#/VOLUME) */
                   const float vdx,            /* Fluid x-drift (VELOCITY) */
                   const float vdy,            /* Fluid y-drift (VELOCITY) */
                   const float vdz,            /* Fluid z-drift (VELOCITY) */
                   const float kT0,            /* Fluid temperature (ENERGY) */
                   const float m0,             /* Fluid p. mass (MASS) */
                   const float r0,             /* Fluid p. radius (LENGTH) */
                   species_t * RESTRICT sp,    /* Species */
                   const float rsp,            /* Species p. radius (LENGTH) */
                   rng_pool_t * RESTRICT rp,   /* Entropy pool */
                   const int interval ) {      /* How often to apply this */
  hard_sphere_t * hs;

  if( n0<0 || kT0<0 || m0<=0 || r0<0 ||
      !sp || sp->m<=0 || rsp<0 ) ERROR(( "Bad args" ));

  MALLOC( hs, 1 );
  hs->twomu_mi = 2.*m0   /(sp->m + m0);
  hs->twomu_mj = 2.*sp->m/(sp->m + m0);
  hs->Kc       = M_PI*(r0+rsp)*(r0+rsp)*sp->g->cvac;
  hs->udx      = vdx / sp->g->cvac;
  hs->udy      = vdy / sp->g->cvac;
  hs->udz      = vdz / sp->g->cvac;
  hs->ut       = sqrtf(kT0/(m0*sp->g->cvac*sp->g->cvac));

  hs->ut2           = kT0/(m0*sp->g->cvac*sp->g->cvac);
  hs->alpha_Kt2ut4  = (8./M_PI)*(hs->Kc*n0*hs->Kc*n0)*(hs->ut2*hs->ut2);
  hs->beta_Kt2ut2   = (4./(12.-3.*M_PI))*(hs->Kc*n0*hs->Kc*n0)*hs->ut2;
  hs->gamma_Kt2     = ((3.*M_PI-8.)/(24.-6*M_PI))*(hs->Kc*n0*hs->Kc*n0);
  hs->ut2          += FLT_MIN;

  REGISTER_OBJECT( hs, checkpt_hard_sphere, restore_hard_sphere, NULL );
  return unary_collision_model( name,
                   (unary_rate_constant_func_t)hard_sphere_fluid_rate_constant,
                   (unary_collision_func_t)    hard_sphere_fluid_collision,
                                hs, sp, rp, interval );
}

collision_op_t *
hard_sphere( const char * RESTRICT name, /* Model name */
             species_t * RESTRICT spi,   /* Species-i */
             const float ri,             /* Species-i p. radius (LENGTH) */
             species_t * RESTRICT spj,   /* Species-j */
             const float rj,             /* Species-j p. radius (LENGTH) */
             rng_pool_t * RESTRICT rp,   /* Entropy pool */
             const double sample,        /* Sampling density */
             const int interval ) {      /* How often to apply this */
  hard_sphere_t * hs;

  if( !spi || spi->m<=0 || ri<0 ||
      !spj || spj->m<=0 || rj<0 || spi->g!=spj->g ) ERROR(( "Bad args" ));

  MALLOC( hs, 1 );
  CLEAR(  hs, 1 );
  hs->twomu_mi = 2*spj->m/(spi->m+spj->m);
  hs->twomu_mj = 2*spi->m/(spi->m+spj->m);
  hs->Kc       = spi->g->cvac*M_PI*(ri+rj)*(ri+rj);

  REGISTER_OBJECT( hs, checkpt_hard_sphere, restore_hard_sphere, NULL );
  return binary_collision_model( name,
                        (binary_rate_constant_func_t)hard_sphere_rate_constant,
                        (binary_collision_func_t)    hard_sphere_collision,
                                 hs, spi, spj, rp, sample, interval );
}

