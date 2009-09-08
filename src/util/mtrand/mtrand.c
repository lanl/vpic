/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Adapted from earlier V4PIC versions.
 *
 */

#include "mtrand.h"
#include "mtrand_conv.h" // For drand53_o, drand53_c0, drand53_c1, ...
#include "../checkpt/checkpt.h"

/*******************************************************
 * Mersenne-Twister 19937 Random Number Generator Core *
 *******************************************************/

#define MT_N          624
#define MT_M          397
#define MT_TWIST(u,v) (((((u)&0x80000000UL)|((v)&0x7fffffffUL))>>1)^((-((v)&1))&0x9908b0dfUL))
#define MT_TEMPER(u)  (u) ^= ( (u) >> 11 );                \
                      (u) ^= ( (u) << 7  ) & 0x9d2c5680UL; \
                      (u) ^= ( (u) << 15 ) & 0xefc60000UL; \
                      (u) ^= ( (u) >> 18 )

struct mt_rng {
  uint32_t next;
  uint32_t state[MT_N];
};

static void
mt_next_state( mt_rng_t * rng ) {
  int j;
  uint32_t * p;

  rng->next = 0;
  p = rng->state;
  for( j=MT_N-MT_M+1; --j; p++ ) p[0] = p[MT_M]      ^ MT_TWIST( p[0], p[1] );
  for( j=MT_M       ; --j; p++ ) p[0] = p[MT_M-MT_N] ^ MT_TWIST( p[0], p[1] );
  /**/                           p[0] = p[MT_M-MT_N] ^ MT_TWIST( p[0], rng->state[0] );
}

#define URAND32( rng, y )                               \
  if( (rng)->next==MT_N ) mt_next_state( rng );         \
  (y) = (rng)->state[ (rng)->next++ ];			\
  MT_TEMPER(y)

/*****************************************************************************
 * Constructors and destructors                                              *
 *****************************************************************************/

/* Though the checkpt/restore functions are not part of the public
   API, they must not be declared static */

void
checkpt_mt_rng( const mt_rng_t * rng ) {
  CHECKPT( rng, 1 );
}

mt_rng_t *
restore_mt_rng( void ) {
  mt_rng_t * rng;
  RESTORE( rng );
  return rng;
}

mt_rng_t *
new_mt_rng( unsigned int seed ) {
  mt_rng_t * rng;
  MALLOC( rng, 1 );
  seed_mt_rng( rng, seed );
  REGISTER_OBJECT( rng, checkpt_mt_rng, restore_mt_rng, NULL );
  return rng;
}

void
delete_mt_rng( mt_rng_t * rng ) {
  if( rng==NULL ) return;
  UNREGISTER_OBJECT( rng );
  FREE( rng );
}

/*****************************************************************************
 * Seeders and related functions                                             *
 *****************************************************************************/

void
seed_mt_rng( mt_rng_t * rng,
             unsigned int seed ) {
  int j;
  rng->next = MT_N;
  rng->state[0] = seed ^ 0x900df00c; // state[0] on srand(1) is goodfood
  for( j=1; j<MT_N; j++ )
    rng->state[j] = 1812433253*(rng->state[j-1]^(rng->state[j-1]>>30))+j;
}

// mt_getsize returns the number of chars required to hold the
// generator's internal state in the format used by mt_getstate
// function. This routine assumes 8-bit "char"s

size_t
get_mt_rng_size( mt_rng_t * rng ) {
  return 4*(MT_N+1);
}

// mt_getstate saves the state of the generator in a machine
// independent format on machines with 8-bits "char"s (true on
// virtually all hardware in the last 30 years) if the char array is
// large enough to hold it.

void
get_mt_rng_state( mt_rng_t * rng,
                  void * _s ) {
  uint8_t * s = (uint8_t *)_s;
  size_t j;

  // Serialize rng->next
  *(s++) = (rng->next & 0x000000ff) >> 0;
  *(s++) = (rng->next & 0x0000ff00) >> 8;
  *(s++) = (rng->next & 0x00ff0000) >> 16;
  *(s++) = (rng->next & 0xff000000) >> 24;

  // Serialize rng->state
  for( j=0; j<MT_N; j++ ) {
    *(s++) = (rng->state[j] & 0x000000ff) >> 0;
    *(s++) = (rng->state[j] & 0x0000ff00) >> 8;
    *(s++) = (rng->state[j] & 0x00ff0000) >> 16;
    *(s++) = (rng->state[j] & 0xff000000) >> 24;
  }
}

// mt_setstate sets the state of the generator. It can take a state
// provided by mt_getstate or a user designed state of at least 5
// bytes. A valid user designed state has at least one non-zero char
// among s[4] through s[min(n,get_mt_rng_size(h))-1]. This routine
// assumes 8-bit chars.

void
set_mt_rng_state( mt_rng_t * rng,
                  const void * _s,
                  size_t n ) {
  const uint8_t * s = (const uint8_t *)_s;
  size_t j, k;

  // Extract rng->next
  rng->next  = ((uint32_t)(s[0])) << 0;
  rng->next |= ((uint32_t)(s[1])) << 8;
  rng->next |= ((uint32_t)(s[2])) << 16;
  rng->next |= ((uint32_t)(s[3])) << 24;
  if( rng->next<0 || rng->next>MT_N ) rng->next=MT_N;

  // Extract rng->state
  for( j=0, k=4; j<MT_N; j++ ) {
    rng->state[j]  = ((uint32_t)(s[k++])) << 0;  if( k==n ) k=4;
    rng->state[j] |= ((uint32_t)(s[k++])) << 8;  if( k==n ) k=4;
    rng->state[j] |= ((uint32_t)(s[k++])) << 16; if( k==n ) k=4;
    rng->state[j] |= ((uint32_t)(s[k++])) << 24; if( k==n ) k=4;
  }
}

/*****************************************************************************
 * Generate integer random numbers                                           *
 *****************************************************************************/

// The irf code-generator is very portable but not maximally efficient
// when dealing with short and char datatypes (for example, on a
// 32-bit system with 8-bit chars, 4 chars could be made for every
// 32-bit rand). Note: any self-respecting compiler will optimize the
// bit shift at compile time

#define INT_RNG( name, type, s )                                        \
  type                                                                  \
  mt_##name( mt_rng_t * rng ) {                                         \
    uint32_t y;                                                         \
    URAND32( rng, y );                                                  \
    return (type)(y>>(s+((CHAR_BIT*sizeof(type)<32)?                    \
                         (32-CHAR_BIT*sizeof(type)):0)));               \
  }                                                                     \
                                                                        \
  void                                                                  \
  mt_##name##_fill( mt_rng_t * rng,                                     \
                    type * x,                                           \
                    size_t n ) {                                        \
    uint32_t y;                                                         \
    for(;n;n--) {                                                       \
      URAND32( rng, y );                                                \
      *(x++)=(type)(y>>(s+((CHAR_BIT*sizeof(type)<32)?                  \
                           (32-CHAR_BIT*sizeof(type)):0)));             \
    }                                                                   \
  }

INT_RNG( crand,  signed char,        1 ) // Force signed chars
INT_RNG( hrand,  signed short int,   1 )
INT_RNG( rand,   signed int,         1 )
INT_RNG( lrand,  signed long int,    1 )

INT_RNG( ucrand, unsigned char,      0 )
INT_RNG( uhrand, unsigned short int, 0 )
INT_RNG( urand,  unsigned int,       0 )
INT_RNG( ulrand, unsigned long int,  0 )

#undef INT_RNG

/*****************************************************************************
 * Generate floating point random numbers                                    *
 *****************************************************************************/

#define FLT_RNG( name, type, which )                                    \
  type                                                                  \
  mt_##name( mt_rng_t * rng ) {                                         \
    uint32_t a;                                                         \
    URAND32( rng, a );                                                  \
    return which(a);                                                    \
  }                                                                     \
                                                                        \
  void                                                                  \
  mt_##name##_fill( mt_rng_t * rng,                                     \
                    type * x,                                           \
                    size_t n ) {                                        \
    uint32_t a;                                                         \
    for(;n;n--) {							\
      URAND32( rng, a );                                                \
      *(x++) = which(a);						\
    }                                                                   \
  }

#define DBL_RNG( name, which )                                     \
  double                                                           \
  mt_##name( mt_rng_t * rng ) {                                    \
    uint32_t a, b;                                                 \
    URAND32( rng, a );                                             \
    URAND32( rng, b );                                             \
    return which(a,b);                                             \
  }                                                                \
                                                                   \
  void                                                             \
  mt_##name##_fill( mt_rng_t * rng,                                \
                    double * x,                                    \
                    size_t n ) {                                   \
    uint32_t a, b;                                                 \
    for(;n;n--) {                                                  \
      URAND32( rng, a );                                           \
      URAND32( rng, b );                                           \
      *(x++) = which(a,b);                                         \
    }                                                              \
  }

FLT_RNG( frand,         float,  frand24_o  )
FLT_RNG( frand_c0,      float,  frand24_c0 )
FLT_RNG( frand_c1,      float,  frand24_c1 )
FLT_RNG( frand_c,       float,  frand24_c  )

FLT_RNG( fast_drand,    double, drand32_o  )
FLT_RNG( fast_drand_c0, double, drand32_c0 )
FLT_RNG( fast_drand_c1, double, drand32_c1 )
FLT_RNG( fast_drand_c,  double, drand32_c  )

DBL_RNG( drand,    drand53_o  )
DBL_RNG( drand_c0, drand53_c0 )
DBL_RNG( drand_c1, drand53_c1 )
DBL_RNG( drand_c,  drand53_c  )

#undef FLT_RNG
#undef DBL_RNG

/*****************************************************************************
 * Generate a normal random number. Range is (-inf,inf)                      *
 * f(x) = exp( -x^2 / 2 ) / sqrt( 2*pi )                                     *
 *****************************************************************************/

/* The ziggurat method is a very efficient method for random number
   generation of variables with monotonically decreasing PDFs.

   For a Gaussian, if (x,y) are picked randomly but uniformly in the
   region 0 <= y <= f(x) = exp(-x^2/2) for x in (-inf,inf) then the
   PDF of x alone is Gaussian.  Rejection testing is based on the
   observation that if we could pick (x,y) in a region uniformly that
   completely covers this and reject the points that are outside, the
   PDF of the remaining x would also be Gaussian distribted.  The goal
   then is to find a region for which random points can be generated
   uniformly efficiently and that tightly covers the original region.
   The ziggurat method constructs such a region then also uses
   monotonicity to improve performance further.

   To be specific, cover f(x) with N regions of equal area enumerated
   0:N-1.  Let regions 0:N-2 be the strips:
     ( -x_{i+1},+x_{i+1} ) x ( f(x_{i+1}), f(x_i) )
   and let x_0 = 0 and x_{i+1} > x_i.  Let region N-1 be the
   rectangular strip:
     ( -x_{N-1},+x_{N-1} ) x ( 0, f(x_{N-1}) )
   plus a tail for |x|>=x_{N-1}, 0<=y<=t(x) where exp(-x^2/2)<=t(x)
   for |x|>x_{N-1}.

   This set of regions (or ziggurat) tightly and completely covers the
   original region.  Since each individual region has the same area,
   picking (x,y) uniformly within that region consists of:

     (1) Picking a region at random uniformly.
     (2) Picking a point uniformly within that region

   Step 1 is trivial and step 2 is trivial for all but region N-1.
   This yields the beginning of the Ziggurat algortihm:

   for(;;)
     i = rand uniformly from [0,1,...N-1]
     if i!=N-1,
       x = rand uniformly from ( -x_{i+1}, +x_{i+1} )
       y = rand uniformly from ( f(x_{i+1}, f(x_i)  ]
       if y < exp(-x^2/2), return x
     else
       handle region N-1

   This can be dramatically optimized if it is noted that, for all but
   the last region, if the x-coordinate less than x_i, the
   y-coordinate is guaranteed to fall under the Gaussian by
   monotonicity!  This yields the more efficient:

   for(;;)
     i = rand uniformly from [0,1,...N-1]
     if i!=N-1,
       x = rand uniformly from ( -x_{i+1}, +x_{i+1} )
       if |x| < x_i, return x
       y = rand uniformly from ( f(x_{i+1}, f(x_i)  ]
       if y < exp(-x^2/2), return x
     else
       handle region N-1

   Let 2v be the total area of each region and let r = x_{N-1}.  The
   last region handled more elegantly by noting r f(r) / v percent of
   the time, a uniform random point in region N-1 comes from the
   rectangular part (2 r f(r) is the rectangle area and 2v is the
   total area).  Letting x_N = v/f(r) then, the quick acceptance test
   can account for the rectangular part of region N-1.  However, in
   the last region when not quickly accepted, we need to compute a
   point in the uniformly in the tail region.  This yields:

   for(;;)
     i = rand uniformly from [0,1,...N-1]
     x = rand uniformly from ( -x_{i+1}, +x_{i+1} )
     if |x| < x_i, return x
     if i!=N-1,
       y = rand uniformly from ( f(x_{i+1}, f(x_i)  ]
     else
       (x,y) = rand uniformly from region N-1 tails
     if y < exp(-x^2/2), return x

   Using an exponential decay that matches with the first derivative
   at f(x) at r is a good choice for the tail:
     g(x) = exp(-r^2/2 ) exp( -r ( |x| - r ) )
          = exp( -r ( |x| - r/2 ) )
   An x with the distribution of g is easily generated directly, and,
   given such an x, a y uniformly distributed on ( 0, g(x) ] yields a
   (x,y) uniformly from the tail region.  This yields:

   for(;;)
     i = rand uniformly from [0,1,...N-1]
     x = rand uniformly from ( -x_{i+1}, +x_{i+1} )
     if |x| < x_i, return x
     if i!=N-1,
       y = rand uniformly from ( f(x_{i+1}, f(x_i)  ]
     else
       x = rand uniformly from [+1,-1] *
           ( r - (1/r) log rand uniformly from (0,1] )
       y = rand uniformly from ( 0, exp( -r ( |x| - r/2 ) ) ]
     if y < exp(-x^2/2), return x

   The last bit of magic is that if dealing with single precision
   randomness, i and x may be extracted from the bits of a single
   unsigned 32-bit random integer.  Specifically, let N=64 and:

     u = 32-bit rand
     s = bit  0    of u (      a  1-bit rand)
     i = bits 1:6  of u (2   * a  6-bit rand)
     k = bit  7    of u (2^7 * a  1-bit rand)
     j = bits 8:31 of u (2^8 * a 24-bit rand)

   Then, an appropriate x (uniform rand with no finite precision
   biases in the interval [-x_{i+1},+x_{i+1}]:
     z = (s?-1:1) (j+2k) / 2^32
     x = x_{i+1} z

   Note that j+2k is a 2^8 times a rand on [0,2^24] p_0 = p_{2^24 } =
   1/2^25 and p_n = 1/2^24 otherwise ("trapezoidal rand").  As a
   result, it is as though z was picked in continuum interval (-1,+1)
   and rounded to the nearest representable point in the desired
   interval from the 2^25+1 point lattice uniformly (even signed zeros
   are picked with sanely).  Due to the uniform choice of sign and the
   number of bits in j, the computation of z is exact and yields a
   full precision float.

   Note that given all the above:
     v     = r f(r) + exp(-r^2/2) / r
     x_N   = v/f(r)
     x_N-1 = r
     x_i   = inverse_f( f(x_{i+1}) + v/x_{i+1} )
     x_0   = 0
   and r can be found iteratively.  Letting:
     N        = 64
     R        = x_{N-1}
     scale    = 1/2^32
     zig_x[i] = x_i
     zig_y[i] = f(x_i)
   yields the freaky efficient algorithm below.  The vast majority of
   the normals are generated with a single u32, some quick table
   lookups and floating point multiplies.

   Similarly considerations hold for the double precision variant.
   There a 64-bit rand is used instead of a 32-bit rand to generate a
   53-bit trapezoid rand, 8-bit index, 1 bit sign and N = 256 and
   scale = 1/2^64. */

double
mt_drandn( mt_rng_t * rng ) {
  uint32_t a, b, i, s;
  double x, y, j;

  static const double sgn[2] = { 1., -1. };

# include "drandn_table.h"

  for(;;) {

    // Extract components of a 64-bit rand 

    URAND32( rng, a );
    URAND32( rng, b );

    s = ( a & 0x00000001 );      //         1-bit uniform   rand 
    i = ( a & 0x000001fe ) >> 1; //         8-bit uniform   rand
    j = 4294967296.*b + (( a & 0xfffff800 ) + (( a & 0x00000400 ) << 1));
                                 // 2^11 ( 53-bit trapezoid rand )

    // Construct |x| and see if we can accept this point 
    // FIXME: COULD PRECOMPUTE SCALE * ZIG_X[i+1] AND/OR
    // EVEN DO THE COMPARE DIRECTLY ON j!

    x = j*(scale*zig_x[i+1]);
    if( x<zig_x[i] ) break; /* Vast majority of the time */
 
    // Construct a y for rejection testing 
 
    URAND32( rng, a );
    URAND32( rng, b );
    y = drand53_c( a, b );
    if( i!=N-1 ) y = zig_y[i]+(zig_y[i+1]-zig_y[i])*y;
    else { /* In tail */
      URAND32( rng, a );
      URAND32( rng, b );
      x  = R - (1./R)*log( drand53_c1( a, b ) );
      y *= exp( -R*( x - 0.5*R ) );
    }

    if( y < exp(-0.5*x*x) ) break;
 }

 return sgn[s]*x; // FIXME: Use copysign, trinary or branch? 
}

void
mt_drandn_fill( mt_rng_t * rng,
                double * x,
                size_t n ) {
  for( ; n; n-- ) *(x++) = mt_drandn( rng );
}

float
mt_frandn( mt_rng_t * rng ) {
  uint32_t a, i, j, s;
  float x, y;

  static const float sgn[2] = { 1.f, -1.f };

# include "frandn_table.h"

  for(;;) {

    // Extract components of a 32-bit rand 

    URAND32( rng, a );
    s = ( a & 0x00000001 );      //        1-bit uniform   rand 
    i = ( a & 0x0000007e ) >> 1; //        6-bit uniform   rand
    j = ( a & 0x00000080 ) << 1; // 2^8 (  1-bit uniform   rand )
    j = ( a & 0xffffff00 ) + j;  // 2^8 ( 24-bit trapezoid rand )

    // Construct |x| and see if we can accept this point 
    // FIXME: COULD PRECOMPUTE SCALE * ZIG_X[i+1] AND/OR
    // EVEN DO THE COMPARE DIRECTLY ON j!

    x = j*(scale*zig_x[i+1]);
    if( x<zig_x[i] ) break; /* Vast majority of the time */
 
    // Construct a y for rejection testing 
 
    URAND32( rng, a );
    y = frand24_c(a);
    if( i!=N-1 ) y = zig_y[i]+(zig_y[i+1]-zig_y[i])*y;
    else { /* In tail */
      URAND32( rng, a );
      x  = R - (1.f/R)*logf( frand24_c1(a) );
      y *= expf( -R*( x - 0.5f*R ) );
    }

    if( y < expf(-0.5f*x*x) ) break;
 }

 return sgn[s]*x; // FIXME: Use copysign, trinary or branch? 
}

void
mt_frandn_fill( mt_rng_t * rng,
                float * x,
                size_t n ) {
  for( ; n; n-- ) *(x++) = mt_frandn( rng );
}

/*****************************************************************************
 * Generate an exponential random number. Range is (0,inf)                   *
 * f(x) = exp(-x) for x>0, 0 otherwise                                       *
 *****************************************************************************/

double
mt_drande( mt_rng_t * rng ) {
  uint32_t a, b;
  URAND32( rng, a );
  URAND32( rng, b );
  return -log( drand53_c1(a,b) );
}

void
mt_drande_fill( mt_rng_t * rng,
                double * x,
                size_t n ) {
  uint32_t a, b;

  for( ; n; n-- ) {
    URAND32( rng, a );
    URAND32( rng, b );
    *(x++) = -log( drand53_c1(a,b) );
  }
}

/*****************************************************************************
 * Fill x with a random permutation of the integers {0,1,2,...,n-2,n-1}      *
 * Routine generally valid for 1 <= n <= min(2^32,int_max)                   *
 *****************************************************************************/

void
mt_randperm( mt_rng_t * rng,
             int * x,
             int n ) {
  uint32_t a, d;
  int i, t, r;

  // Create the initial permutation

  for( i=0; i<n; i++ ) x[i] = i;

  // Apply a random swap to each element of the permutation such that
  // any of the n! permutations could be generated with equal
  // probability. Note: The method used to pick a random number on
  // [0...t-1] ([0...n-i-1]) uses the highest quality (most
  // significant) bits of a URAND32 to determine an _exactly_
  // uniformly distributed rand. A simpler "mod n" method uses the
  // least significant bits and is not exactly uniform unless 2^32 is
  // an exact multiple of t. Note: the d calculation is done with
  // double precision because 2^32 cannot be represented with 32-bit
  // ints and not all 32-bit integers have an exact single precision
  // representation.

  for( i=0; i<n-1; i++ ) {
    t = n-i;
    d = (uint32_t)(4294967296./(double)t);
    do {
      URAND32( rng, a );
      r = a/d;
    } while( r>=t );
    r += i;
    t    = x[i];
    x[i] = x[r];
    x[r] = t;
  }
}

/*****************************************************************************
 * Shuffle the array x. x has n members. Each member has size s              *
 * This routine works for 1 <= n <= min(2^32,size_t_max)                     *
 *****************************************************************************/

void
mt_shuffle( mt_rng_t * rng,
            void * x,
            size_t n,
            size_t s ) {
  uint32_t a, d;
  size_t i, t, r;
  char *xi, *xr, c;

  // See randperm comment
  for( i=0; i<n-1; i++ ) {
    t = n-i;
    d = (size_t)((double)4294967296./(double)t);
    do {
      URAND32( rng, a );
      r = a/d;
    } while( r>=t );
    r += i;
    xi = ((char *)x) + s*i;
    xr = ((char *)x) + s*r;
    for( t=s; t; t-- ) {
      c       = *xi;
      *(xi++) = *xr;
      *(xr++) = c;
    }
  }
}
