#define IN_rng
#include "rng_private.h"
#include "../checkpt/checkpt.h"

/* Private API **************************************************************/

void
checkpt_rng( const rng_t * RESTRICT r ) {
  CHECKPT_ALIGNED( r, 1, 128 );
}

rng_t *
restore_rng( void ) {
  rng_t * r;
  RESTORE_ALIGNED( r );
  return r;
}

/* Public API ***************************************************************/

/* Structors */

rng_t *
new_rng( int seed ) {
  rng_t * r;
  MALLOC_ALIGNED( r, 1, 128 );
  seed_rng( r, seed );
  REGISTER_OBJECT( r, checkpt_rng, restore_rng, NULL );
  return r;
}

void
delete_rng( rng_t * r ) {
  if( !r ) return;
  UNREGISTER_OBJECT( r );
  FREE_ALIGNED( r );
}

/* Seeders */

#define u(n) r->state.u32[(n)]

static void
adjust_rng( rng_t * RESTRICT r ) {
  static const uint32_t parity[4] =
    { SFMT_PARITY0, SFMT_PARITY1, SFMT_PARITY2, SFMT_PARITY3 };
  uint32_t bit = 0;
  int n;

  /* Check if the generator is on a cycle with the desired period. */

  for( n=0; n<4; n++  ) bit ^= ( u(n) & parity[n] );
  for( n=16; n; n>>=1 ) bit ^= ( bit >> n );
  bit &= 1;
  if( bit==1 ) return;

  /* Nope.  Adjust the generator */

  for( n=0; n<4; n++ )
    for( bit=1; bit; bit<<=1 )
      if( (bit & parity[n]) ) { u(n) ^= bit; return; }
}

rng_t *
seed_rng( rng_t * RESTRICT r,
          int seed ) {
  int n;
  if( !r ) ERROR(( "Bad args" ));
  u(0) = (uint32_t)seed;
  for( n=1; n<SFMT_N32; n++ )
    u(n) = ((uint32_t)1812433253) * (u(n-1)^(u(n-1)>>30)) + n;
  adjust_rng( r );
  r->n = SFMT_NC;
  return r;
}

#undef u

/* Uniform integer generators */

#define _( type, prefix, state_prefix, is_signed )              \
type                                                            \
prefix##rand( rng_t * RESTRICT r ) {                            \
  type x;                                                       \
  RNG_NEXT( x, type, r, state_prefix, is_signed );              \
  return x;                                                     \
}                                                               \
                                                                \
type *                                                          \
prefix##rand_fill( rng_t * RESTRICT r,                          \
                   type  * RESTRICT x,                          \
                   size_t str_ele,                              \
                   size_t n_ele ) {                             \
  size_t n;                                                     \
  if( !n_ele ) return x;                                        \
  if( !r || !x ) ERROR(( "Bad args" ));                         \
  for( n=0; n<n_ele; n++ )                                      \
    RNG_NEXT( x[n*str_ele], type, r, state_prefix, is_signed ); \
  return x;                                                     \
}

_( char,    c,   uc,  1 ) _( unsigned char,  uc,  uc,  0 )
_( short,   h,   uh,  1 ) _( unsigned short, uh,  uh,  0 )
_( int,     i,   ui,  1 ) _( unsigned int,   ui,  ui,  0 )
_( long,    l,   ul,  1 ) _( unsigned long,  ul,  ul,  0 )
_( int8_t,  i8,  u8,  1 ) _( uint8_t,        u8,  u8,  0 )
_( int16_t, i16, u16, 1 ) _( uint16_t,       u16, u16, 0 )
_( int32_t, i32, u32, 1 ) _( uint32_t,       u32, u32, 0 )
_( int64_t, i64, u64, 1 ) _( uint64_t,       u64, u64, 0 )

#undef _

/* Uniform floating point generators */

#define _( type, prefix, variant, state_type, state_prefix )    \
type                                                            \
prefix##rand##variant( rng_t * RESTRICT r ) {                   \
  state_type u;                                                 \
  RNG_NEXT( u, state_type, r, state_prefix, 0 );                \
  return conv_##prefix##rand##variant( u );                     \
}                                                               \
                                                                \
type *                                                          \
prefix##rand##variant##_fill( rng_t * RESTRICT r,               \
                              type  * RESTRICT x,               \
                              size_t str_ele,                   \
                              size_t n_ele ) {                  \
  state_type u;                                                 \
  size_t n;                                                     \
  if( !n_ele ) return x;                                        \
  if( !r || !x ) ERROR(( "Bad args" ));                         \
  for( n=0; n<n_ele; n++ ) {                                    \
    RNG_NEXT( u, state_type, r, state_prefix, 0 );              \
    x[n*str_ele] = conv_##prefix##rand##variant( u );           \
  }                                                             \
  return x;                                                     \
}

_( float, f, ,    uint32_t, u32 ) _( double, d, ,    uint64_t, u64 )
_( float, f, _c0, uint32_t, u32 ) _( double, d, _c0, uint64_t, u64 )
_( float, f, _c1, uint32_t, u32 ) _( double, d, _c1, uint64_t, u64 )
_( float, f, _c,  uint32_t, u32 ) _( double, d, _c,  uint64_t, u64 )

#undef _

/* Normal generators */

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

#include "frandn_table.h"

float
frandn( rng_t * RESTRICT r ) {
  uint32_t a, i, j, s;
  float x, y;

  static const float scale = 1.f/4.294967296e+09f;
  static const float sgn[2] = { 1.f, -1.f };

  for(;;) {

    // Extract components of a 32-bit rand

#   if FRANDN_N!=64
#   error "frandn_table.h does not match frandn()"
#   endif

    RNG_NEXT( a, uint32_t, r, u32, 0 );
    s = ( a &   (uint32_t)0x01  );      //        1-bit uniform   rand
    i = ( a &   (uint32_t)0x7e  ) >> 1; //        6-bit uniform   rand
    j = ( a &   (uint32_t)0x80  ) << 1; // 2^8 (  1-bit uniform   rand )
    j = ( a & (~(uint32_t)0xff) ) + j;  // 2^8 ( 24-bit trapezoid rand )

    // Construct |x| and see if we can accept this point

    x = j*(scale*frandn_zig_x[i+1]);
    if( LIKELY( x<frandn_zig_x[i] ) ) break; // Vast majority of the time

    // Construct a y for rejection testing

    RNG_NEXT( a, uint32_t, r, u32, 0 );
    y = conv_frand_c(a);
    if( LIKELY( i!=FRANDN_N-1 ) )
      y = frandn_zig_y[i]+(frandn_zig_y[i+1]-frandn_zig_y[i])*y;
    else { // In tail
      RNG_NEXT( a, uint32_t, r, u32, 0 );
      x  = FRANDN_R - (1.f/FRANDN_R)*logf( conv_frand_c1(a) );
      y *= expf( -FRANDN_R*( x - 0.5f*FRANDN_R ) );
    }

    if( y < expf(-0.5f*x*x) ) break;
 }

 return sgn[s]*x; // FIXME: Use copysign, trinary or branch?
}

float *
frandn_fill( rng_t * RESTRICT r,
             float * RESTRICT x,
             size_t str_ele,
             size_t n_ele ) {
  size_t n;
  if( !n_ele ) return x;
  if( !r || !x ) ERROR(( "Bad args" ));
  for( n=0; n<n_ele; n++ ) x[ n*str_ele ] = frandn( r );
  return x;
}

#include "drandn_table.h"

double
drandn( rng_t * RESTRICT r ) {
  uint64_t a, i, j, s;
  double x, y;

  static const double scale = 1./1.8446744073709551616e+19;
  static const double sgn[2] = { 1., -1. };

  for(;;) {

    // Extract components of a 64-bit rand

#   if DRANDN_N!=256
#   error "drandn_table.h does not match drandn"
#   endif

    RNG_NEXT( a, uint64_t, r, u64, 0 );
    s =   a &   (uint64_t)0x001;        //         1-bit uniform   rand
    i = ( a &   (uint64_t)0x1fe ) >> 1; //         8-bit uniform   rand
    j = ( a &   (uint64_t)0x400 ) << 1; // 2^11 (  1-bit uniform   rand )
    j = ( a & (~(uint64_t)0x3ff)) + j;  // 2^11 ( 53-bit trapezoid rand )

    // Construct |x| and see if we can accept this point

    x = j*(scale*drandn_zig_x[i+1]);
    if( LIKELY( x<drandn_zig_x[i] ) ) break; // Vast majority of the time

    // Construct a y for rejection testing

    RNG_NEXT( a, uint64_t, r, u64, 0 );
    y = conv_drand_c(a);
    if( LIKELY( i!=DRANDN_N-1 ) )
      y = drandn_zig_y[i]+(drandn_zig_y[i+1]-drandn_zig_y[i])*y;
    else { // In tail
      RNG_NEXT( a, uint64_t, r, u64, 0 );
      x  = DRANDN_R - (1./DRANDN_R)*log( conv_drand_c1(a) );
      y *= exp( -DRANDN_R*( x - 0.5*DRANDN_R ) );
    }

    if( y < exp(-0.5*x*x) ) break;
 }

 return sgn[s]*x; // FIXME: Use copysign, trinary or branch?
}

double *
drandn_fill( rng_t  * RESTRICT r,
             double * RESTRICT x,
             size_t str_ele,
             size_t n_ele ) {
  size_t n;
  if( !n_ele ) return x;
  if( !r || !x ) ERROR(( "Bad args" ));
  for( n=0; n<n_ele; n++ ) x[ n*str_ele ] = drandn( r );
  return x;
}

/* Exponential random numbers */

/* Uses the transformation method */

float
frande( rng_t * RESTRICT r ) {
  uint32_t a;
  RNG_NEXT( a, uint32_t, r, u32, 0 );
  return -logf( conv_frand_c1(a) );
}

float *
frande_fill( rng_t * RESTRICT r,
             float * RESTRICT x,
             size_t str_ele,
             size_t n_ele ) {
  size_t n;
  if( !n_ele ) return x;
  if( !r || !x ) ERROR(( "Bad args" ));
  for( n=0; n<n_ele; n++ ) x[ n*str_ele ] = frande( r );
  return x;
}

double
drande( rng_t * RESTRICT r ) {
  uint64_t a;
  RNG_NEXT( a, uint64_t, r, u64, 0 );
  return -log( conv_drand_c1(a) );
}

double *
drande_fill( rng_t  * RESTRICT r,
             double * RESTRICT x,
             size_t str_ele,
             size_t n_ele ) {
  size_t n;
  if( !n_ele ) return x;
  if( !r || !x ) ERROR(( "Bad args" ));
  for( n=0; n<n_ele; n++ ) x[ n*str_ele ] = drande( r );
  return x;
}

/* Miscellaneous */

int *
randperm( rng_t * RESTRICT r,
          int   * RESTRICT x,
          int n ) {
  uint32_t a, d;
  int t, i, j;

  if( !n ) return x;
  if( !r || !x ) ERROR(( "Bad args" ));

  // Create the initial permutation

  for( i=0; i<n; i++ ) x[i] = i;

  // Apply a random swap to each element of the permutation such that
  // any of the n! permutations could be generated with equal
  // probability. Note: The method used to pick a random number on
  // [0...t-1] ([0...n-i-1]) uses the highest quality (most
  // significant) bits of a URAND32 to determine an _exactly_
  // uniformly distributed rand.  A simpler "mod n" method uses the
  // least significant bits and is not exactly uniform unless 2^32-1
  // is an exact multiple of t.

  for( i=0; i<n-1; i++ ) {
    d = UINT32_MAX / (uint32_t)(n-i);
    do {
      RNG_NEXT( a, uint32_t, r, u32, 0 );
      j = i + (int)(a/d);
    } while( UNLIKELY( j>=n ) );
    t = x[i], x[i] = x[j], x[j] = t;
  }

  return x;
}

void *
shuffle( rng_t * RESTRICT r,
         void  * RESTRICT _x,
         size_t sz_ele,
         size_t str_ele,
         size_t n_ele ) {
  uint32_t a, d;
  size_t i, j, k;

  if( n_ele<2 || !sz_ele ) return _x;
  if( !r || !_x || sz_ele>str_ele ) ERROR(( "Bad args" ));

  // See randperm comment

# define SHUFFLE(type) do {                     \
    type * RESTRICT x = (type  *)_x, t;         \
    for( i=0; i<n_ele-1; i++ ) {                \
      d = UINT64_MAX / (uint64_t)(n_ele-i);     \
      do {                                      \
        RNG_NEXT( a, uint64_t, r, u64, 0 );     \
        j = i + (size_t)(a/d);                  \
      } while( UNLIKELY( j>=n_ele ) );          \
      SWAP_IJ;                                  \
    }                                           \
    return _x;                                  \
  } while(0)


  if( (str_ele % sz_ele)==0 ) {

#   define SWAP_IJ                                                      \
    t = x[i*str_ele], x[i*str_ele] = x[j*str_ele], x[j*str_ele] = t

    switch( sz_ele ) {
    case 1:  str_ele /= sz_ele; SHUFFLE( uint8_t  );
    case 2:  str_ele /= sz_ele; SHUFFLE( uint16_t );
    case 4:  str_ele /= sz_ele; SHUFFLE( uint32_t );
    case 8:  str_ele /= sz_ele; SHUFFLE( uint64_t );
    default: break;
    }

#   undef SWAP_IJ

  }

# define SWAP_IJ                                                      \
  for( k=0; k<sz_ele; k++ )                                           \
    t = x[i*str_ele+k], x[i*str_ele+k] = x[j*str_ele+k], x[j*str_ele+k] = t

  SHUFFLE( uint8_t );

# undef SWAP_IJ

  return _x;
}

