#ifndef _rng_spe_h_
#define _rng_spe_h_
#ifdef __SPU__

/* The SPU implementation gives you much of the functionality of the
   regular implementation but most of these are declared as inline. */

#define IN_rng
#include "rng_private.h"

STATIC_INLINE rng_t * ALIGNED(128)
mfc_get_rng( MEM_PTR( rng_t, 128 ) _r ) {
  rng_t * r;
  SPU_MALLOC( r, 1, 128 );
  mfc_get( r, _r, sizeof(*r), 31, 0, 0 );
  mfc_write_tag_mask( 1<<31 );
  mfc_read_tag_status_all();
  return r;
}

STATIC_INLINE void
mfc_put_rng( rng_t * ALIGNED(128) r,
             MEM_PTR( rng_t, 128 ) _r ) {
  mfc_put( r, _r, sizeof(*r), 31, 0, 0 );
  mfc_write_tag_mask( 1<<31 );
  mfc_read_tag_status_all();
}

/* Uniform integer generators */

#define _( type, prefix, state_prefix, is_signed ) \
STATIC_INLINE type                                 \
prefix##rand( rng_t * RESTRICT r ) {               \
  type x;                                          \
  RNG_NEXT( x, type, r, state_prefix, is_signed ); \
  return x;                                        \
}

#ifdef NEED_crand
_( char,    c,   uc,  1 ) 
#endif
#ifdef NEED_hrand
_( short,   h,   uh,  1 ) 
#endif
#ifdef NEED_irand
_( int,     i,   ui,  1 ) 
#endif
#ifdef NEED_lrand
_( long,    l,   ul,  1 ) 
#endif
#ifdef NEED_i8rand
_( int8_t,  i8,  u8,  1 ) 
#endif
#ifdef NEED_i16rand
_( int16_t, i16, u16, 1 ) 
#endif
#ifdef NEED_i32rand
_( int32_t, i32, u32, 1 ) 
#endif
#ifdef NEED_i64rand
_( int64_t, i64, u64, 1 ) 
#endif
#ifdef NEED_ucrand
_( unsigned char,  uc,  uc,  0 )
#endif
#ifdef NEED_uhrand
_( unsigned short, uh,  uh,  0 )
#endif
#ifdef NEED_uirand
_( unsigned int,   ui,  ui,  0 )
#endif
#ifdef NEED_ulrand
_( unsigned long,  ul,  ul,  0 )
#endif
#ifdef NEED_u8rand
_( uint8_t,        u8,  u8,  0 )
#endif
#ifdef NEED_u16rand
_( uint16_t,       u16, u16, 0 )
#endif
#ifdef NEED_u32rand
_( uint32_t,       u32, u32, 0 )
#endif
#ifdef NEED_64rand
_( uint64_t,       u64, u64, 0 )
#endif

#undef _

/* Uniform floating point generators */

#define _( type, prefix, variant, state_type, state_prefix )    \
STATIC_INLINE type                                              \
prefix##rand##variant( rng_t * RESTRICT r ) {                   \
  state_type u;                                                 \
  RNG_NEXT( u, state_type, r, state_prefix, 0 );                \
  return conv_##prefix##rand##variant( u );                     \
}

#ifdef NEED_frand
_( float, f, ,    uint32_t, u32 )
#endif
#ifdef NEED_frand_c0
_( float, f, _c0, uint32_t, u32 )
#endif
#ifdef NEED_frand_c1
_( float, f, _c1, uint32_t, u32 )
#endif
#ifdef NEED_frand_c
_( float, f, _c,  uint32_t, u32 )
#endif
#ifdef NEED_drand
_( double, d, ,    uint64_t, u64 )
#endif
#ifdef NEED_drand_c0
_( double, d, _c0, uint64_t, u64 )
#endif
#ifdef NEED_drand_c1
_( double, d, _c1, uint64_t, u64 )
#endif
#ifdef NEED_drand_c
_( double, d, _c,  uint64_t, u64 )
#endif

#undef _

/* Normal floating point generators */

#ifdef NEED_frandn
#include "frandn_table.h"
STATIC_INLINE float
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
    // FIXME: COULD PRECOMPUTE SCALE * ZIG_X[i+1] AND/OR
    // EVEN DO THE COMPARE DIRECTLY ON j!

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
#endif

#ifdef NEED_drandn
#include "drandn_table.h"
STATIC_INLINE double
drandn( rng_t * RESTRICT r ) {
  uint64_t a, i, j, s;
  double x, y;

  static const double scale = 1./1.8446744073709551616e+19;
  static const double sgn[2] = { 1., -1. };

  for(;;) {
    
    // Extract components of a 64-bit rand 

#   if DRANDN_N!=256
#   error "drandn_table.h does not math drandn()"
#   endif
    
    RNG_NEXT( a, uint64_t, r, u64, 0 );
    s =   a &   (uint64_t)0x001;        //         1-bit uniform   rand 
    i = ( a &   (uint64_t)0x1fe ) >> 1; //         8-bit uniform   rand
    j = ( a &   (uint64_t)0x400 ) << 1; // 2^11 (  1-bit uniform   rand )
    j = ( a & (~(uint64_t)0x3ff)) + j;  // 2^11 ( 53-bit trapezoid rand )

    // Construct |x| and see if we can accept this point 
    // FIXME: COULD PRECOMPUTE SCALE * ZIG_X[i+1] AND/OR
    // EVEN DO THE COMPARE DIRECTLY ON j!

    x = j*(scale*drandn_zig_x[i+1]);
    if( LIKELY( x<drandn_zig_x[i] ) ) break; // Vast majority of the time
 
    // Construct a y for rejection testing 
 
    RNG_NEXT( a, uint64_t, r, u64, 0 );
    y = conv_drand_c(a);
    if( LIKELY( i!=DRANDN_N-1 ) )
      y = drandn_zig_y[i]+(drandn_zig_y[i+1]-drandn_zig_y[i])*y;
    else { // In tail 
      RNG_NEXT( a, uint64_t, r, u64, 0 );
      x  = FRANDN_R - (1./FRANDN_R)*log( conv_drand_c1(a) );
      y *= exp( -FRANDN_R*( x - 0.5*FRANDN_R ) );
    }

    if( y < exp(-0.5*x*x) ) break;
 }

 return sgn[s]*x; // FIXME: Use copysign, trinary or branch? 
}
#endif

/* Exponential random numbers */

#ifdef NEED_frande
STATIC_INLINE float
frande( rng_t * RESTRICT r ) {
  uint32_t a;
  RNG_NEXT( a, uint32_t, r, u32, 0 );
  return -logf( conv_frand_c1(a) );
}
#endif

#ifdef NEED_drande
STATIC_INLINE double
drande( rng_t * RESTRICT r ) {
  uint64_t a;
  RNG_NEXT( a, uint64_t, r, u64, 0 );
  return -log( conv_drand_c1(a) );
}
#endif

/* Miscellaneous */

#ifdef NEED_randperm
STATIC_INLINE int *
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
#endif

#ifdef NEED_shuffle
STATIC_INLINE void *
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
#endif

#endif /* __SPU__ */
#endif /* _rng_spe_h_ */
