#ifndef _rng_private_h_
#define _rng_private_h_

#ifndef IN_rng
#error "Do not include rng_private.h; use rng.h"
#endif

#include "rng.h"

#if defined(__SSE2__) /* Use SSE-2 accelerated version */

#include <emmintrin.h>

typedef struct sfmt_128 {
  __m128i u;
} sfmt_128_t;

#else /* Use portable version */

typedef struct sfmt_128 {
  uint32_t u0, u1, u2, u3;
} sfmt_128_t;

#endif

/* Set the default random number generator to 11213.  This has the
   median period (all have practically inexhaustible periods with
   no risk for sequence overlap) and the median equidistribution
   (350 dim at 32-bit precision, 175 dim at 64-bit precision). */

#ifndef SFMT_E
#define SFMT_E 11213
#endif

enum sfmt_parameters {

  /* Parameter sets are defined by:

         SFMT_M   = ... < floor(SFMT_E/128) ...,
         SFMT_L1  = ... < 32 ...,
         SFMT_L2  = ... < 4  ...,
         SFMT_R1  = ... < 32 ...,
         SFMT_R2  = ... < 3  ...
       # define SFMT_MASK0   ((uint32_t)...)
       # define SFMT_MASK1   ((uint32_t)...)
       # define SFMT_MASK2   ((uint32_t)...)
       # define SFMT_MASK3   ((uint32_t)...)
       # define SFMT_PARITY0 ((uint32_t)...)
       # define SFMT_PARITY1 ((uint32_t)...)
       # define SFMT_PARITY2 ((uint32_t)...)
       # define SFMT_PARITY3 ((uint32_t)...)

     The masks and parities technically can't be part of the enum
     because enums hate potentially large unsigned quantities */

# if SFMT_E==607 /* Verified */

  SFMT_M =   2, SFMT_L1 = 15, SFMT_L2 = 3, SFMT_R1 = 13, SFMT_R2 = 3,
# define SFMT_MASK0   ((uint32_t)0xfdff37ff)
# define SFMT_MASK1   ((uint32_t)0xef7f3f7d)
# define SFMT_MASK2   ((uint32_t)0xff777b7d)
# define SFMT_MASK3   ((uint32_t)0x7ff7fb2f)
# define SFMT_PARITY0 ((uint32_t)0x00000001)
# define SFMT_PARITY1 ((uint32_t)0x00000000)
# define SFMT_PARITY2 ((uint32_t)0x00000000)
# define SFMT_PARITY3 ((uint32_t)0x5986f054)

# elif SFMT_E==1279 /* Verified */

  SFMT_M =   7, SFMT_L1 = 14, SFMT_L2 = 3, SFMT_R1 =  5, SFMT_R2 = 1,
# define SFMT_MASK0   ((uint32_t)0xf7fefffd)
# define SFMT_MASK1   ((uint32_t)0x7fefcfff)
# define SFMT_MASK2   ((uint32_t)0xaff3ef3f)
# define SFMT_MASK3   ((uint32_t)0xb5ffff7f)
# define SFMT_PARITY0 ((uint32_t)0x00000001)
# define SFMT_PARITY1 ((uint32_t)0x00000000)
# define SFMT_PARITY2 ((uint32_t)0x00000000)
# define SFMT_PARITY3 ((uint32_t)0x20000000)

# elif SFMT_E==2281 /* Verified */

  SFMT_M =  12, SFMT_L1 = 19, SFMT_L2 = 1, SFMT_R1 =  5, SFMT_R2 = 1,
# define SFMT_MASK0   ((uint32_t)0xbff7ffbf)
# define SFMT_MASK1   ((uint32_t)0xfdfffffe)
# define SFMT_MASK2   ((uint32_t)0xf7ffef7f)
# define SFMT_MASK3   ((uint32_t)0xf2f7cbbf)
# define SFMT_PARITY0 ((uint32_t)0x00000001)
# define SFMT_PARITY1 ((uint32_t)0x00000000)
# define SFMT_PARITY2 ((uint32_t)0x00000000)
# define SFMT_PARITY3 ((uint32_t)0x41dfa600)

# elif SFMT_E==4253 /* Verified */

  SFMT_M =  17, SFMT_L1 = 20, SFMT_L2 = 1, SFMT_R1 =  7, SFMT_R2 = 1,
# define SFMT_MASK0   ((uint32_t)0x9f7bffff)
# define SFMT_MASK1   ((uint32_t)0x9fffff5f)
# define SFMT_MASK2   ((uint32_t)0x3efffffb)
# define SFMT_MASK3   ((uint32_t)0xfffff7bb)
# define SFMT_PARITY0 ((uint32_t)0xa8000001)
# define SFMT_PARITY1 ((uint32_t)0xaf5390a3)
# define SFMT_PARITY2 ((uint32_t)0xb740b3f8)
# define SFMT_PARITY3 ((uint32_t)0x6c11486d)

# elif SFMT_E==11213 /* Verified */

  SFMT_M =  68, SFMT_L1 = 14, SFMT_L2 = 3, SFMT_R1 =  7, SFMT_R2 = 3,
# define SFMT_MASK0   ((uint32_t)0xeffff7fb)
# define SFMT_MASK1   ((uint32_t)0xffffffef)
# define SFMT_MASK2   ((uint32_t)0xdfdfbfff)
# define SFMT_MASK3   ((uint32_t)0x7fffdbfd)
# define SFMT_PARITY0 ((uint32_t)0x00000001)
# define SFMT_PARITY1 ((uint32_t)0x00000000)
# define SFMT_PARITY2 ((uint32_t)0xe8148000)
# define SFMT_PARITY3 ((uint32_t)0xd0c7afa3)

# elif SFMT_E==19937 /* Verified */

  SFMT_M = 122, SFMT_L1 = 18, SFMT_L2 = 1, SFMT_R1 = 11, SFMT_R2 = 1,
# define SFMT_MASK0   ((uint32_t)0xdfffffef)
# define SFMT_MASK1   ((uint32_t)0xddfecb7f)
# define SFMT_MASK2   ((uint32_t)0xbffaffff)
# define SFMT_MASK3   ((uint32_t)0xbffffff6)
# define SFMT_PARITY0 ((uint32_t)0x00000001)
# define SFMT_PARITY1 ((uint32_t)0x00000000)
# define SFMT_PARITY2 ((uint32_t)0x00000000)
# define SFMT_PARITY3 ((uint32_t)0x13c9e684)

# elif SFMT_E==44497 /* Verified */

  SFMT_M = 330, SFMT_L1 =  5, SFMT_L2 = 3, SFMT_R1 =  9, SFMT_R2 = 3,
# define SFMT_MASK0   ((uint32_t)0xeffffffb)
# define SFMT_MASK1   ((uint32_t)0xdfbebfff)
# define SFMT_MASK2   ((uint32_t)0xbfbf7bef)
# define SFMT_MASK3   ((uint32_t)0x9ffd7bff)
# define SFMT_PARITY0 ((uint32_t)0x00000001)
# define SFMT_PARITY1 ((uint32_t)0x00000000)
# define SFMT_PARITY2 ((uint32_t)0xa3ac4000)
# define SFMT_PARITY3 ((uint32_t)0xecc1327a)

  /* Note: SFMT_E==86243 not supported because SFMT_L2 is too large
     for this implementation */

# elif SFMT_E==132049 /* Verified */

  SFMT_M = 110, SFMT_L1 = 19, SFMT_L2 = 1, SFMT_R1 = 21, SFMT_R2 = 1,
# define SFMT_MASK0   ((uint32_t)0xffffbb5f)
# define SFMT_MASK1   ((uint32_t)0xfb6ebf95)
# define SFMT_MASK2   ((uint32_t)0xfffefffa)
# define SFMT_MASK3   ((uint32_t)0xcff77fff)
# define SFMT_PARITY0 ((uint32_t)0x00000001)
# define SFMT_PARITY1 ((uint32_t)0x00000000)
# define SFMT_PARITY2 ((uint32_t)0xcb520000)
# define SFMT_PARITY3 ((uint32_t)0xc7e91c7d)

# elif SFMT_E==216091 /* Verified */

  SFMT_M = 627, SFMT_L1 = 11, SFMT_L2 = 3, SFMT_R1 = 10, SFMT_R2 = 1,
# define SFMT_MASK0   ((uint32_t)0xbff7bff7)
# define SFMT_MASK1   ((uint32_t)0xbfffffff)
# define SFMT_MASK2   ((uint32_t)0xbffffa7f)
# define SFMT_MASK3   ((uint32_t)0xffddfbfb)
# define SFMT_PARITY0 ((uint32_t)0xf8000001)
# define SFMT_PARITY1 ((uint32_t)0x89e80709)
# define SFMT_PARITY2 ((uint32_t)0x3bd2b64b)
# define SFMT_PARITY3 ((uint32_t)0x0c64b1e4)

# else
# error "Unsupported SFMT exponent"
# endif

  /* Some useful derived quantities */

  SFMT_N   = SFMT_E/128 + 1, /* Number of 128-bit vectors in state */
  SFMT_NM  = SFMT_N-SFMT_M, 
  SFMT_L2A = 8*SFMT_L2,
  SFMT_R2A = 8*SFMT_R2,
  SFMT_L2B = 32-SFMT_L2A,
  SFMT_R2B = 32-SFMT_R2A,

  SFMT_NC  = SFMT_N*sizeof(sfmt_128_t),
  SFMT_NH  = SFMT_NC/sizeof(unsigned short),
  SFMT_NI  = SFMT_NC/sizeof(unsigned int),
  SFMT_NL  = SFMT_NC/sizeof(unsigned long),
  SFMT_N8  = SFMT_NC/sizeof(uint8_t),
  SFMT_N16 = SFMT_NC/sizeof(uint16_t),
  SFMT_N32 = SFMT_NC/sizeof(uint32_t),
  SFMT_N64 = SFMT_NC/sizeof(uint64_t)
};

struct rng {
  union {
    sfmt_128_t sfmt[ SFMT_N ]; /* Actual randgen state */
    /* Other ways to index into the state */
    unsigned char  uc[ SFMT_NC ];
    unsigned short uh[ SFMT_NH ];
    unsigned int   ui[ SFMT_NI ];
    unsigned long  ul[ SFMT_NL ];
    uint8_t   u8[  SFMT_N8 ];
    uint16_t u16[ SFMT_N16 ];
    uint32_t u32[ SFMT_N32 ];
    uint64_t u64[ SFMT_N64 ];
  } state;
  uint32_t n;      /* Next unextracted byte */
  uint32_t pad[3]; /* 16-byte align */
};

#if defined(__SSE2__)

# define DECL_SFMT                                                 \
  __m128i a_u, mask = _mm_setr_epi32( SFMT_MASK0, SFMT_MASK1,      \
                                      SFMT_MASK2, SFMT_MASK3 )

# define SFMT( a, b, c, d )                                        \
  a_u = a.u;                                                       \
  a.u = _mm_xor_si128(   a_u, _mm_xor_si128(                       \
        _mm_xor_si128(   _mm_slli_si128( a_u, SFMT_L2 ),           \
          _mm_and_si128( _mm_srli_epi32( b.u, SFMT_R1 ), mask ) ), \
        _mm_xor_si128(   _mm_srli_si128( c.u, SFMT_R2 ),           \
                         _mm_slli_epi32( d.u, SFMT_L1 ) ) ) )

#else

# define DECL_SFMT                                                      \
  uint32_t x0, x1, x2, x3, y0, y1, y2, y3

# define SFMT( a, b, c, d )                                             \
  x0 = ( a.u0 << SFMT_L2A );                                            \
  x1 = ( a.u1 << SFMT_L2A ) | ( a.u0 >> SFMT_L2B );                     \
  x2 = ( a.u2 << SFMT_L2A ) | ( a.u1 >> SFMT_L2B );                     \
  x3 = ( a.u3 << SFMT_L2A ) | ( a.u2 >> SFMT_L2B );                     \
  y0 = ( c.u0 >> SFMT_R2A ) | ( c.u1 << SFMT_R2B );                     \
  y1 = ( c.u1 >> SFMT_R2A ) | ( c.u2 << SFMT_R2B );                     \
  y2 = ( c.u2 >> SFMT_R2A ) | ( c.u3 << SFMT_R2B );                     \
  y3 = ( c.u3 >> SFMT_R2A );                                            \
  a.u0 ^= (x0 ^ ((b.u0>>SFMT_R1)&SFMT_MASK0)) ^ (y0 ^ (d.u0<<SFMT_L1)); \
  a.u1 ^= (x1 ^ ((b.u1>>SFMT_R1)&SFMT_MASK1)) ^ (y1 ^ (d.u1<<SFMT_L1)); \
  a.u2 ^= (x2 ^ ((b.u2>>SFMT_R1)&SFMT_MASK2)) ^ (y2 ^ (d.u2<<SFMT_L1)); \
  a.u3 ^= (x3 ^ ((b.u3>>SFMT_R1)&SFMT_MASK3)) ^ (y3 ^ (d.u3<<SFMT_L1))

#endif

STATIC_INLINE void
sfmt_next( sfmt_128_t * RESTRICT sfmt ) {
  DECL_SFMT;
  int n;

  SFMT( sfmt[0], sfmt[ SFMT_M   ], sfmt[ SFMT_N-2 ], sfmt[ SFMT_N-1 ] );
  SFMT( sfmt[1], sfmt[ SFMT_M+1 ], sfmt[ SFMT_N-1 ], sfmt[ 0        ] );
  for( n=2; n<SFMT_NM; n++ ) {
    SFMT( sfmt[n], sfmt[n+SFMT_M],  sfmt[n-2], sfmt[n-1] );
  }
  for( ; n<SFMT_N; n++ ) {
    SFMT( sfmt[n], sfmt[n-SFMT_NM], sfmt[n-2], sfmt[n-1] );
  }
}

# undef SFMT
# undef DECL_SFMT

/* Note that SFMT_NC is sizeof(r->state.p[0]) aligned */

#define RNG_NEXT( a, t, r, p, rs ) do {                                   \
    uint32_t _n = ((r)->n +   ((uint32_t)sizeof((r)->state.p[0]))-1 ) &   \
      /**/                 (~(((uint32_t)sizeof((r)->state.p[0]))-1));    \
    if( _n >= SFMT_NC ) sfmt_next( (r)->state.sfmt ), _n = 0;             \
    (a) = ((r)->state.p[ _n/(uint32_t)sizeof((r)->state.p[0]) ] >> (rs)); \
    (r)->n =             _n+(uint32_t)sizeof((r)->state.p[0]);            \
  } while(0)

/* Integer to random floating point conversions

   These have very rigorous interpretations.  Consider the drand
   conversions.  Imagine a real random number is generated on [0,1).
   The result is then rounded to the 2^53+1 point uniform lattice that
   covers this interval (endpoints inclusive).  The c0 variant
   corresponds to rounding the real random down to the nearest lattice
   point (as such, the result 1 will never occur).  The c1 variant
   corresponds to rounding the real random up to the nearest lattice
   point (as such, the result 0 will never occur).  The c variant
   corresponds to rounding to even the real random (as such, both 0
   and 1 could occur).  The open variant corresponds to rounding a
   2^52+1 point uniform lattice and rounding to the midpoint of the
   lattice (we lose have a bit of precision because the smallest
   representable number less than 1 is 1-eps/2 and for this midpoint
   of the lattice, the lattice spacing over the whole interval must be
   2*(eps/2), rather than eps/2.
    
   Similar issues apply for the frand. */

#define conv_frand(u32)    ((((u32)>>9 )+0.5f     )*(2.f/16777216.f       ))
#define conv_frand_c0(u32) (( (u32)>>8            )*(1.f/16777216.f       ))
#define conv_frand_c1(u32) ((((u32)>>8 )+1        )*(1.f/16777216.f       ))
#define conv_frand_c(u32)  ((((u32)>>8 )+((u32)&1))*(1.f/16777216.f       ))

#define conv_drand(u64)    ((((u64)>>12)+0.5      )*(2. /9007199254740992.))
#define conv_drand_c0(u64) (( (u64)>>11           )*(1. /9007199254740992.))
#define conv_drand_c1(u64) ((((u64)>>11)+1        )*(1. /9007199254740992.))
#define conv_drand_c(u64)  ((((u64)>>11)+((u64)&1))*(1. /9007199254740992.))

#endif /* _rng_private_h_ */
