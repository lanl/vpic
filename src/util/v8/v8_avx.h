#ifndef _v8_avx_h_
#define _v8_avx_h_

#ifndef IN_v8_h
#error "Do not include v8_avx.h directly; use v8.h"
#endif

#include <immintrin.h>
#include <math.h>

#define V8_ACCELERATION
#define V8_AVX_ACCELERATION

#ifndef ALIGNED
#define ALIGNED( n )
#endif

#define ALWAYS_INLINE __attribute__( ( always_inline ) )

// Why does GNU not define this function?
// #ifdef __GNUC__
#ifndef __INTEL_COMPILER
#define _mm256_set_m128( va, vb )                                              \
    _mm256_insertf128_ps( _mm256_castps128_ps256( vb ), va, 1 )
#endif

namespace v8
{
class v8;
class v8int;
class v8float;

////////////////
// v8 base class

class v8
{
    friend class v8int;
    friend class v8float;

    // v8 miscellaneous friends

    friend inline int any( const v8& a ) ALWAYS_INLINE;
    friend inline int all( const v8& a ) ALWAYS_INLINE;

    template <int n>
    friend inline v8 splat( const v8& a ) ALWAYS_INLINE;

    template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
    friend inline v8 shuffle( const v8& a ) ALWAYS_INLINE;

    friend inline void swap( v8& a, v8& b ) ALWAYS_INLINE;
    friend inline void transpose( v8& a0, v8& a1, v8& a2, v8& a3, v8& a4,
                                  v8& a5, v8& a6, v8& a7 ) ALWAYS_INLINE;

    // v8int miscellaneous friends

    friend inline v8 czero( const v8int& c, const v8& a ) ALWAYS_INLINE;
    friend inline v8 notczero( const v8int& c, const v8& a ) ALWAYS_INLINE;
    friend inline v8 merge( const v8int& c, const v8& a,
                            const v8& b ) ALWAYS_INLINE;

    // v8 memory manipulation friends

    friend inline void load_8x1( const void* ALIGNED( 16 ) p,
                                 v8& a ) ALWAYS_INLINE;
    friend inline void store_8x1( const v8& a,
                                  void* ALIGNED( 16 ) p ) ALWAYS_INLINE;
    friend inline void stream_8x1( const v8& a,
                                   void* ALIGNED( 16 ) p ) ALWAYS_INLINE;
    friend inline void clear_8x1( void* ALIGNED( 16 ) dst ) ALWAYS_INLINE;
    friend inline void copy_8x1( void* ALIGNED( 16 ) dst,
                                 const void* ALIGNED( 16 ) src ) ALWAYS_INLINE;
    friend inline void swap_8x1( void* ALIGNED( 16 ) a,
                                 void* ALIGNED( 16 ) b ) ALWAYS_INLINE;

    // v8 transposed memory manipulation friends
    // Note: Half aligned values are permissible in the 8x2_tr variants.

    friend inline void load_8x1_tr( const void* a0, const void* a1,
                                    const void* a2, const void* a3,
                                    const void* a4, const void* a5,
                                    const void* a6, const void* a7,
                                    v8& a ) ALWAYS_INLINE;

    friend inline void
    load_8x2_tr( const void* ALIGNED( 8 ) a0, const void* ALIGNED( 8 ) a1,
                 const void* ALIGNED( 8 ) a2, const void* ALIGNED( 8 ) a3,
                 const void* ALIGNED( 8 ) a4, const void* ALIGNED( 8 ) a5,
                 const void* ALIGNED( 8 ) a6, const void* ALIGNED( 8 ) a7,
                 v8& a, v8& b ) ALWAYS_INLINE;

    friend inline void
    load_8x3_tr( const void* ALIGNED( 16 ) a0, const void* ALIGNED( 16 ) a1,
                 const void* ALIGNED( 16 ) a2, const void* ALIGNED( 16 ) a3,
                 const void* ALIGNED( 16 ) a4, const void* ALIGNED( 16 ) a5,
                 const void* ALIGNED( 16 ) a6, const void* ALIGNED( 16 ) a7,
                 v8& a, v8& b, v8& c ) ALWAYS_INLINE;

    friend inline void
    load_8x4_tr( const void* ALIGNED( 16 ) a0, const void* ALIGNED( 16 ) a1,
                 const void* ALIGNED( 16 ) a2, const void* ALIGNED( 16 ) a3,
                 const void* ALIGNED( 16 ) a4, const void* ALIGNED( 16 ) a5,
                 const void* ALIGNED( 16 ) a6, const void* ALIGNED( 16 ) a7,
                 v8& a, v8& b, v8& c, v8& d ) ALWAYS_INLINE;

    friend inline void
    load_8x8_tr( const void* ALIGNED( 16 ) a0, const void* ALIGNED( 16 ) a1,
                 const void* ALIGNED( 16 ) a2, const void* ALIGNED( 16 ) a3,
                 const void* ALIGNED( 16 ) a4, const void* ALIGNED( 16 ) a5,
                 const void* ALIGNED( 16 ) a6, const void* ALIGNED( 16 ) a7,
                 v8& a, v8& b, v8& c, v8& d, v8& e, v8& f, v8& g,
                 v8& h ) ALWAYS_INLINE;

    friend inline void store_8x1_tr( const v8& a, void* a0, void* a1, void* a2,
                                     void* a3, void* a4, void* a5, void* a6,
                                     void* a7 ) ALWAYS_INLINE;

    friend inline void
    store_8x2_tr( const v8& a, const v8& b, void* ALIGNED( 8 ) a0,
                  void* ALIGNED( 8 ) a1, void* ALIGNED( 8 ) a2,
                  void* ALIGNED( 8 ) a3, void* ALIGNED( 8 ) a4,
                  void* ALIGNED( 8 ) a5, void* ALIGNED( 8 ) a6,
                  void* ALIGNED( 8 ) a7 ) ALWAYS_INLINE;

    friend inline void
    store_8x3_tr( const v8& a, const v8& b, const v8& c, void* ALIGNED( 16 ) a0,
                  void* ALIGNED( 16 ) a1, void* ALIGNED( 16 ) a2,
                  void* ALIGNED( 16 ) a3, void* ALIGNED( 16 ) a4,
                  void* ALIGNED( 16 ) a5, void* ALIGNED( 16 ) a6,
                  void* ALIGNED( 16 ) a7 ) ALWAYS_INLINE;

    friend inline void store_8x4_tr(
        const v8& a, const v8& b, const v8& c, const v8& d,
        void* ALIGNED( 16 ) a0, void* ALIGNED( 16 ) a1, void* ALIGNED( 16 ) a2,
        void* ALIGNED( 16 ) a3, void* ALIGNED( 16 ) a4, void* ALIGNED( 16 ) a5,
        void* ALIGNED( 16 ) a6, void* ALIGNED( 16 ) a7 ) ALWAYS_INLINE;

    friend inline void store_8x8_tr(
        const v8& a, const v8& b, const v8& c, const v8& d, const v8& e,
        const v8& f, const v8& g, const v8& h, void* ALIGNED( 16 ) a0,
        void* ALIGNED( 16 ) a1, void* ALIGNED( 16 ) a2, void* ALIGNED( 16 ) a3,
        void* ALIGNED( 16 ) a4, void* ALIGNED( 16 ) a5, void* ALIGNED( 16 ) a6,
        void* ALIGNED( 16 ) a7 ) ALWAYS_INLINE;

  protected:
    union {
        int i[8];
        float f[8];
        __m256 v;
    };

  public:
    v8() {} // Default constructor

    v8( const v8& a ) // Copy constructor
    {
        v = a.v;
    }

    ~v8() {} // Default destructor
};

// v8 miscellaneous functions

inline int any( const v8& a )
{
    return a.i[0] || a.i[1] || a.i[2] || a.i[3] || a.i[4] || a.i[5] || a.i[6] ||
           a.i[7];
}

inline int all( const v8& a )
{
    return a.i[0] && a.i[1] && a.i[2] && a.i[3] && a.i[4] && a.i[5] && a.i[6] &&
           a.i[7];
}

template <int n>
inline v8 splat( const v8& a )
{
    v8 b;

    b.v = _mm256_set1_ps( a.v[n] );

    return b;
}

template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
inline v8 shuffle( const v8& a )
{
    v8 b;

    b.i[0] = a.i[i0];
    b.i[1] = a.i[i1];
    b.i[2] = a.i[i2];
    b.i[3] = a.i[i3];
    b.i[4] = a.i[i4];
    b.i[5] = a.i[i5];
    b.i[6] = a.i[i6];
    b.i[7] = a.i[i7];

    return b;
}

inline void swap( v8& a, v8& b )
{
    __m256 a_v = a.v;

    a.v = b.v;

    b.v = a_v;
}

inline void transpose( v8& a0, v8& a1, v8& a2, v8& a3, v8& a4, v8& a5, v8& a6,
                       v8& a7 )
{
    __m256 t0, t1, t2, t3, t4, t5, t6, t7;

    __m256 u0, u1, u2, u3, u4, u5, u6, u7;

    t0 = _mm256_unpacklo_ps( a0.v, a1.v );
    t1 = _mm256_unpackhi_ps( a0.v, a1.v );
    t2 = _mm256_unpacklo_ps( a2.v, a3.v );
    t3 = _mm256_unpackhi_ps( a2.v, a3.v );
    t4 = _mm256_unpacklo_ps( a4.v, a5.v );
    t5 = _mm256_unpackhi_ps( a4.v, a5.v );
    t6 = _mm256_unpacklo_ps( a6.v, a7.v );
    t7 = _mm256_unpackhi_ps( a6.v, a7.v );

    u0 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u1 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    u2 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u3 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    u4 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u5 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    u6 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u7 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 3, 2, 3, 2 ) );

    a0.v = _mm256_permute2f128_ps( u0, u4, 0x20 );
    a1.v = _mm256_permute2f128_ps( u1, u5, 0x20 );
    a2.v = _mm256_permute2f128_ps( u2, u6, 0x20 );
    a3.v = _mm256_permute2f128_ps( u3, u7, 0x20 );
    a4.v = _mm256_permute2f128_ps( u0, u4, 0x31 );
    a5.v = _mm256_permute2f128_ps( u1, u5, 0x31 );
    a6.v = _mm256_permute2f128_ps( u2, u6, 0x31 );
    a7.v = _mm256_permute2f128_ps( u3, u7, 0x31 );
}

// v8 memory manipulation functions

inline void load_8x1( const void* ALIGNED( 16 ) p, v8& a )
{
    a.i[0] = ( (const int* ALIGNED( 16 ))p )[0];
    a.i[1] = ( (const int* ALIGNED( 16 ))p )[1];
    a.i[2] = ( (const int* ALIGNED( 16 ))p )[2];
    a.i[3] = ( (const int* ALIGNED( 16 ))p )[3];
    a.i[4] = ( (const int* ALIGNED( 16 ))p )[4];
    a.i[5] = ( (const int* ALIGNED( 16 ))p )[5];
    a.i[6] = ( (const int* ALIGNED( 16 ))p )[6];
    a.i[7] = ( (const int* ALIGNED( 16 ))p )[7];
}

inline void store_8x1( const v8& a, void* ALIGNED( 16 ) p )
{
    ( (int* ALIGNED( 16 ))p )[0] = a.i[0];
    ( (int* ALIGNED( 16 ))p )[1] = a.i[1];
    ( (int* ALIGNED( 16 ))p )[2] = a.i[2];
    ( (int* ALIGNED( 16 ))p )[3] = a.i[3];
    ( (int* ALIGNED( 16 ))p )[4] = a.i[4];
    ( (int* ALIGNED( 16 ))p )[5] = a.i[5];
    ( (int* ALIGNED( 16 ))p )[6] = a.i[6];
    ( (int* ALIGNED( 16 ))p )[7] = a.i[7];
}

inline void stream_8x1( const v8& a, void* ALIGNED( 16 ) p )
{
    ( (int* ALIGNED( 16 ))p )[0] = a.i[0];
    ( (int* ALIGNED( 16 ))p )[1] = a.i[1];
    ( (int* ALIGNED( 16 ))p )[2] = a.i[2];
    ( (int* ALIGNED( 16 ))p )[3] = a.i[3];
    ( (int* ALIGNED( 16 ))p )[4] = a.i[4];
    ( (int* ALIGNED( 16 ))p )[5] = a.i[5];
    ( (int* ALIGNED( 16 ))p )[6] = a.i[6];
    ( (int* ALIGNED( 16 ))p )[7] = a.i[7];
}

inline void clear_8x1( void* ALIGNED( 16 ) p )
{
    ( (int* ALIGNED( 16 ))p )[0] = 0;
    ( (int* ALIGNED( 16 ))p )[1] = 0;
    ( (int* ALIGNED( 16 ))p )[2] = 0;
    ( (int* ALIGNED( 16 ))p )[3] = 0;
    ( (int* ALIGNED( 16 ))p )[4] = 0;
    ( (int* ALIGNED( 16 ))p )[5] = 0;
    ( (int* ALIGNED( 16 ))p )[6] = 0;
    ( (int* ALIGNED( 16 ))p )[7] = 0;
}

// FIXME: Ordering semantics
inline void copy_8x1( void* ALIGNED( 16 ) dst, const void* ALIGNED( 16 ) src )
{
    ( (int* ALIGNED( 16 ))dst )[0] = ( (const int* ALIGNED( 16 ))src )[0];
    ( (int* ALIGNED( 16 ))dst )[1] = ( (const int* ALIGNED( 16 ))src )[1];
    ( (int* ALIGNED( 16 ))dst )[2] = ( (const int* ALIGNED( 16 ))src )[2];
    ( (int* ALIGNED( 16 ))dst )[3] = ( (const int* ALIGNED( 16 ))src )[3];
    ( (int* ALIGNED( 16 ))dst )[4] = ( (const int* ALIGNED( 16 ))src )[4];
    ( (int* ALIGNED( 16 ))dst )[5] = ( (const int* ALIGNED( 16 ))src )[5];
    ( (int* ALIGNED( 16 ))dst )[6] = ( (const int* ALIGNED( 16 ))src )[6];
    ( (int* ALIGNED( 16 ))dst )[7] = ( (const int* ALIGNED( 16 ))src )[7];
}

inline void swap_8x1( void* ALIGNED( 16 ) a, void* ALIGNED( 16 ) b )
{
    int t;

    t = ( (int* ALIGNED( 16 ))a )[0];
    ( (int* ALIGNED( 16 ))a )[0] = ( (int* ALIGNED( 16 ))b )[0];
    ( (int* ALIGNED( 16 ))b )[0] = t;

    t = ( (int* ALIGNED( 16 ))a )[1];
    ( (int* ALIGNED( 16 ))a )[1] = ( (int* ALIGNED( 16 ))b )[1];
    ( (int* ALIGNED( 16 ))b )[1] = t;

    t = ( (int* ALIGNED( 16 ))a )[2];
    ( (int* ALIGNED( 16 ))a )[2] = ( (int* ALIGNED( 16 ))b )[2];
    ( (int* ALIGNED( 16 ))b )[2] = t;

    t = ( (int* ALIGNED( 16 ))a )[3];
    ( (int* ALIGNED( 16 ))a )[3] = ( (int* ALIGNED( 16 ))b )[3];
    ( (int* ALIGNED( 16 ))b )[3] = t;

    t = ( (int* ALIGNED( 16 ))a )[4];
    ( (int* ALIGNED( 16 ))a )[4] = ( (int* ALIGNED( 16 ))b )[4];
    ( (int* ALIGNED( 16 ))b )[4] = t;

    t = ( (int* ALIGNED( 16 ))a )[5];
    ( (int* ALIGNED( 16 ))a )[5] = ( (int* ALIGNED( 16 ))b )[5];
    ( (int* ALIGNED( 16 ))b )[5] = t;

    t = ( (int* ALIGNED( 16 ))a )[6];
    ( (int* ALIGNED( 16 ))a )[6] = ( (int* ALIGNED( 16 ))b )[6];
    ( (int* ALIGNED( 16 ))b )[6] = t;

    t = ( (int* ALIGNED( 16 ))a )[7];
    ( (int* ALIGNED( 16 ))a )[7] = ( (int* ALIGNED( 16 ))b )[7];
    ( (int* ALIGNED( 16 ))b )[7] = t;
}

// v8 transposed memory manipulation functions

inline void load_8x1_tr( const void* a0, const void* a1, const void* a2,
                         const void* a3, const void* a4, const void* a5,
                         const void* a6, const void* a7, v8& a )
{
    a.i[0] = ( (const int*)a0 )[0];
    a.i[1] = ( (const int*)a1 )[0];
    a.i[2] = ( (const int*)a2 )[0];
    a.i[3] = ( (const int*)a3 )[0];
    a.i[4] = ( (const int*)a4 )[0];
    a.i[5] = ( (const int*)a5 )[0];
    a.i[6] = ( (const int*)a6 )[0];
    a.i[7] = ( (const int*)a7 )[0];
}

inline void
load_8x2_tr( const void* ALIGNED( 8 ) a0, const void* ALIGNED( 8 ) a1,
             const void* ALIGNED( 8 ) a2, const void* ALIGNED( 8 ) a3,
             const void* ALIGNED( 8 ) a4, const void* ALIGNED( 8 ) a5,
             const void* ALIGNED( 8 ) a6, const void* ALIGNED( 8 ) a7, v8& a,
             v8& b )
{
    __m128 zero;
    __m128 t0, t1, t2, t3;
    __m256 u0, u1;

    zero = _mm_setzero_ps();

    t0 = _mm_loadh_pi( _mm_loadl_pi( zero, (__m64*)a0 ), (__m64*)a1 );
    t1 = _mm_loadh_pi( _mm_loadl_pi( zero, (__m64*)a2 ), (__m64*)a3 );
    t2 = _mm_loadh_pi( _mm_loadl_pi( zero, (__m64*)a4 ), (__m64*)a5 );
    t3 = _mm_loadh_pi( _mm_loadl_pi( zero, (__m64*)a6 ), (__m64*)a7 );

    u0 = _mm256_set_m128( t2, t0 );
    u1 = _mm256_set_m128( t3, t1 );

    a.v = _mm256_shuffle_ps( u0, u1, _MM_SHUFFLE( 2, 0, 2, 0 ) );
    b.v = _mm256_shuffle_ps( u0, u1, _MM_SHUFFLE( 3, 1, 3, 1 ) );
}

inline void
load_8x3_tr( const void* ALIGNED( 16 ) a0, const void* ALIGNED( 16 ) a1,
             const void* ALIGNED( 16 ) a2, const void* ALIGNED( 16 ) a3,
             const void* ALIGNED( 16 ) a4, const void* ALIGNED( 16 ) a5,
             const void* ALIGNED( 16 ) a6, const void* ALIGNED( 16 ) a7, v8& a,
             v8& b, v8& c )
{
    a.i[0] = ( (const int* ALIGNED( 16 ))a0 )[0];
    b.i[0] = ( (const int* ALIGNED( 16 ))a0 )[1];
    c.i[0] = ( (const int* ALIGNED( 16 ))a0 )[2];

    a.i[1] = ( (const int* ALIGNED( 16 ))a1 )[0];
    b.i[1] = ( (const int* ALIGNED( 16 ))a1 )[1];
    c.i[1] = ( (const int* ALIGNED( 16 ))a1 )[2];

    a.i[2] = ( (const int* ALIGNED( 16 ))a2 )[0];
    b.i[2] = ( (const int* ALIGNED( 16 ))a2 )[1];
    c.i[2] = ( (const int* ALIGNED( 16 ))a2 )[2];

    a.i[3] = ( (const int* ALIGNED( 16 ))a3 )[0];
    b.i[3] = ( (const int* ALIGNED( 16 ))a3 )[1];
    c.i[3] = ( (const int* ALIGNED( 16 ))a3 )[2];

    a.i[4] = ( (const int* ALIGNED( 16 ))a4 )[0];
    b.i[4] = ( (const int* ALIGNED( 16 ))a4 )[1];
    c.i[4] = ( (const int* ALIGNED( 16 ))a4 )[2];

    a.i[5] = ( (const int* ALIGNED( 16 ))a5 )[0];
    b.i[5] = ( (const int* ALIGNED( 16 ))a5 )[1];
    c.i[5] = ( (const int* ALIGNED( 16 ))a5 )[2];

    a.i[6] = ( (const int* ALIGNED( 16 ))a6 )[0];
    b.i[6] = ( (const int* ALIGNED( 16 ))a6 )[1];
    c.i[6] = ( (const int* ALIGNED( 16 ))a6 )[2];

    a.i[7] = ( (const int* ALIGNED( 16 ))a7 )[0];
    b.i[7] = ( (const int* ALIGNED( 16 ))a7 )[1];
    c.i[7] = ( (const int* ALIGNED( 16 ))a7 )[2];
}

inline void
load_8x4_tr( const void* ALIGNED( 16 ) a0, const void* ALIGNED( 16 ) a1,
             const void* ALIGNED( 16 ) a2, const void* ALIGNED( 16 ) a3,
             const void* ALIGNED( 16 ) a4, const void* ALIGNED( 16 ) a5,
             const void* ALIGNED( 16 ) a6, const void* ALIGNED( 16 ) a7, v8& a,
             v8& b, v8& c, v8& d )
{
    __m256 tmp0, tmp1, tmp2, tmp3;

    a.v = _mm256_set_m128( _mm_load_ps( (const float*)a4 ),
                           _mm_load_ps( (const float*)a0 ) );
    b.v = _mm256_set_m128( _mm_load_ps( (const float*)a5 ),
                           _mm_load_ps( (const float*)a1 ) );
    c.v = _mm256_set_m128( _mm_load_ps( (const float*)a6 ),
                           _mm_load_ps( (const float*)a2 ) );
    d.v = _mm256_set_m128( _mm_load_ps( (const float*)a7 ),
                           _mm_load_ps( (const float*)a3 ) );

    tmp0 = _mm256_shuffle_ps( a.v, b.v, 0x44 );
    tmp2 = _mm256_shuffle_ps( a.v, b.v, 0xEE );
    tmp1 = _mm256_shuffle_ps( c.v, d.v, 0x44 );
    tmp3 = _mm256_shuffle_ps( c.v, d.v, 0xEE );

    a.v = _mm256_shuffle_ps( tmp0, tmp1, 0x88 );
    b.v = _mm256_shuffle_ps( tmp0, tmp1, 0xDD );
    c.v = _mm256_shuffle_ps( tmp2, tmp3, 0x88 );
    d.v = _mm256_shuffle_ps( tmp2, tmp3, 0xDD );
}

inline void
load_8x8_tr( const void* ALIGNED( 16 ) a0, const void* ALIGNED( 16 ) a1,
             const void* ALIGNED( 16 ) a2, const void* ALIGNED( 16 ) a3,
             const void* ALIGNED( 16 ) a4, const void* ALIGNED( 16 ) a5,
             const void* ALIGNED( 16 ) a6, const void* ALIGNED( 16 ) a7, v8& a,
             v8& b, v8& c, v8& d, v8& e, v8& f, v8& g, v8& h )
{
    __m256 t0, t1, t2, t3, t4, t5, t6, t7;

    __m256 u0, u1, u2, u3, u4, u5, u6, u7;

    a.v = _mm256_load_ps( (const float*)a0 );
    b.v = _mm256_load_ps( (const float*)a1 );
    c.v = _mm256_load_ps( (const float*)a2 );
    d.v = _mm256_load_ps( (const float*)a3 );
    e.v = _mm256_load_ps( (const float*)a4 );
    f.v = _mm256_load_ps( (const float*)a5 );
    g.v = _mm256_load_ps( (const float*)a6 );
    h.v = _mm256_load_ps( (const float*)a7 );

    t0 = _mm256_unpacklo_ps( a.v, b.v );
    t1 = _mm256_unpackhi_ps( a.v, b.v );
    t2 = _mm256_unpacklo_ps( c.v, d.v );
    t3 = _mm256_unpackhi_ps( c.v, d.v );
    t4 = _mm256_unpacklo_ps( e.v, f.v );
    t5 = _mm256_unpackhi_ps( e.v, f.v );
    t6 = _mm256_unpacklo_ps( g.v, h.v );
    t7 = _mm256_unpackhi_ps( g.v, h.v );

    u0 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u1 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    u2 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u3 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    u4 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u5 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    u6 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u7 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 3, 2, 3, 2 ) );

    a.v = _mm256_permute2f128_ps( u0, u4, 0x20 );
    b.v = _mm256_permute2f128_ps( u1, u5, 0x20 );
    c.v = _mm256_permute2f128_ps( u2, u6, 0x20 );
    d.v = _mm256_permute2f128_ps( u3, u7, 0x20 );
    e.v = _mm256_permute2f128_ps( u0, u4, 0x31 );
    f.v = _mm256_permute2f128_ps( u1, u5, 0x31 );
    g.v = _mm256_permute2f128_ps( u2, u6, 0x31 );
    h.v = _mm256_permute2f128_ps( u3, u7, 0x31 );
}

inline void store_8x1_tr( const v8& a, void* a0, void* a1, void* a2, void* a3,
                          void* a4, void* a5, void* a6, void* a7 )
{
    ( (int*)a0 )[0] = a.i[0];
    ( (int*)a1 )[0] = a.i[1];
    ( (int*)a2 )[0] = a.i[2];
    ( (int*)a3 )[0] = a.i[3];
    ( (int*)a4 )[0] = a.i[4];
    ( (int*)a5 )[0] = a.i[5];
    ( (int*)a6 )[0] = a.i[6];
    ( (int*)a7 )[0] = a.i[7];
}

inline void store_8x2_tr( const v8& a, const v8& b, void* ALIGNED( 8 ) a0,
                          void* ALIGNED( 8 ) a1, void* ALIGNED( 8 ) a2,
                          void* ALIGNED( 8 ) a3, void* ALIGNED( 8 ) a4,
                          void* ALIGNED( 8 ) a5, void* ALIGNED( 8 ) a6,
                          void* ALIGNED( 8 ) a7 )
{
    __m256 u0, u1;
    __m128 t0, t1, t2, t3;

    u0 = _mm256_unpacklo_ps( a.v, b.v );
    u1 = _mm256_unpackhi_ps( a.v, b.v );

    t0 = _mm256_extractf128_ps( u0, 0 );
    t1 = _mm256_extractf128_ps( u1, 0 );
    t2 = _mm256_extractf128_ps( u0, 1 );
    t3 = _mm256_extractf128_ps( u1, 1 );

    _mm_storel_pi( (__m64*)a0, t0 );
    _mm_storeh_pi( (__m64*)a1, t0 );

    _mm_storel_pi( (__m64*)a2, t1 );
    _mm_storeh_pi( (__m64*)a3, t1 );

    _mm_storel_pi( (__m64*)a4, t2 );
    _mm_storeh_pi( (__m64*)a5, t2 );

    _mm_storel_pi( (__m64*)a6, t3 );
    _mm_storeh_pi( (__m64*)a7, t3 );
}

inline void store_8x3_tr( const v8& a, const v8& b, const v8& c,
                          void* ALIGNED( 16 ) a0, void* ALIGNED( 16 ) a1,
                          void* ALIGNED( 16 ) a2, void* ALIGNED( 16 ) a3,
                          void* ALIGNED( 16 ) a4, void* ALIGNED( 16 ) a5,
                          void* ALIGNED( 16 ) a6, void* ALIGNED( 16 ) a7 )
{
    ( (int* ALIGNED( 16 ))a0 )[0] = a.i[0];
    ( (int* ALIGNED( 16 ))a0 )[1] = b.i[0];
    ( (int* ALIGNED( 16 ))a0 )[2] = c.i[0];

    ( (int* ALIGNED( 16 ))a1 )[0] = a.i[1];
    ( (int* ALIGNED( 16 ))a1 )[1] = b.i[1];
    ( (int* ALIGNED( 16 ))a1 )[2] = c.i[1];

    ( (int* ALIGNED( 16 ))a2 )[0] = a.i[2];
    ( (int* ALIGNED( 16 ))a2 )[1] = b.i[2];
    ( (int* ALIGNED( 16 ))a2 )[2] = c.i[2];

    ( (int* ALIGNED( 16 ))a3 )[0] = a.i[3];
    ( (int* ALIGNED( 16 ))a3 )[1] = b.i[3];
    ( (int* ALIGNED( 16 ))a3 )[2] = c.i[3];

    ( (int* ALIGNED( 16 ))a4 )[0] = a.i[4];
    ( (int* ALIGNED( 16 ))a4 )[1] = b.i[4];
    ( (int* ALIGNED( 16 ))a4 )[2] = c.i[4];

    ( (int* ALIGNED( 16 ))a5 )[0] = a.i[5];
    ( (int* ALIGNED( 16 ))a5 )[1] = b.i[5];
    ( (int* ALIGNED( 16 ))a5 )[2] = c.i[5];

    ( (int* ALIGNED( 16 ))a6 )[0] = a.i[6];
    ( (int* ALIGNED( 16 ))a6 )[1] = b.i[6];
    ( (int* ALIGNED( 16 ))a6 )[2] = c.i[6];

    ( (int* ALIGNED( 16 ))a7 )[0] = a.i[7];
    ( (int* ALIGNED( 16 ))a7 )[1] = b.i[7];
    ( (int* ALIGNED( 16 ))a7 )[2] = c.i[7];
}

inline void store_8x4_tr( const v8& a, const v8& b, const v8& c, const v8& d,
                          void* ALIGNED( 16 ) a0, void* ALIGNED( 16 ) a1,
                          void* ALIGNED( 16 ) a2, void* ALIGNED( 16 ) a3,
                          void* ALIGNED( 16 ) a4, void* ALIGNED( 16 ) a5,
                          void* ALIGNED( 16 ) a6, void* ALIGNED( 16 ) a7 )
{
    __m256 u0, u1, u2, u3;
    __m256 t0, t1, t2, t3;
    __m128 s0, s1, s2, s3, s4, s5, s6, s7;

    u0 = _mm256_unpacklo_ps( a.v, b.v );
    u1 = _mm256_unpacklo_ps( c.v, d.v );
    u2 = _mm256_unpackhi_ps( a.v, b.v );
    u3 = _mm256_unpackhi_ps( c.v, d.v );

    t0 = _mm256_shuffle_ps( u0, u1, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    t1 = _mm256_shuffle_ps( u0, u1, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    t2 = _mm256_shuffle_ps( u2, u3, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    t3 = _mm256_shuffle_ps( u2, u3, _MM_SHUFFLE( 3, 2, 3, 2 ) );

    s0 = _mm256_extractf128_ps( t0, 0 );
    s1 = _mm256_extractf128_ps( t1, 0 );
    s2 = _mm256_extractf128_ps( t2, 0 );
    s3 = _mm256_extractf128_ps( t3, 0 );

    s4 = _mm256_extractf128_ps( t0, 1 );
    s5 = _mm256_extractf128_ps( t1, 1 );
    s6 = _mm256_extractf128_ps( t2, 1 );
    s7 = _mm256_extractf128_ps( t3, 1 );

    _mm_store_ps( (float*)a0, s0 );
    _mm_store_ps( (float*)a1, s1 );
    _mm_store_ps( (float*)a2, s2 );
    _mm_store_ps( (float*)a3, s3 );
    _mm_store_ps( (float*)a4, s4 );
    _mm_store_ps( (float*)a5, s5 );
    _mm_store_ps( (float*)a6, s6 );
    _mm_store_ps( (float*)a7, s7 );
}

inline void store_8x8_tr( const v8& a, const v8& b, const v8& c, const v8& d,
                          const v8& e, const v8& f, const v8& g, const v8& h,
                          void* ALIGNED( 16 ) a0, void* ALIGNED( 16 ) a1,
                          void* ALIGNED( 16 ) a2, void* ALIGNED( 16 ) a3,
                          void* ALIGNED( 16 ) a4, void* ALIGNED( 16 ) a5,
                          void* ALIGNED( 16 ) a6, void* ALIGNED( 16 ) a7 )
{
    __m256 t0, t1, t2, t3, t4, t5, t6, t7;

    __m256 u0, u1, u2, u3, u4, u5, u6, u7;

    t0 = _mm256_unpacklo_ps( a.v, b.v );
    t1 = _mm256_unpackhi_ps( a.v, b.v );
    t2 = _mm256_unpacklo_ps( c.v, d.v );
    t3 = _mm256_unpackhi_ps( c.v, d.v );
    t4 = _mm256_unpacklo_ps( e.v, f.v );
    t5 = _mm256_unpackhi_ps( e.v, f.v );
    t6 = _mm256_unpacklo_ps( g.v, h.v );
    t7 = _mm256_unpackhi_ps( g.v, h.v );

    u0 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u1 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    u2 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u3 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    u4 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u5 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 3, 2, 3, 2 ) );
    u6 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 1, 0, 1, 0 ) );
    u7 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 3, 2, 3, 2 ) );

    t0 = _mm256_permute2f128_ps( u0, u4, 0x20 );
    t1 = _mm256_permute2f128_ps( u1, u5, 0x20 );
    t2 = _mm256_permute2f128_ps( u2, u6, 0x20 );
    t3 = _mm256_permute2f128_ps( u3, u7, 0x20 );
    t4 = _mm256_permute2f128_ps( u0, u4, 0x31 );
    t5 = _mm256_permute2f128_ps( u1, u5, 0x31 );
    t6 = _mm256_permute2f128_ps( u2, u6, 0x31 );
    t7 = _mm256_permute2f128_ps( u3, u7, 0x31 );

    _mm256_store_ps( (float*)a0, t0 );
    _mm256_store_ps( (float*)a1, t1 );
    _mm256_store_ps( (float*)a2, t2 );
    _mm256_store_ps( (float*)a3, t3 );
    _mm256_store_ps( (float*)a4, t4 );
    _mm256_store_ps( (float*)a5, t5 );
    _mm256_store_ps( (float*)a6, t6 );
    _mm256_store_ps( (float*)a7, t7 );
}

//////////////
// v8int class

class v8int : public v8
{
    // v8int prefix unary operator friends

    friend inline v8int operator+( const v8int& a ) ALWAYS_INLINE;
    friend inline v8int operator-( const v8int& a ) ALWAYS_INLINE;
    friend inline v8int operator~( const v8int& a ) ALWAYS_INLINE;
    friend inline v8int operator!( const v8int& a ) ALWAYS_INLINE;
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v8int prefix increment / decrement operator friends

    friend inline v8int operator++( v8int& a ) ALWAYS_INLINE;
    friend inline v8int operator--( v8int& a ) ALWAYS_INLINE;

    // v8int postfix increment / decrement operator friends

    friend inline v8int operator++( v8int& a, int ) ALWAYS_INLINE;
    friend inline v8int operator--( v8int& a, int ) ALWAYS_INLINE;

    // v8int binary operator friends

    friend inline v8int operator+( const v8int& a,
                                   const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator-( const v8int& a,
                                   const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator*(const v8int& a, const v8int& b)ALWAYS_INLINE;
    friend inline v8int operator/( const v8int& a,
                                   const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator%( const v8int& a,
                                   const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator^( const v8int& a,
                                   const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator&(const v8int& a, const v8int& b)ALWAYS_INLINE;
    friend inline v8int operator|( const v8int& a,
                                   const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator<<( const v8int& a,
                                    const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator>>( const v8int& a,
                                    const v8int& b ) ALWAYS_INLINE;

    // v8int logical operator friends

    friend inline v8int operator<( const v8int& a,
                                   const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator>( const v8int& a,
                                   const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator==( const v8int& a,
                                    const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator!=( const v8int& a,
                                    const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator<=( const v8int& a,
                                    const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator>=( const v8int& a,
                                    const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator&&( const v8int& a,
                                    const v8int& b ) ALWAYS_INLINE;
    friend inline v8int operator||( const v8int& a,
                                    const v8int& b ) ALWAYS_INLINE;

    // v8int miscellaneous friends

    friend inline v8int abs( const v8int& a ) ALWAYS_INLINE;
    friend inline v8 czero( const v8int& c, const v8& a ) ALWAYS_INLINE;
    friend inline v8 notczero( const v8int& c, const v8& a ) ALWAYS_INLINE;
    // FIXME: cswap, notcswap!
    friend inline v8 merge( const v8int& c, const v8& t,
                            const v8& f ) ALWAYS_INLINE;

    // v8float unary operator friends

    friend inline v8int operator!( const v8float& a ) ALWAYS_INLINE;

    // v8float logical operator friends

    friend inline v8int operator<( const v8float& a,
                                   const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator>( const v8float& a,
                                   const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator==( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator!=( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator<=( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator>=( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator&&( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator||( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;

    // v8float miscellaneous friends

    friend inline v8float clear_bits( const v8int& m,
                                      const v8float& a ) ALWAYS_INLINE;
    friend inline v8float set_bits( const v8int& m,
                                    const v8float& a ) ALWAYS_INLINE;
    friend inline v8float toggle_bits( const v8int& m,
                                       const v8float& a ) ALWAYS_INLINE;

  public:
    // v8int constructors / destructors

    v8int() {} // Default constructor

    v8int( const v8int& a ) // Copy constructor
    {
        v = a.v;
    }

    v8int( const v8& a ) // Init from mixed
    {
        v = a.v;
    }

    v8int( int a ) // Init from scalar
    {
        union {
            int i;
            float f;
        } u;
        u.i = a;
        v = _mm256_set1_ps( u.f );
    }

    v8int( int i0, int i1, int i2, int i3, int i4, int i5, int i6,
           int i7 ) // Init from scalars
    {
        union {
            int i;
            float f;
        } u0, u1, u2, u3, u4, u5, u6, u7;

        u0.i = i0;
        u1.i = i1;
        u2.i = i2;
        u3.i = i3;
        u4.i = i4;
        u5.i = i5;
        u6.i = i6;
        u7.i = i7;

        v = _mm256_setr_ps( u0.f, u1.f, u2.f, u3.f, u4.f, u5.f, u6.f, u7.f );
    }

    ~v8int() {} // Destructor

    // v8int assignment operators

#define ASSIGN( op )                                                           \
    inline v8int& operator op( const v8int& b )                                \
    {                                                                          \
        i[0] op b.i[0];                                                        \
        i[1] op b.i[1];                                                        \
        i[2] op b.i[2];                                                        \
        i[3] op b.i[3];                                                        \
        i[4] op b.i[4];                                                        \
        i[5] op b.i[5];                                                        \
        i[6] op b.i[6];                                                        \
        i[7] op b.i[7];                                                        \
        return *this;                                                          \
    }

    inline v8int& operator=( const v8int& b )
    {
        v = b.v;
        return *this;
    }

    ASSIGN( += )
    ASSIGN( -= )
    ASSIGN( *= )
    ASSIGN( /= )
    ASSIGN( %= )

    inline v8int& operator^=( const v8int& b )
    {
        v = _mm256_xor_ps( v, b.v );
        return *this;
    }

    inline v8int& operator&=( const v8int& b )
    {
        v = _mm256_and_ps( v, b.v );
        return *this;
    }

    inline v8int& operator|=( const v8int& b )
    {
        v = _mm256_or_ps( v, b.v );
        return *this;
    }

    ASSIGN( <<= )
    ASSIGN( >>= )

#undef ASSIGN

    // v8int member access operator

    inline int& operator[]( int n ) { return i[n]; }

    inline int operator()( int n ) { return i[n]; }
};

// v8int prefix unary operators

#define PREFIX_UNARY( op )                                                     \
    inline v8int operator op( const v8int& a )                                 \
    {                                                                          \
        v8int b;                                                               \
        b.i[0] = ( op a.i[0] );                                                \
        b.i[1] = ( op a.i[1] );                                                \
        b.i[2] = ( op a.i[2] );                                                \
        b.i[3] = ( op a.i[3] );                                                \
        b.i[4] = ( op a.i[4] );                                                \
        b.i[5] = ( op a.i[5] );                                                \
        b.i[6] = ( op a.i[6] );                                                \
        b.i[7] = ( op a.i[7] );                                                \
        return b;                                                              \
    }

inline v8int operator+( const v8int& a )
{
    v8int b;

    b.v = a.v;

    return b;
}

PREFIX_UNARY( -)

inline v8int operator!( const v8int& a )
{
    v8int b;

    b.i[0] = -( !a.i[0] );
    b.i[1] = -( !a.i[1] );
    b.i[2] = -( !a.i[2] );
    b.i[3] = -( !a.i[3] );
    b.i[4] = -( !a.i[4] );
    b.i[5] = -( !a.i[5] );
    b.i[6] = -( !a.i[6] );
    b.i[7] = -( !a.i[7] );

    return b;
}

inline v8int operator~( const v8int& a )
{
    v8int b;

    union {
        int i;
        float f;
    } u;

    u.i = -1;

    b.v = _mm256_xor_ps( a.v, _mm256_set1_ps( u.f ) );

    return b;
}

#undef PREFIX_UNARY

// v8int prefix increment / decrement

#define PREFIX_INCDEC( op )                                                    \
    inline v8int operator op( v8int& a )                                       \
    {                                                                          \
        v8int b;                                                               \
        b.i[0] = ( op a.i[0] );                                                \
        b.i[1] = ( op a.i[1] );                                                \
        b.i[2] = ( op a.i[2] );                                                \
        b.i[3] = ( op a.i[3] );                                                \
        b.i[4] = ( op a.i[4] );                                                \
        b.i[5] = ( op a.i[5] );                                                \
        b.i[6] = ( op a.i[6] );                                                \
        b.i[7] = ( op a.i[7] );                                                \
        return b;                                                              \
    }

PREFIX_INCDEC( ++)
PREFIX_INCDEC( --)

#undef PREFIX_INCDEC

// v8int postfix increment / decrement

#define POSTFIX_INCDEC( op )                                                   \
    inline v8int operator op( v8int& a, int )                                  \
    {                                                                          \
        v8int b;                                                               \
        b.i[0] = ( a.i[0] op );                                                \
        b.i[1] = ( a.i[1] op );                                                \
        b.i[2] = ( a.i[2] op );                                                \
        b.i[3] = ( a.i[3] op );                                                \
        b.i[4] = ( a.i[4] op );                                                \
        b.i[5] = ( a.i[5] op );                                                \
        b.i[6] = ( a.i[6] op );                                                \
        b.i[7] = ( a.i[7] op );                                                \
        return b;                                                              \
    }

POSTFIX_INCDEC( ++)
POSTFIX_INCDEC( --)

#undef POSTFIX_INCDEC

// v8int binary operators

#define BINARY( op )                                                           \
    inline v8int operator op( const v8int& a, const v8int& b )                 \
    {                                                                          \
        v8int c;                                                               \
        c.i[0] = a.i[0] op b.i[0];                                             \
        c.i[1] = a.i[1] op b.i[1];                                             \
        c.i[2] = a.i[2] op b.i[2];                                             \
        c.i[3] = a.i[3] op b.i[3];                                             \
        c.i[4] = a.i[4] op b.i[4];                                             \
        c.i[5] = a.i[5] op b.i[5];                                             \
        c.i[6] = a.i[6] op b.i[6];                                             \
        c.i[7] = a.i[7] op b.i[7];                                             \
        return c;                                                              \
    }

BINARY( +)
BINARY( -)
BINARY( * )
BINARY( / )
BINARY( % )

inline v8int operator^( const v8int& a, const v8int& b )
{
    v8int c;

    c.v = _mm256_xor_ps( a.v, b.v );

    return c;
}

inline v8int operator&( const v8int& a, const v8int& b )
{
    v8int c;

    c.v = _mm256_and_ps( a.v, b.v );

    return c;
}

inline v8int operator|( const v8int& a, const v8int& b )
{
    v8int c;

    c.v = _mm256_or_ps( a.v, b.v );

    return c;
}

BINARY( << )
BINARY( >> )

#undef BINARY

// v8int logical operators

#define LOGICAL( op )                                                          \
    inline v8int operator op( const v8int& a, const v8int& b )                 \
    {                                                                          \
        v8int c;                                                               \
        c.i[0] = -( a.i[0] op b.i[0] );                                        \
        c.i[1] = -( a.i[1] op b.i[1] );                                        \
        c.i[2] = -( a.i[2] op b.i[2] );                                        \
        c.i[3] = -( a.i[3] op b.i[3] );                                        \
        c.i[4] = -( a.i[4] op b.i[4] );                                        \
        c.i[5] = -( a.i[5] op b.i[5] );                                        \
        c.i[6] = -( a.i[6] op b.i[6] );                                        \
        c.i[7] = -( a.i[7] op b.i[7] );                                        \
        return c;                                                              \
    }

LOGICAL( < )
LOGICAL( > )
LOGICAL( == )
LOGICAL( != )
LOGICAL( <= )
LOGICAL( >= )
LOGICAL( &&)
LOGICAL( || )

#undef LOGICAL

// v8int miscellaneous functions

inline v8int abs( const v8int& a )
{
    v8int b;

    b.i[0] = ( a.i[0] >= 0 ) ? a.i[0] : -a.i[0];
    b.i[1] = ( a.i[1] >= 0 ) ? a.i[1] : -a.i[1];
    b.i[2] = ( a.i[2] >= 0 ) ? a.i[2] : -a.i[2];
    b.i[3] = ( a.i[3] >= 0 ) ? a.i[3] : -a.i[3];
    b.i[4] = ( a.i[4] >= 0 ) ? a.i[4] : -a.i[4];
    b.i[5] = ( a.i[5] >= 0 ) ? a.i[5] : -a.i[5];
    b.i[6] = ( a.i[6] >= 0 ) ? a.i[6] : -a.i[6];
    b.i[7] = ( a.i[7] >= 0 ) ? a.i[7] : -a.i[7];

    return b;
}

inline v8 czero( const v8int& c, const v8& a )
{
    v8 b;

    b.v = _mm256_andnot_ps( c.v, a.v );

    return b;
}

inline v8 notczero( const v8int& c, const v8& a )
{
    v8 b;

    b.v = _mm256_and_ps( c.v, a.v );

    return b;
}

inline v8 merge( const v8int& c, const v8& t, const v8& f )
{
    __m256 c_v = c.v;

    v8 tf;

    tf.v =
        _mm256_or_ps( _mm256_andnot_ps( c_v, f.v ), _mm256_and_ps( c_v, t.v ) );

    return tf;
}

////////////////
// v8float class

class v8float : public v8
{
    // v8float prefix unary operator friends

    friend inline v8float operator+( const v8float& a ) ALWAYS_INLINE;
    friend inline v8float operator-( const v8float& a ) ALWAYS_INLINE;
    friend inline v8float operator~( const v8float& a ) ALWAYS_INLINE;
    friend inline v8int operator!( const v8float& a ) ALWAYS_INLINE;
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v8float prefix increment / decrement operator friends

    friend inline v8float operator++( v8float& a ) ALWAYS_INLINE;
    friend inline v8float operator--( v8float& a ) ALWAYS_INLINE;

    // v8float postfix increment / decrement operator friends

    friend inline v8float operator++( v8float& a, int ) ALWAYS_INLINE;
    friend inline v8float operator--( v8float& a, int ) ALWAYS_INLINE;

    // v8float binary operator friends

    friend inline v8float operator+( const v8float& a,
                                     const v8float& b ) ALWAYS_INLINE;
    friend inline v8float operator-( const v8float& a,
                                     const v8float& b ) ALWAYS_INLINE;
    friend inline v8float operator*(const v8float& a,
                                    const v8float& b)ALWAYS_INLINE;
    friend inline v8float operator/( const v8float& a,
                                     const v8float& b ) ALWAYS_INLINE;

    // v8float logical operator friends

    friend inline v8int operator<( const v8float& a,
                                   const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator>( const v8float& a,
                                   const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator==( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator!=( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator<=( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator>=( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator&&( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;
    friend inline v8int operator||( const v8float& a,
                                    const v8float& b ) ALWAYS_INLINE;

    // v8float math library friends

#define CMATH_FR1( fn )                                                        \
    friend inline v8float fn( const v8float& a ) ALWAYS_INLINE
#define CMATH_FR2( fn )                                                        \
    friend inline v8float fn( const v8float& a, const v8float& b ) ALWAYS_INLINE

    CMATH_FR1( acos );
    CMATH_FR1( asin );
    CMATH_FR1( atan );
    CMATH_FR2( atan2 );
    CMATH_FR1( ceil );
    CMATH_FR1( cos );
    CMATH_FR1( cosh );
    CMATH_FR1( exp );
    CMATH_FR1( fabs );
    CMATH_FR1( floor );
    CMATH_FR2( fmod );
    CMATH_FR1( log );
    CMATH_FR1( log10 );
    CMATH_FR2( pow );
    CMATH_FR1( sin );
    CMATH_FR1( sinh );
    CMATH_FR1( sqrt );
    CMATH_FR1( tan );
    CMATH_FR1( tanh );

    CMATH_FR2( copysign );

#undef CMATH_FR1
#undef CMATH_FR2

    // v8float miscellaneous friends

    friend inline v8float rsqrt_approx( const v8float& a ) ALWAYS_INLINE;
    friend inline v8float rsqrt( const v8float& a ) ALWAYS_INLINE;
    friend inline v8float rcp_approx( const v8float& a ) ALWAYS_INLINE;
    friend inline v8float rcp( const v8float& a ) ALWAYS_INLINE;
    friend inline v8float fma( const v8float& a, const v8float& b,
                               const v8float& c ) ALWAYS_INLINE;
    friend inline v8float fms( const v8float& a, const v8float& b,
                               const v8float& c ) ALWAYS_INLINE;
    friend inline v8float fnms( const v8float& a, const v8float& b,
                                const v8float& c ) ALWAYS_INLINE;
    friend inline v8float clear_bits( const v8int& m,
                                      const v8float& a ) ALWAYS_INLINE;
    friend inline v8float set_bits( const v8int& m,
                                    const v8float& a ) ALWAYS_INLINE;
    friend inline v8float toggle_bits( const v8int& m,
                                       const v8float& a ) ALWAYS_INLINE;
    friend inline void increment_8x1( float* ALIGNED( 16 ) p,
                                      const v8float& a ) ALWAYS_INLINE;
    friend inline void decrement_8x1( float* ALIGNED( 16 ) p,
                                      const v8float& a ) ALWAYS_INLINE;
    friend inline void scale_8x1( float* ALIGNED( 16 ) p,
                                  const v8float& a ) ALWAYS_INLINE;

  public:
    // v8float constructors / destructors

    v8float() {} // Default constructor

    v8float( const v8float& a ) // Copy constructor
    {
        v = a.v;
    }

    v8float( const v8& a ) // Init from mixed
    {
        v = a.v;
    }

    v8float( float a ) // Init from scalar
    {
        v = _mm256_set1_ps( a );
    }

    v8float( float f0, float f1, float f2, float f3, float f4, float f5,
             float f6, float f7 ) // Init from scalars
    {
        v = _mm256_setr_ps( f0, f1, f2, f3, f4, f5, f6, f7 );
    }

    ~v8float() {} // Destructor

    // v8float assignment operators

#define ASSIGN( op, intrin )                                                   \
    inline v8float& operator op( const v8float& b )                            \
    {                                                                          \
        v = intrin( v, b.v );                                                  \
        return *this;                                                          \
    }

    inline v8float& operator=( const v8float& b )
    {
        v = b.v;
        return *this;
    }

    ASSIGN( +=, _mm256_add_ps )
    ASSIGN( -=, _mm256_sub_ps )
    ASSIGN( *=, _mm256_mul_ps )
    ASSIGN( /=, _mm256_div_ps )

#undef ASSIGN

    // v8float member access operator

    inline float& operator[]( int n ) { return f[n]; }

    inline float operator()( int n ) { return f[n]; }
};

// v8float prefix unary operators

inline v8float operator+( const v8float& a )
{
    v8float b;

    b.v = a.v;

    return b;
}

inline v8float operator-( const v8float& a )
{
    v8float b;

    b.v = _mm256_sub_ps( _mm256_setzero_ps(), a.v );

    return b;
}

inline v8int operator!( const v8float& a )
{
    v8int b;

    b.v = _mm256_cmp_ps( _mm256_setzero_ps(), a.v, _CMP_EQ_OS );

    return b;
}

// v8float prefix increment / decrement operators

inline v8float operator++( v8float& a )
{
    v8float b;
    __m256 t = _mm256_add_ps( a.v, _mm256_set1_ps( 1 ) );

    a.v = t;
    b.v = t;

    return b;
}

inline v8float operator--( v8float& a )
{
    v8float b;
    __m256 t = _mm256_sub_ps( a.v, _mm256_set1_ps( 1 ) );

    a.v = t;
    b.v = t;

    return b;
}

// v8float postfix increment / decrement operators

inline v8float operator++( v8float& a, int )
{
    v8float b;
    __m256 a_v = a.v;

    a.v = _mm256_add_ps( a_v, _mm256_set1_ps( 1 ) );
    b.v = a_v;

    return b;
}

inline v8float operator--( v8float& a, int )
{
    v8float b;
    __m256 a_v = a.v;

    a.v = _mm256_sub_ps( a_v, _mm256_set1_ps( 1 ) );
    b.v = a_v;

    return b;
}

// v8float binary operators

#define BINARY( op, intrin )                                                   \
    inline v8float operator op( const v8float& a, const v8float& b )           \
    {                                                                          \
        v8float c;                                                             \
        c.v = intrin( a.v, b.v );                                              \
        return c;                                                              \
    }

BINARY( +, _mm256_add_ps )
BINARY( -, _mm256_sub_ps )
BINARY( *, _mm256_mul_ps )
BINARY( /, _mm256_div_ps )

#undef BINARY

// v8float logical operators

#define LOGICAL( op, intrin, flag )                                            \
    inline v8int operator op( const v8float& a, const v8float& b )             \
    {                                                                          \
        v8int c;                                                               \
        c.v = intrin( a.v, b.v, flag );                                        \
        return c;                                                              \
    }

LOGICAL( <, _mm256_cmp_ps, _CMP_LT_OS )
LOGICAL( >, _mm256_cmp_ps, _CMP_GT_OS )
LOGICAL( ==, _mm256_cmp_ps, _CMP_EQ_OS )
LOGICAL( !=, _mm256_cmp_ps, _CMP_NEQ_OS )
LOGICAL( <=, _mm256_cmp_ps, _CMP_LE_OS )
LOGICAL( >=, _mm256_cmp_ps, _CMP_GE_OS )

inline v8int operator&&( const v8float& a, const v8float& b )
{
    v8int c;
    __m256 vzero = _mm256_setzero_ps();

    c.v = _mm256_and_ps( _mm256_cmp_ps( a.v, vzero, _CMP_NEQ_OS ),
                         _mm256_cmp_ps( b.v, vzero, _CMP_NEQ_OS ) );

    return c;
}

inline v8int operator||( const v8float& a, const v8float& b )
{
    v8int c;
    __m256 vzero = _mm256_setzero_ps();

    c.v = _mm256_or_ps( _mm256_cmp_ps( a.v, vzero, _CMP_NEQ_OS ),
                        _mm256_cmp_ps( b.v, vzero, _CMP_NEQ_OS ) );

    return c;
}

#undef LOGICAL

// v8float math library functions

#define CMATH_FR1( fn )                                                        \
    inline v8float fn( const v8float& a )                                      \
    {                                                                          \
        v8float b;                                                             \
        b.f[0] = ::fn( a.f[0] );                                               \
        b.f[1] = ::fn( a.f[1] );                                               \
        b.f[2] = ::fn( a.f[2] );                                               \
        b.f[3] = ::fn( a.f[3] );                                               \
        b.f[4] = ::fn( a.f[4] );                                               \
        b.f[5] = ::fn( a.f[5] );                                               \
        b.f[6] = ::fn( a.f[6] );                                               \
        b.f[7] = ::fn( a.f[7] );                                               \
        return b;                                                              \
    }

#define CMATH_FR2( fn )                                                        \
    inline v8float fn( const v8float& a, const v8float& b )                    \
    {                                                                          \
        v8float c;                                                             \
        c.f[0] = ::fn( a.f[0], b.f[0] );                                       \
        c.f[1] = ::fn( a.f[1], b.f[1] );                                       \
        c.f[2] = ::fn( a.f[2], b.f[2] );                                       \
        c.f[3] = ::fn( a.f[3], b.f[3] );                                       \
        c.f[4] = ::fn( a.f[4], b.f[4] );                                       \
        c.f[5] = ::fn( a.f[5], b.f[5] );                                       \
        c.f[6] = ::fn( a.f[6], b.f[6] );                                       \
        c.f[7] = ::fn( a.f[7], b.f[7] );                                       \
        return c;                                                              \
    }

CMATH_FR1( acos )
CMATH_FR1( asin )
CMATH_FR1( atan ) CMATH_FR2( atan2 ) CMATH_FR1( ceil ) CMATH_FR1( cos )
    CMATH_FR1( cosh ) CMATH_FR1( exp )
    /*CMATH_FR1(fabs)*/ CMATH_FR1( floor ) CMATH_FR2( fmod ) CMATH_FR1( log )
        CMATH_FR1( log10 ) CMATH_FR2( pow ) CMATH_FR1( sin ) CMATH_FR1( sinh )
    /*CMATH_FR1(sqrt)*/ CMATH_FR1( tan ) CMATH_FR1( tanh )

        inline v8float fabs( const v8float& a )
{
    v8float b;

    b.v = _mm256_andnot_ps( _mm256_set1_ps( -0.f ), a.v );

    return b;
}

inline v8float sqrt( const v8float& a )
{
    v8float b;

    b.v = _mm256_sqrt_ps( a.v );

    return b;
}

inline v8float copysign( const v8float& a, const v8float& b )
{
    v8float c;
    __m256 t = _mm256_set1_ps( -0.f );

    c.v = _mm256_or_ps( _mm256_and_ps( t, b.v ), _mm256_andnot_ps( t, a.v ) );

    return c;
}

#undef CMATH_FR1
#undef CMATH_FR2

// v8float miscellaneous functions

inline v8float rsqrt_approx( const v8float& a )
{
    v8float b;

    b.v = _mm256_rsqrt_ps( a.v );

    return b;
}

#if 0
  inline v8float rsqrt( const v8float &a )
  {
    v8float b;

    b.f[0] = ::sqrt( 1.0f / a.f[0] );
    b.f[1] = ::sqrt( 1.0f / a.f[1] );
    b.f[2] = ::sqrt( 1.0f / a.f[2] );
    b.f[3] = ::sqrt( 1.0f / a.f[3] );
    b.f[4] = ::sqrt( 1.0f / a.f[4] );
    b.f[5] = ::sqrt( 1.0f / a.f[5] );
    b.f[6] = ::sqrt( 1.0f / a.f[6] );
    b.f[7] = ::sqrt( 1.0f / a.f[7] );

    return b;
  }
#endif

inline v8float rsqrt( const v8float& a )
{
    v8float b;
    __m256 a_v = a.v, b_v;

    b_v = _mm256_rsqrt_ps( a_v );

    // Note: It is quicker to just call div_ps and sqrt_ps if more
    // refinement desired!
    b.v = _mm256_add_ps(
        b_v, _mm256_mul_ps(
                 _mm256_set1_ps( 0.5f ),
                 _mm256_sub_ps(
                     b_v, _mm256_mul_ps(
                              a_v, _mm256_mul_ps(
                                       b_v, _mm256_mul_ps( b_v, b_v ) ) ) ) ) );

    return b;
}

#if 0
  inline v8float rsqrt( const v8float &a )
  {
    v8float b;

    b.v = _mm256_div_ps( _mm256_set1_ps( 1.0f ), _mm256_sqrt_ps( a.v ) );

    return b;
  }
#endif

#if 0
  inline v8float rsqrt( const v8float &a )
  {
    v8float b;

    for( int j = 0; j < 8; j++ )
      b.f[j] = ::sqrt( 1.0f / a.f[j] );

    return b;
  }
#endif

inline v8float rcp_approx( const v8float& a )
{
    v8float b;

    b.v = _mm256_rcp_ps( a.v );

    return b;
}

#if 0
  inline v8float rcp( const v8float &a )
  {
    v8float b;

    b.f[0] = 1.0f / a.f[0];
    b.f[1] = 1.0f / a.f[1];
    b.f[2] = 1.0f / a.f[2];
    b.f[3] = 1.0f / a.f[3];
    b.f[4] = 1.0f / a.f[4];
    b.f[5] = 1.0f / a.f[5];
    b.f[6] = 1.0f / a.f[6];
    b.f[7] = 1.0f / a.f[7];

    return b;
  }
#endif

inline v8float rcp( const v8float& a )
{
    v8float b;
    __m256 a_v = a.v, b_v;

    b_v = _mm256_rcp_ps( a_v );
    b.v = _mm256_sub_ps( _mm256_add_ps( b_v, b_v ),
                         _mm256_mul_ps( a_v, _mm256_mul_ps( b_v, b_v ) ) );

    return b;
}

#if 0
  inline v8float rcp( const v8float &a )
  {
    v8float b;

    b.v = _mm256_div_ps( _mm256_set1_ps( 1.0f ), a.v );

    return b;
  }
#endif

inline v8float fma( const v8float& a, const v8float& b, const v8float& c )
{
    v8float d;

    d.v = _mm256_add_ps( _mm256_mul_ps( a.v, b.v ), c.v );

    // d.v = _mm256_fmadd_ps( a.v, b.v, c.v );

    return d;
}

inline v8float fms( const v8float& a, const v8float& b, const v8float& c )
{
    v8float d;

    d.v = _mm256_sub_ps( _mm256_mul_ps( a.v, b.v ), c.v );

    // d.v = _mm256_fmsub_ps( a.v, b.v, c.v );

    return d;
}

inline v8float fnms( const v8float& a, const v8float& b, const v8float& c )
{
    v8float d;

    d.v = _mm256_sub_ps( c.v, _mm256_mul_ps( a.v, b.v ) );

    // d.v = _mm256_fnmadd_ps( a.v, b.v, c.v );

    return d;
}

inline v8float clear_bits( const v8int& m, const v8float& a )
{
    v8float b;

    b.v = _mm256_andnot_ps( m.v, a.v );

    return b;
}

inline v8float set_bits( const v8int& m, const v8float& a )
{
    v8float b;

    b.v = _mm256_or_ps( m.v, a.v );

    return b;
}

inline v8float toggle_bits( const v8int& m, const v8float& a )
{
    v8float b;

    b.v = _mm256_xor_ps( m.v, a.v );

    return b;
}

inline void increment_8x1( float* ALIGNED( 16 ) p, const v8float& a )
{
    _mm256_store_ps( p, _mm256_add_ps( _mm256_load_ps( p ), a.v ) );
}

inline void decrement_8x1( float* ALIGNED( 16 ) p, const v8float& a )
{
    _mm256_store_ps( p, _mm256_sub_ps( _mm256_load_ps( p ), a.v ) );
}

inline void scale_8x1( float* ALIGNED( 16 ) p, const v8float& a )
{
    _mm256_store_ps( p, _mm256_mul_ps( _mm256_load_ps( p ), a.v ) );
}

} // namespace v8

#endif // _v8_avx_h_
