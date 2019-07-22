#ifndef _v4_neon_h_
#define _v4_neon_h_

#ifndef IN_v4_h
#error "Do not include v4_neon.h directly; use v4.h"
#endif

#include <arm_neon.h>
#include <math.h>

#define V4_ACCELERATION
#define V4_NEON_ACCELERATION

#ifndef ALIGNED
#define ALIGNED(n)
#endif

// This does not work with gcc 5.3.1 and the -fopenmp-simd
// flag.  Does not seem to work with -fopenmp either.  Not
// sure why.  It does work with the Intel compiler.  Need
// to try later versions of gcc.
// #define ALWAYS_VECTORIZE _Pragma( "omp simd" )

// #define ALWAYS_VECTORIZE _Pragma( "simd" )

#define ALWAYS_VECTORIZE \
  _Pragma( "simd" ) \
  _Pragma( "vector aligned" )

#define ALWAYS_INLINE __attribute__((always_inline))

namespace v4
{
  class v4;
  class v4int;
  class v4float;

  ////////////////
  // v4 base class

  class v4
  {
    friend class v4int;
    friend class v4float;

    // v4 miscellaneous friends

    friend inline int any( const v4 &a ) ALWAYS_INLINE;
    friend inline int all( const v4 &a ) ALWAYS_INLINE;

    template<int n>
    friend inline v4 splat( const v4 &a ) ALWAYS_INLINE;

    template<int i0, int i1, int i2, int i3>
    friend inline v4 shuffle( const v4 &a ) ALWAYS_INLINE;

    friend inline void swap( v4 &a, v4 &b ) ALWAYS_INLINE;
    friend inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 ) ALWAYS_INLINE;

    // v4int miscellaneous friends

    friend inline v4    czero( const v4int &c, const v4 &a ) ALWAYS_INLINE;
    friend inline v4 notczero( const v4int &c, const v4 &a ) ALWAYS_INLINE;
    friend inline v4 merge( const v4int &c, const v4 &a, const v4 &b ) ALWAYS_INLINE;

    // v4 memory manipulation friends

    friend inline void   load_4x1( const void * ALIGNED(16) p, v4 &a ) ALWAYS_INLINE;
    friend inline void  store_4x1( const v4 &a, void * ALIGNED(16) p ) ALWAYS_INLINE;
    friend inline void stream_4x1( const v4 &a, void * ALIGNED(16) p ) ALWAYS_INLINE;
    friend inline void  clear_4x1( void * ALIGNED(16) dst ) ALWAYS_INLINE;
    friend inline void   copy_4x1( void * ALIGNED(16) dst,
                                   const void * ALIGNED(16) src ) ALWAYS_INLINE;
    friend inline void   swap_4x1( void * ALIGNED(16) a, void * ALIGNED(16) b ) ALWAYS_INLINE;

    // v4 transposed memory manipulation friends

    friend inline void load_4x1_tr( const void *a0, const void *a1,
                                    const void *a2, const void *a3,
                                    v4 &a ) ALWAYS_INLINE;

    friend inline void load_4x2_tr( const void * ALIGNED(8) a0,
                                    const void * ALIGNED(8) a1,
                                    const void * ALIGNED(8) a2,
                                    const void * ALIGNED(8) a3,
                                    v4 &a, v4 &b ) ALWAYS_INLINE;

    friend inline void load_4x3_tr( const void * ALIGNED(16) a0,
                                    const void * ALIGNED(16) a1,
                                    const void * ALIGNED(16) a2,
                                    const void * ALIGNED(16) a3,
                                    v4 &a, v4 &b, v4 &c ) ALWAYS_INLINE;

    friend inline void load_4x4_tr( const void * ALIGNED(16) a0,
                                    const void * ALIGNED(16) a1,
                                    const void * ALIGNED(16) a2,
                                    const void * ALIGNED(16) a3,
                                    v4 &a, v4 &b, v4 &c, v4 &d ) ALWAYS_INLINE;

    friend inline void store_4x1_tr( const v4 &a,
                                     void *a0, void *a1, void *a2, void *a3 ) ALWAYS_INLINE;

    friend inline void store_4x2_tr( const v4 &a, const v4 &b,
                                     void * ALIGNED(8) a0,
                                     void * ALIGNED(8) a1,
                                     void * ALIGNED(8) a2,
                                     void * ALIGNED(8) a3 ) ALWAYS_INLINE;

    friend inline void store_4x3_tr( const v4 &a, const v4 &b, const v4 &c,
                                     void * ALIGNED(16) a0,
                                     void * ALIGNED(16) a1,
                                     void * ALIGNED(16) a2,
                                     void * ALIGNED(16) a3 ) ALWAYS_INLINE;

    friend inline void store_4x4_tr( const v4 &a, const v4 &b,
                                     const v4 &c, const v4 &d,
                                     void * ALIGNED(16) a0,
                                     void * ALIGNED(16) a1,
                                     void * ALIGNED(16) a2,
                                     void * ALIGNED(16) a3 ) ALWAYS_INLINE;

  protected:

    union
    {
      int         i[4];
      float       f[4];
      int32x4_t   vsi;
      uint32x4_t  vui;
      float32x4_t v;
    };

  public:

    v4() {}                    // Default constructor

    v4( const v4 &a )          // Copy constructor
    {
      v = a.v;

      // ALWAYS_VECTORIZE
      // for( int j = 0; j < 4; j++ )
      //   i[j] = a.i[j];
    }

    ~v4() {}                   // Default destructor
  };

  // v4 miscellaneous functions

  inline int any( const v4 &a )
  {
    return a.i[0] || a.i[1] || a.i[2] || a.i[3];
  }

  inline int all( const v4 &a )
  {
    return a.i[0] && a.i[1] && a.i[2] && a.i[3];
  }

  template<int n>
  inline v4 splat( const v4 & a )
  {
    v4 b;
    // __m128 a_v = a.v;

    // b.v = _mm_shuffle_ps( a_v, a_v, ( n*permute<1,1,1,1>::value ) );


    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      b.i[j] = a.i[n];

    return b;
  }

  template<int i0, int i1, int i2, int i3>
  inline v4 shuffle( const v4 & a )
  {
    v4 b;
    // __m128 a_v = a.v;

    // b.v = _mm_shuffle_ps( a_v, a_v, ( permute<i0,i1,i2,i3>::value ) );

    b.i[0] = a.i[i0];
    b.i[1] = a.i[i1];
    b.i[2] = a.i[i2];
    b.i[3] = a.i[i3];

    return b;
  }

  #define sw(x,y) x^=y, y^=x, x^=y

  inline void swap( v4 &a, v4 &b )
  {
    // __m128 a_v = a.v;

    // a.v = b.v;

    // b.v = a_v;

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      sw( a.i[j], b.i[j] );
  }

  inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 )
  {
    float32x4_t r, s, t, u;

    r = vtrn1q_f32( a0.v, a1.v );
    s = vtrn2q_f32( a0.v, a1.v );

    t = vtrn1q_f32( a2.v, a3.v );
    u = vtrn2q_f32( a2.v, a3.v );

    a0.v = vtrn1q_f64( r, t );
    a2.v = vtrn2q_f64( r, t );

    a1.v = vtrn1q_f64( s, u );
    a3.v = vtrn2q_f64( s, u );
  }

  #if 0
  inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 )
  {
    float32x4_t a0_v, a2_v, t, u;

    //-----------------------------------------------------------------
    float32x2_t a0_vh = vget_high_f32( a0.v );
    float32x2_t a1_vh = vget_high_f32( a1.v );

    float32x2x2_t res_a0a1_h = vzip_f32( a0_vh, a1_vh );

    t = vcombine_f32( res_a0a1_h.val[0], res_a0a1_h.val[1] );
    //-----------------------------------------------------------------
    // t    = _mm_unpackhi_ps( a0.v, a1.v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t a0_vl = vget_low_f32( a0.v );
    float32x2_t a1_vl = vget_low_f32( a1.v );

    float32x2x2_t res_a0a1_l = vzip_f32( a0_vl, a1_vl );

    a0_v = vcombine_f32( res_a0a1_l.val[0], res_a0a1_l.val[1] );
    //-----------------------------------------------------------------
    // a0_v = _mm_unpacklo_ps( a0.v, a1.v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t a2_vh = vget_high_f32( a2.v );
    float32x2_t a3_vh = vget_high_f32( a3.v );

    float32x2x2_t res_a2a3_h = vzip_f32( a2_vh, a3_vh );

    u = vcombine_f32( res_a2a3_h.val[0], res_a2a3_h.val[1] );
    //-----------------------------------------------------------------
    // u    = _mm_unpackhi_ps( a2.v, a3.v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t a2_vl = vget_low_f32( a2.v );
    float32x2_t a3_vl = vget_low_f32( a3.v );

    float32x2x2_t res_a2a3_l = vzip_f32( a2_vl, a3_vl );

    a2_v = vcombine_f32( res_a2a3_l.val[0], res_a2a3_l.val[1] );
    //-----------------------------------------------------------------
    // a2_v = _mm_unpacklo_ps( a2.v, a3.v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    a0.v[0] = a0_v[0];
    a0.v[1] = a0_v[1];
    a0.v[2] = a2_v[0];
    a0.v[3] = a2_v[1];
    //-----------------------------------------------------------------
    // a0.v = _mm_movelh_ps( a0_v, a2_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    a1.v[0] = a0_v[2];
    a1.v[1] = a0_v[3];
    a1.v[2] = a2_v[2];
    a1.v[3] = a2_v[3];
    //-----------------------------------------------------------------
    // a1.v = _mm_movehl_ps( a2_v, a0_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    a2.v[0] = t[0];
    a2.v[1] = t[1];
    a2.v[2] = u[0];
    a2.v[3] = u[1];
    //-----------------------------------------------------------------
    // a2.v = _mm_movelh_ps( t, u );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    a3.v[0] = t[2];
    a3.v[1] = t[3];
    a3.v[2] = u[2];
    a3.v[3] = u[3];
    //-----------------------------------------------------------------
    // a3.v = _mm_movehl_ps( u, t );
    //-----------------------------------------------------------------

    // sw( a0.i[1],a1.i[0] ); sw( a0.i[2],a2.i[0] ); sw( a0.i[3],a3.i[0] );
    //                        sw( a1.i[2],a2.i[1] ); sw( a1.i[3],a3.i[1] );
    //                                               sw( a2.i[3],a3.i[2] );
  }
  #endif

  #if 0
  // Portable version.
  inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 )
  {
    sw( a0.i[1],a1.i[0] ); sw( a0.i[2],a2.i[0] ); sw( a0.i[3],a3.i[0] );
                           sw( a1.i[2],a2.i[1] ); sw( a1.i[3],a3.i[1] );
                                                  sw( a2.i[3],a3.i[2] );
  }
  #endif

  #undef sw

  // v4 memory manipulation functions

  inline void load_4x1( const void * ALIGNED(16) p,
                        v4 &a )
  {
    a.v = vld1q_f32( ( float * ) p );

    // a.v = _mm_load_ps( ( float * ) p );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   a.i[j] = ((const int * ALIGNED(16))p)[j];
  }

  inline void store_4x1( const v4 &a,
                         void * ALIGNED(16) p )
  {
    vst1q_f32( ( float * ) p, a.v );

    // _mm_store_ps( ( float * ) p, a.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   ((int * ALIGNED(16))p)[j] = a.i[j];
  }

  inline void stream_4x1( const v4 &a,
                          void * ALIGNED(16) p )
  {
    // _mm_stream_ps( ( float * ) p, a.v );

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      ((int * ALIGNED(16))p)[j] = a.i[j];
  }

  inline void clear_4x1( void * ALIGNED(16) p )
  {
    vst1q_f32( ( float * ) p, vdupq_n_f32( 0.0f ) );

    // _mm_store_ps( ( float * ) p, _mm_setzero_ps() );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   ((int * ALIGNED(16))p)[j] = 0;
  }

  // FIXME: Ordering semantics
  inline void copy_4x1( void * ALIGNED(16) dst,
                        const void * ALIGNED(16) src )
  {
    vst1q_f32( ( float * ) dst, vld1q_f32( ( const float * ) src ) );

    // _mm_store_ps( ( float * ) dst, _mm_load_ps( ( const float * ) src ) );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   ((int * ALIGNED(16))dst)[j] = ((const int * ALIGNED(16))src)[j];
  }

  inline void swap_4x1( void * ALIGNED(16) a,
                        void * ALIGNED(16) b )
  {
    float32x4_t t = vld1q_f32( ( float * ) a );

    vst1q_f32( ( float * ) a, vld1q_f32( ( float * ) b ) );
    vst1q_f32( ( float * ) b, t );

    // __m128 t = _mm_load_ps( ( float * ) a );

    // _mm_store_ps( ( float * ) a, _mm_load_ps( ( float * ) b ) );
    // _mm_store_ps( ( float * ) b, t );

    // int t;

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    // {
    //   t = ((int * ALIGNED(16))a)[j];
    //   ((int * ALIGNED(16))a)[j] = ((int * ALIGNED(16))b)[j];
    //   ((int * ALIGNED(16))b)[j] = t;
    // }
  }

  // v4 transposed memory manipulation functions

  inline void load_4x1_tr( const void *a0,
			   const void *a1,
                           const void *a2,
			   const void *a3,
                           v4 &a )
  {
    // a.v = _mm_setr_ps( ( (const float *) a0 )[0],
    //                    ( (const float *) a1 )[0],
    //                    ( (const float *) a2 )[0],
    //                    ( (const float *) a3 )[0] );

    // Not correct.
    // float32x4x4_t mat = vld4q_f32( (const float *) a0 );

    // a.v = mat.val[0];

    a.i[0] = ((const int *)a0)[0];
    a.i[1] = ((const int *)a1)[0];
    a.i[2] = ((const int *)a2)[0];
    a.i[3] = ((const int *)a3)[0];
  }

  inline void load_4x2_tr( const void * ALIGNED(8) a0,
                           const void * ALIGNED(8) a1,
                           const void * ALIGNED(8) a2,
                           const void * ALIGNED(8) a3,
                           v4 &a,
			   v4 &b )
  {
    // __m128 a_v, b_v, t;

    // b_v = _mm_setzero_ps();

    // t   = _mm_loadh_pi( _mm_loadl_pi( b_v, (__m64 *)a0 ), (__m64 *)a1 );
    // b_v = _mm_loadh_pi( _mm_loadl_pi( b_v, (__m64 *)a2 ), (__m64 *)a3 );

    // a_v = _mm_shuffle_ps( t, b_v, 0x88 );
    // b_v = _mm_shuffle_ps( t, b_v, 0xdd );

    // a.v = a_v;
    // b.v = b_v;

    // Not correct.
    // float32x4x4_t mat = vld4q_f32( (const float *) a0 );

    // a.v = mat.val[0];
    // b.v = mat.val[1];

    a.i[0] = ((const int * ALIGNED(8))a0)[0];
    b.i[0] = ((const int * ALIGNED(8))a0)[1];

    a.i[1] = ((const int * ALIGNED(8))a1)[0];
    b.i[1] = ((const int * ALIGNED(8))a1)[1];

    a.i[2] = ((const int * ALIGNED(8))a2)[0];
    b.i[2] = ((const int * ALIGNED(8))a2)[1];

    a.i[3] = ((const int * ALIGNED(8))a3)[0];
    b.i[3] = ((const int * ALIGNED(8))a3)[1];
  }

  inline void load_4x3_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
                           v4 &a,
			   v4 &b,
			   v4 &c )
  {
    // __m128 a_v, b_v, c_v, t, u;

    // t   = _mm_load_ps( (const float *)a0 );
    // b_v = _mm_load_ps( (const float *)a1 );
    // c_v = _mm_load_ps( (const float *)a2 );
    // u   = _mm_load_ps( (const float *)a3 );

    // a_v = _mm_unpacklo_ps( t, b_v );
    // b_v = _mm_unpackhi_ps( t, b_v );
    // t   = _mm_unpacklo_ps( c_v, u );
    // u   = _mm_unpackhi_ps( c_v, u );

    // c_v = _mm_movelh_ps( b_v, u );
    // b_v = _mm_movehl_ps( t, a_v );
    // a_v = _mm_movelh_ps( a_v, t );

    // a.v = a_v;
    // b.v = b_v;
    // c.v = c_v;

    // Not correct.
    // float32x4x4_t mat = vld4q_f32( (const float *) a0 );

    // a.v = mat.val[0];
    // b.v = mat.val[1];
    // c.v = mat.val[2];

    a.i[0] = ((const int * ALIGNED(16))a0)[0];
    b.i[0] = ((const int * ALIGNED(16))a0)[1];
    c.i[0] = ((const int * ALIGNED(16))a0)[2];

    a.i[1] = ((const int * ALIGNED(16))a1)[0];
    b.i[1] = ((const int * ALIGNED(16))a1)[1];
    c.i[1] = ((const int * ALIGNED(16))a1)[2];

    a.i[2] = ((const int * ALIGNED(16))a2)[0];
    b.i[2] = ((const int * ALIGNED(16))a2)[1];
    c.i[2] = ((const int * ALIGNED(16))a2)[2];

    a.i[3] = ((const int * ALIGNED(16))a3)[0];
    b.i[3] = ((const int * ALIGNED(16))a3)[1];
    c.i[3] = ((const int * ALIGNED(16))a3)[2]; 
  }

  #if 0
  inline void load_4x4_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
                           v4 &a,
			   v4 &b,
			   v4 &c,
			   v4 &d )
  {
    // Not correct.
    float32x4x4_t mat = vld4q_f32( (const float *) a0 );

    a.v = mat.val[0];
    b.v = mat.val[1];
    c.v = mat.val[2];
    d.v = mat.val[3];

    // a.i[0] = ((const int * ALIGNED(16))a0)[0];
    // b.i[0] = ((const int * ALIGNED(16))a0)[1];
    // c.i[0] = ((const int * ALIGNED(16))a0)[2];
    // d.i[0] = ((const int * ALIGNED(16))a0)[3];

    // a.i[1] = ((const int * ALIGNED(16))a1)[0];
    // b.i[1] = ((const int * ALIGNED(16))a1)[1];
    // c.i[1] = ((const int * ALIGNED(16))a1)[2];
    // d.i[1] = ((const int * ALIGNED(16))a1)[3];

    // a.i[2] = ((const int * ALIGNED(16))a2)[0];
    // b.i[2] = ((const int * ALIGNED(16))a2)[1];
    // c.i[2] = ((const int * ALIGNED(16))a2)[2];
    // d.i[2] = ((const int * ALIGNED(16))a2)[3];

    // a.i[3] = ((const int * ALIGNED(16))a3)[0];
    // b.i[3] = ((const int * ALIGNED(16))a3)[1];
    // c.i[3] = ((const int * ALIGNED(16))a3)[2];
    // d.i[3] = ((const int * ALIGNED(16))a3)[3];
  }
  #endif

  #if 1
  inline void load_4x4_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
                           v4 &a,
			   v4 &b,
			   v4 &c,
			   v4 &d )
  {
    //-----------------------------------------------------------------
    float32x4_t a_v, b_v, c_v, d_v, t, u;
    //-----------------------------------------------------------------
    // __m128 a_v, b_v, c_v, d_v, t, u;
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    a_v = vld1q_f32( (const float *) a0 );
    b_v = vld1q_f32( (const float *) a1 );
    c_v = vld1q_f32( (const float *) a2 );
    d_v = vld1q_f32( (const float *) a3 );
    //-----------------------------------------------------------------
    // a_v = _mm_load_ps( (const float *) a0 );
    // b_v = _mm_load_ps( (const float *) a1 );
    // c_v = _mm_load_ps( (const float *) a2 );
    // d_v = _mm_load_ps( (const float *) a3 );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t a_vh = vget_high_f32( a_v );
    float32x2_t b_vh = vget_high_f32( b_v );

    float32x2x2_t res_ab_h = vzip_f32( a_vh, b_vh );

    t = vcombine_f32( res_ab_h.val[0], res_ab_h.val[1] );
    //-----------------------------------------------------------------
    // t   = _mm_unpackhi_ps( a_v, b_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t c_vh = vget_high_f32( c_v );
    float32x2_t d_vh = vget_high_f32( d_v );

    float32x2x2_t res_cd_h = vzip_f32( c_vh, d_vh );

    u = vcombine_f32( res_cd_h.val[0], res_cd_h.val[1] );
    //-----------------------------------------------------------------
    // u   = _mm_unpackhi_ps( c_v, d_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t a_vl = vget_low_f32( a_v );
    float32x2_t b_vl = vget_low_f32( b_v );

    float32x2x2_t res_ab_l = vzip_f32( a_vl, b_vl );

    a_v = vcombine_f32( res_ab_l.val[0], res_ab_l.val[1] );
    //-----------------------------------------------------------------
    // a_v = _mm_unpacklo_ps( a_v, b_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t c_vl = vget_low_f32( c_v );
    float32x2_t d_vl = vget_low_f32( d_v );

    float32x2x2_t res_cd_l = vzip_f32( c_vl, d_vl );

    c_v = vcombine_f32( res_cd_l.val[0], res_cd_l.val[1] );
    //-----------------------------------------------------------------
    // c_v = _mm_unpacklo_ps( c_v, d_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    a.v[0] = a_v[0];
    a.v[1] = a_v[1];
    a.v[2] = c_v[0];
    a.v[3] = c_v[1];
    //-----------------------------------------------------------------
    // a.v = _mm_movelh_ps( a_v, c_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    c.v[0] = t[0];
    c.v[1] = t[1];
    c.v[2] = u[0];
    c.v[3] = u[1];
    //-----------------------------------------------------------------
    // c.v = _mm_movelh_ps( t, u );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    b.v[0] = a_v[2];
    b.v[1] = a_v[3];
    b.v[2] = c_v[2];
    b.v[3] = c_v[3];
    //-----------------------------------------------------------------
    // b.v = _mm_movehl_ps( c_v, a_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    d.v[0] = t[2];
    d.v[1] = t[3];
    d.v[2] = u[2];
    d.v[3] = u[3];
    //-----------------------------------------------------------------
    // d.v = _mm_movehl_ps( u, t );
    //-----------------------------------------------------------------

    // a.i[0] = ((const int * ALIGNED(16))a0)[0];
    // b.i[0] = ((const int * ALIGNED(16))a0)[1];
    // c.i[0] = ((const int * ALIGNED(16))a0)[2];
    // d.i[0] = ((const int * ALIGNED(16))a0)[3];

    // a.i[1] = ((const int * ALIGNED(16))a1)[0];
    // b.i[1] = ((const int * ALIGNED(16))a1)[1];
    // c.i[1] = ((const int * ALIGNED(16))a1)[2];
    // d.i[1] = ((const int * ALIGNED(16))a1)[3];

    // a.i[2] = ((const int * ALIGNED(16))a2)[0];
    // b.i[2] = ((const int * ALIGNED(16))a2)[1];
    // c.i[2] = ((const int * ALIGNED(16))a2)[2];
    // d.i[2] = ((const int * ALIGNED(16))a2)[3];

    // a.i[3] = ((const int * ALIGNED(16))a3)[0];
    // b.i[3] = ((const int * ALIGNED(16))a3)[1];
    // c.i[3] = ((const int * ALIGNED(16))a3)[2];
    // d.i[3] = ((const int * ALIGNED(16))a3)[3];
  }
  #endif

  inline void store_4x1_tr( const v4 &a,
                            void *a0,
			    void *a1,
                            void *a2,
			    void *a3 )
  {
    ( (int *) a0 )[0] = a.i[0];
    ( (int *) a1 )[0] = a.i[1];
    ( (int *) a2 )[0] = a.i[2];
    ( (int *) a3 )[0] = a.i[3];
  }

  inline void store_4x2_tr( const v4 &a,
			    const v4 &b,
                            void * ALIGNED(8) a0,
			    void * ALIGNED(8) a1,
                            void * ALIGNED(8) a2,
			    void * ALIGNED(8) a3 )
  {
    // __m128 a_v = a.v, b_v = b.v, t;

    // t = _mm_unpacklo_ps( a_v, b_v ); // a0 b0 a1 b1 -> t

    // _mm_storel_pi( (__m64 *)a0, t ); // a0 b0       -> a0
    // _mm_storeh_pi( (__m64 *)a1, t ); // a1 b1       -> a1

    // t = _mm_unpackhi_ps( a_v, b_v ); // a2 b2 a3 b3 -> t

    // _mm_storel_pi( (__m64 *)a2, t ); // a2 b2       -> a2
    // _mm_storeh_pi( (__m64 *)a3, t ); // a3 b3       -> a3

    ( ( int * ALIGNED(8) ) a0 )[0] = a.i[0];
    ( ( int * ALIGNED(8) ) a0 )[1] = b.i[0];

    ( ( int * ALIGNED(8) ) a1 )[0] = a.i[1];
    ( ( int * ALIGNED(8) ) a1 )[1] = b.i[1];

    ( ( int * ALIGNED(8) ) a2 )[0] = a.i[2];
    ( ( int * ALIGNED(8) ) a2 )[1] = b.i[2];

    ( ( int * ALIGNED(8) ) a3 )[0] = a.i[3];
    ( ( int * ALIGNED(8) ) a3 )[1] = b.i[3];
  }

  inline void store_4x3_tr( const v4 &a,
			    const v4 &b,
			    const v4 &c,
                            void * ALIGNED(16) a0,
			    void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2,
			    void * ALIGNED(16) a3 )
  {
    // __m128 a_v = a.v, b_v = b.v, t;

    // t = _mm_unpacklo_ps( a_v, b_v ); // a0 b0 a1 b1 -> t

    // _mm_storel_pi( (__m64 *)a0, t ); // a0 b0       -> a0
    // _mm_storeh_pi( (__m64 *)a1, t ); // a1 b1       -> a1

    // t = _mm_unpackhi_ps( a_v, b_v ); // a2 b2 a3 b3 -> t

    // _mm_storel_pi( (__m64 *)a2, t ); // a2 b2       -> a2
    // _mm_storeh_pi( (__m64 *)a3, t ); // a3 b3       -> a3

    // ((float *)a0)[2] = c.f[0];
    // ((float *)a1)[2] = c.f[1];
    // ((float *)a2)[2] = c.f[2];
    // ((float *)a3)[2] = c.f[3];

    ( ( int * ALIGNED(16) ) a0 )[0] = a.i[0];
    ( ( int * ALIGNED(16) ) a0 )[1] = b.i[0];
    ( ( int * ALIGNED(16) ) a0 )[2] = c.i[0];

    ( ( int * ALIGNED(16) ) a1 )[0] = a.i[1];
    ( ( int * ALIGNED(16) ) a1 )[1] = b.i[1];
    ( ( int * ALIGNED(16) ) a1 )[2] = c.i[1];

    ( ( int * ALIGNED(16) ) a2 )[0] = a.i[2];
    ( ( int * ALIGNED(16) ) a2 )[1] = b.i[2];
    ( ( int * ALIGNED(16) ) a2 )[2] = c.i[2];

    ( ( int * ALIGNED(16) ) a3 )[0] = a.i[3];
    ( ( int * ALIGNED(16) ) a3 )[1] = b.i[3];
    ( ( int * ALIGNED(16) ) a3 )[2] = c.i[3];
  }

  #if 0
  inline void store_4x4_tr( const v4 &a,
			    const v4 &b,
			    const v4 &c,
			    const v4 &d,
                            void * ALIGNED(16) a0,
			    void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2,
			    void * ALIGNED(16) a3 )
  {
    // Not correct.
    float32x4x4_t mat;

    mat.val[0] = a.v;
    mat.val[1] = b.v;
    mat.val[2] = c.v;
    mat.val[3] = d.v;

    vst4q_f32( (const float *) a0, mat );

    // ((int * ALIGNED(16))a0)[0] = a.i[0];
    // ((int * ALIGNED(16))a0)[1] = b.i[0];
    // ((int * ALIGNED(16))a0)[2] = c.i[0];
    // ((int * ALIGNED(16))a0)[3] = d.i[0];

    // ((int * ALIGNED(16))a1)[0] = a.i[1];
    // ((int * ALIGNED(16))a1)[1] = b.i[1];
    // ((int * ALIGNED(16))a1)[2] = c.i[1];
    // ((int * ALIGNED(16))a1)[3] = d.i[1];

    // ((int * ALIGNED(16))a2)[0] = a.i[2];
    // ((int * ALIGNED(16))a2)[1] = b.i[2];
    // ((int * ALIGNED(16))a2)[2] = c.i[2];
    // ((int * ALIGNED(16))a2)[3] = d.i[2];

    // ((int * ALIGNED(16))a3)[0] = a.i[3];
    // ((int * ALIGNED(16))a3)[1] = b.i[3];
    // ((int * ALIGNED(16))a3)[2] = c.i[3];
    // ((int * ALIGNED(16))a3)[3] = d.i[3];
  }
  #endif

  #if 1
  inline void store_4x4_tr( const v4 &a,
			    const v4 &b,
			    const v4 &c,
			    const v4 &d,
                            void * ALIGNED(16) a0,
			    void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2,
			    void * ALIGNED(16) a3 )
  {
    //-----------------------------------------------------------------
    float32x4_t a_v = a.v, b_v = b.v, c_v = c.v, d_v = d.v, t, u;
    //-----------------------------------------------------------------
    // __m128 a_v = a.v, b_v = b.v, c_v = c.v, d_v = d.v, t, u;
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t a_vh = vget_high_f32( a_v );
    float32x2_t b_vh = vget_high_f32( b_v );

    float32x2x2_t res_ab_h = vzip_f32( a_vh, b_vh );

    t = vcombine_f32( res_ab_h.val[0], res_ab_h.val[1] );
    //-----------------------------------------------------------------
    // t   = _mm_unpackhi_ps( a_v, b_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t a_vl = vget_low_f32( a_v );
    float32x2_t b_vl = vget_low_f32( b_v );

    float32x2x2_t res_ab_l = vzip_f32( a_vl, b_vl );

    a_v = vcombine_f32( res_ab_l.val[0], res_ab_l.val[1] );
    //-----------------------------------------------------------------
    // a_v = _mm_unpacklo_ps( a_v, b_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t c_vh = vget_high_f32( c_v );
    float32x2_t d_vh = vget_high_f32( d_v );

    float32x2x2_t res_cd_h = vzip_f32( c_vh, d_vh );

    u = vcombine_f32( res_cd_h.val[0], res_cd_h.val[1] );
    //-----------------------------------------------------------------
    // u   = _mm_unpackhi_ps( c_v, d_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    float32x2_t c_vl = vget_low_f32( c_v );
    float32x2_t d_vl = vget_low_f32( d_v );

    float32x2x2_t res_cd_l = vzip_f32( c_vl, d_vl );

    c_v = vcombine_f32( res_cd_l.val[0], res_cd_l.val[1] );
    //-----------------------------------------------------------------
    // c_v = _mm_unpacklo_ps( c_v, d_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    b_v[0] = a_v[2];
    b_v[1] = a_v[3];
    b_v[2] = c_v[2];
    b_v[3] = c_v[3];
    //-----------------------------------------------------------------
    // b_v = _mm_movehl_ps( c_v, a_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    a_v[0] = a_v[0];
    a_v[1] = a_v[1];
    a_v[2] = c_v[0];
    a_v[3] = c_v[1];
    //-----------------------------------------------------------------
    // a_v = _mm_movelh_ps( a_v, c_v );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    c_v[0] = t[0];
    c_v[1] = t[1];
    c_v[2] = u[0];
    c_v[3] = u[1];
    //-----------------------------------------------------------------
    // c_v = _mm_movelh_ps( t, u );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    d_v[0] = t[2];
    d_v[1] = t[3];
    d_v[2] = u[2];
    d_v[3] = u[3];
    //-----------------------------------------------------------------
    // d_v = _mm_movehl_ps( u, t );
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    vst1q_f32( (float *) a0, a_v );
    vst1q_f32( (float *) a1, b_v );
    vst1q_f32( (float *) a2, c_v );
    vst1q_f32( (float *) a3, d_v );
    //-----------------------------------------------------------------
    // _mm_store_ps( (float *) a0, a_v );
    // _mm_store_ps( (float *) a1, b_v );
    // _mm_store_ps( (float *) a2, c_v );
    // _mm_store_ps( (float *) a3, d_v );
    //-----------------------------------------------------------------

    // ((int * ALIGNED(16))a0)[0] = a.i[0];
    // ((int * ALIGNED(16))a0)[1] = b.i[0];
    // ((int * ALIGNED(16))a0)[2] = c.i[0];
    // ((int * ALIGNED(16))a0)[3] = d.i[0];

    // ((int * ALIGNED(16))a1)[0] = a.i[1];
    // ((int * ALIGNED(16))a1)[1] = b.i[1];
    // ((int * ALIGNED(16))a1)[2] = c.i[1];
    // ((int * ALIGNED(16))a1)[3] = d.i[1];

    // ((int * ALIGNED(16))a2)[0] = a.i[2];
    // ((int * ALIGNED(16))a2)[1] = b.i[2];
    // ((int * ALIGNED(16))a2)[2] = c.i[2];
    // ((int * ALIGNED(16))a2)[3] = d.i[2];

    // ((int * ALIGNED(16))a3)[0] = a.i[3];
    // ((int * ALIGNED(16))a3)[1] = b.i[3];
    // ((int * ALIGNED(16))a3)[2] = c.i[3];
    // ((int * ALIGNED(16))a3)[3] = d.i[3];
  }
  #endif

  //////////////
  // v4int class

  class v4int : public v4
  {
    // v4int prefix unary operator friends

    friend inline v4int operator  +( const v4int & a ) ALWAYS_INLINE;
    friend inline v4int operator  -( const v4int & a ) ALWAYS_INLINE;
    friend inline v4int operator  ~( const v4int & a ) ALWAYS_INLINE;
    friend inline v4int operator  !( const v4int & a ) ALWAYS_INLINE;
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v4int prefix increment / decrement operator friends

    friend inline v4int operator ++( v4int & a ) ALWAYS_INLINE;
    friend inline v4int operator --( v4int & a ) ALWAYS_INLINE;

    // v4int postfix increment / decrement operator friends

    friend inline v4int operator ++( v4int & a, int ) ALWAYS_INLINE;
    friend inline v4int operator --( v4int & a, int ) ALWAYS_INLINE;

    // v4int binary operator friends

    friend inline v4int operator  +( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator  -( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator  *( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator  /( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator  %( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator  ^( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator  &( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator  |( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator <<( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator >>( const v4int &a, const v4int &b ) ALWAYS_INLINE;

    // v4int logical operator friends

    friend inline v4int operator  <( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator  >( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator ==( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator !=( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator <=( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator >=( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator &&( const v4int &a, const v4int &b ) ALWAYS_INLINE;
    friend inline v4int operator ||( const v4int &a, const v4int &b ) ALWAYS_INLINE;

    // v4int miscellaneous friends

    friend inline v4int abs( const v4int &a ) ALWAYS_INLINE;
    friend inline v4    czero( const v4int &c, const v4 &a ) ALWAYS_INLINE;
    friend inline v4 notczero( const v4int &c, const v4 &a ) ALWAYS_INLINE;
    // FIXME: cswap, notcswap!
    friend inline v4 merge( const v4int &c, const v4 &t, const v4 &f ) ALWAYS_INLINE;

    // v4float unary operator friends

    friend inline v4int operator  !( const v4float & a ) ALWAYS_INLINE;

    // v4float logical operator friends

    friend inline v4int operator  <( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator  >( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator ==( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator !=( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator <=( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator >=( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator &&( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator ||( const v4float &a, const v4float &b ) ALWAYS_INLINE;

    // v4float miscellaneous friends

    friend inline v4float  clear_bits( const v4int &m, const v4float &a ) ALWAYS_INLINE;
    friend inline v4float    set_bits( const v4int &m, const v4float &a ) ALWAYS_INLINE;
    friend inline v4float toggle_bits( const v4int &m, const v4float &a ) ALWAYS_INLINE;

  public:

    // v4int constructors / destructors

    v4int() {}                                // Default constructor

    v4int( const v4int &a )                   // Copy constructor
    {
      v = a.v;

      // ALWAYS_VECTORIZE
      // for( int j = 0; j < 4; j++ )
      //   i[j] = a.i[j];
    }

    v4int( const v4 &a )                      // Init from mixed
    {
      v = a.v;

      // ALWAYS_VECTORIZE
      // for( int j = 0; j < 4; j++ )
      //   i[j] = a.i[j];
    }

    v4int( int a )                            // Init from scalar
    {
      union
      {
        int i;
        float f;
      } u;

      u.i = a;
      v   = vdupq_n_f32( u.f );

      // ALWAYS_VECTORIZE
      // for( int j = 0; j < 4; j++ )
      //   i[j] = a;
    }

    v4int( int i0, int i1, int i2, int i3 )   // Init from scalars
    {
      // union
      // {
      //   int i;
      //   float f;
      // } u0, u1, u2, u3;

      // u0.i = i0;
      // u1.i = i1;
      // u2.i = i2;
      // u3.i = i3;

      // v = _mm_setr_ps( u0.f, u1.f, u2.f, u3.f );

      i[0] = i0;
      i[1] = i1;
      i[2] = i2;
      i[3] = i3;
    }

    ~v4int() {}                               // Destructor

    // v4int assignment operators

    #define ASSIGN(op)                            \
    inline v4int &operator op( const v4int &b )   \
    {                                             \
      ALWAYS_VECTORIZE                            \
      for( int j = 0; j < 4; j++ )                \
        i[j] op b.i[j];                           \
      return *this;                               \
    }

    ASSIGN(+=)
    ASSIGN(-=)
    ASSIGN(*=)
    ASSIGN(/=)
    ASSIGN(%=)
    ASSIGN(<<=)
    ASSIGN(>>=)
    // ASSIGN( =)
    // ASSIGN(^=)
    // ASSIGN(&=)
    // ASSIGN(|=)

    inline v4int &operator =( const v4int &b )
    {
      v = b.v;

      return *this;
    }

    inline v4int &operator ^=( const v4int &b )
    {
      vsi = veorq_s32( vsi, b.vsi );

      return *this;
    }

    inline v4int &operator &=( const v4int &b )
    {
      vsi = vandq_s32( vsi, b.vsi );

      return *this;
    }

    inline v4int &operator |=( const v4int &b )
    {
      vsi = vorrq_s32( vsi, b.vsi );

      return *this;
    }

    #undef ASSIGN

    // v4int member access operator

    inline int &operator []( int n )
    {
      return i[n];
    }

    inline int  operator ()( int n )
    {
      return i[n];
    }
  };

  // v4int prefix unary operators

  #define PREFIX_UNARY(op)                      \
  inline v4int operator op( const v4int & a )   \
  {                                             \
    v4int b;                                    \
    ALWAYS_VECTORIZE                            \
    for( int j = 0; j < 4; j++ )                \
      b.i[j] = ( op a.i[j] );                   \
    return b;                                   \
  }

  PREFIX_UNARY(+)
  PREFIX_UNARY(-)

  inline v4int operator !( const v4int & a )
  {
    v4int b;

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      b.i[j] = - ( !a.i[j] );

    return b;
  }

  PREFIX_UNARY(~)

  #undef PREFIX_UNARY

  // v4int prefix increment / decrement

  #define PREFIX_INCDEC(op)                     \
  inline v4int operator op( v4int & a )         \
  {                                             \
    v4int b;                                    \
    ALWAYS_VECTORIZE                            \
    for( int j = 0; j < 4; j++ )                \
      b.i[j] = ( op a.i[j] );                   \
    return b;                                   \
  }

  PREFIX_INCDEC(++)
  PREFIX_INCDEC(--)

  #undef PREFIX_INCDEC

  // v4int postfix increment / decrement

  #define POSTFIX_INCDEC(op)                   \
  inline v4int operator op( v4int & a, int )   \
  {                                            \
    v4int b;                                   \
    ALWAYS_VECTORIZE                           \
    for( int j = 0; j < 4; j++ )               \
      b.i[j] = ( a.i[j] op );                  \
    return b;                                  \
  }

  POSTFIX_INCDEC(++)
  POSTFIX_INCDEC(--)

  #undef POSTFIX_INCDEC

  // v4int binary operators

  #define BINARY(op)                                            \
  inline v4int operator op( const v4int &a, const v4int &b )    \
  {                                                             \
    v4int c;                                                    \
    ALWAYS_VECTORIZE                                            \
    for( int j = 0; j < 4; j++ )                                \
      c.i[j] = a.i[j] op b.i[j];                                \
    return c;                                                   \
  }

  BINARY(+)
  BINARY(-)
  BINARY(*)
  BINARY(/)
  BINARY(%)
  BINARY(<<)
  BINARY(>>)
  // BINARY(^)
  // BINARY(&)
  // BINARY(|)

  inline v4int operator ^( const v4int &a, const v4int &b )
  {
    v4int c;

    c.vsi = veorq_s32( a.vsi, b.vsi );

    return c;
  }

  inline v4int operator &( const v4int &a, const v4int &b )
  {
    v4int c;

    c.vsi = vandq_s32( a.vsi, b.vsi );

    return c;
  }

  inline v4int operator |( const v4int &a, const v4int &b )
  {
    v4int c;

    c.vsi = vorrq_s32( a.vsi, b.vsi );

    return c;
  }

  #undef BINARY

  // v4int logical operators

  #define LOGICAL(op)                                          \
  inline v4int operator op( const v4int &a, const v4int &b )   \
  {                                                            \
    v4int c;                                                   \
    ALWAYS_VECTORIZE                                           \
    for( int j = 0; j < 4; j++ )                               \
      c.i[j] = - ( a.i[j] op b.i[j] );                         \
    return c;                                                  \
  }

  LOGICAL(<)
  LOGICAL(>)
  LOGICAL(==)
  LOGICAL(!=)
  LOGICAL(<=)
  LOGICAL(>=)
  LOGICAL(&&)
  LOGICAL(||)

  #undef LOGICAL

  // v4int miscellaneous functions

  inline v4int abs( const v4int &a )
  {
    v4int b;

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      b.i[j] = ( a.i[j] >= 0 ) ? a.i[j] : -a.i[j];

    return b;
  }

  inline v4 czero( const v4int &c, const v4 &a )
  {
    v4 b;

    b.vsi = vbicq_s32( c.vsi, a.vsi );

    // b.v = _mm_andnot_ps( c.v, a.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.i[j] = a.i[j] & ~c.i[j];

    return b;
  }

  inline v4 notczero( const v4int &c, const v4 &a )
  {
    v4 b;

    b.vsi = vandq_s32( c.vsi, a.vsi );

    // b.v = _mm_and_ps( c.v, a.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.i[j] = a.i[j] & c.i[j];

    return b;
  }

  inline v4 merge( const v4int &c, const v4 &t, const v4 &f )
  {
    v4 m;

    m.vsi = vorrq_s32( vbicq_s32( c.vsi, f.vsi ),
                       vandq_s32( c.vsi, t.vsi ) );

    // __m128 c_v = c.v;
    // v4 tf;

    // tf.v = _mm_or_ps( _mm_andnot_ps( c_v, f.v ),
    //                   _mm_and_ps( c_v, t.v ) );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   m.i[j] = ( f.i[j] & ~c.i[j] ) | ( t.i[j] & c.i[j] );

    return m;
  }

  ////////////////
  // v4float class

  class v4float : public v4
  {
    // v4float prefix unary operator friends

    friend inline v4float operator  +( const v4float &a ) ALWAYS_INLINE;
    friend inline v4float operator  -( const v4float &a ) ALWAYS_INLINE;
    friend inline v4float operator  ~( const v4float &a ) ALWAYS_INLINE;
    friend inline v4int   operator  !( const v4float &a ) ALWAYS_INLINE;
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v4float prefix increment / decrement operator friends

    friend inline v4float operator ++( v4float &a ) ALWAYS_INLINE;
    friend inline v4float operator --( v4float &a ) ALWAYS_INLINE;

    // v4float postfix increment / decrement operator friends

    friend inline v4float operator ++( v4float &a, int ) ALWAYS_INLINE;
    friend inline v4float operator --( v4float &a, int ) ALWAYS_INLINE;

    // v4float binary operator friends

    friend inline v4float operator  +( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4float operator  -( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4float operator  *( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4float operator  /( const v4float &a, const v4float &b ) ALWAYS_INLINE;

    // v4float logical operator friends

    friend inline v4int operator  <( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator  >( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator ==( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator !=( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator <=( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator >=( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator &&( const v4float &a, const v4float &b ) ALWAYS_INLINE;
    friend inline v4int operator ||( const v4float &a, const v4float &b ) ALWAYS_INLINE;

    // v4float math library friends

    #define CMATH_FR1(fn) friend inline v4float fn( const v4float &a ) ALWAYS_INLINE
    #define CMATH_FR2(fn) friend inline v4float fn( const v4float &a,  \
                                                   const v4float &b ) ALWAYS_INLINE

    CMATH_FR1(acos);  CMATH_FR1(asin);  CMATH_FR1(atan); CMATH_FR2(atan2);
    CMATH_FR1(ceil);  CMATH_FR1(cos);   CMATH_FR1(cosh); CMATH_FR1(exp);
    CMATH_FR1(fabs);  CMATH_FR1(floor); CMATH_FR2(fmod); CMATH_FR1(log);
    CMATH_FR1(log10); CMATH_FR2(pow);   CMATH_FR1(sin);  CMATH_FR1(sinh);
    CMATH_FR1(sqrt);  CMATH_FR1(tan);   CMATH_FR1(tanh);

    CMATH_FR2(copysign);

    #undef CMATH_FR1
    #undef CMATH_FR2

    // v4float miscellaneous friends

    friend inline v4float rsqrt_approx( const v4float &a ) ALWAYS_INLINE;
    friend inline v4float rsqrt       ( const v4float &a ) ALWAYS_INLINE;
    friend inline v4float rcp_approx( const v4float &a ) ALWAYS_INLINE;
    friend inline v4float rcp       ( const v4float &a ) ALWAYS_INLINE;
    friend inline v4float fma ( const v4float &a, const v4float &b, const v4float &c ) ALWAYS_INLINE;
    friend inline v4float fms ( const v4float &a, const v4float &b, const v4float &c ) ALWAYS_INLINE;
    friend inline v4float fnms( const v4float &a, const v4float &b, const v4float &c ) ALWAYS_INLINE;
    friend inline v4float  clear_bits( const v4int &m, const v4float &a ) ALWAYS_INLINE;
    friend inline v4float    set_bits( const v4int &m, const v4float &a ) ALWAYS_INLINE;
    friend inline v4float toggle_bits( const v4int &m, const v4float &a ) ALWAYS_INLINE;
    friend inline void increment_4x1( float * ALIGNED(16) p, const v4float &a ) ALWAYS_INLINE;
    friend inline void decrement_4x1( float * ALIGNED(16) p, const v4float &a ) ALWAYS_INLINE;
    friend inline void     scale_4x1( float * ALIGNED(16) p, const v4float &a ) ALWAYS_INLINE;
    friend inline void trilinear( v4float &wl, v4float &wh ) ALWAYS_INLINE;

  public:

    // v4float constructors / destructors

    v4float() {}                                        // Default constructor

    v4float( const v4float &a )                         // Copy constructor
    {
      v = a.v;

      // ALWAYS_VECTORIZE
      // for( int j = 0; j < 4; j++ )
      //   f[j] = a.f[j];
    }

    v4float( const v4 &a )                              // Init from mixed
    {
      v = a.v;

      // ALWAYS_VECTORIZE
      // for( int j = 0; j < 4; j++ )
      //   f[j] = a.f[j];
    }

    v4float( float a )                                  // Init from scalar
    {
      v = vdupq_n_f32( a );

      // ALWAYS_VECTORIZE
      // for( int j = 0; j < 4; j++ )
      //   f[j] = a;
    }

    v4float( float f0, float f1, float f2, float f3 )   // Init from scalars
    {
      // v = _mm_setr_ps( f0, f1, f2, f3 );

      f[0] = f0;
      f[1] = f1;
      f[2] = f2;
      f[3] = f3;
    }

    ~v4float() {}                                       // Destructor

    // v4float assignment operators

    #define ASSIGN(op,intrin)                           \
    inline v4float &operator op( const v4float &b )     \
    {                                                   \
      v = intrin( v, b.v );                             \
      return *this;                                     \
    }

    ASSIGN( +=, vaddq_f32 )
    ASSIGN( -=, vsubq_f32 )
    ASSIGN( *=, vmulq_f32 )
    ASSIGN( /=, vdivq_f32 )

    inline v4float &operator =( const v4float &b )
    {
      v = b.v;

      return *this;
    }

    #undef ASSIGN

    // #define ASSIGN(op)                                  \
    // inline v4float &operator op( const v4float &b )     \
    // {                                                   \
    //   ALWAYS_VECTORIZE                                  \
    //   for( int j = 0; j < 4; j++ )                      \
    //     f[j] op b.f[j];                                 \
    //   return *this;                                     \
    // }

    // ASSIGN(=)
    // ASSIGN(+=)
    // ASSIGN(-=)
    // ASSIGN(*=)
    // ASSIGN(/=)

    // #undef ASSIGN

    // v4float member access operator

    inline float &operator []( int n )
    {
      return f[n];
    }

    inline float  operator ()( int n )
    {
      return f[n];
    }
  };

  // v4float prefix unary operators

  inline v4float operator +( const v4float &a )
  {
    v4float b;

    b.v = a.v;

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.f[j] = +a.f[j];

    return b;
  }

  inline v4float operator -( const v4float &a )
  {
    v4float b;

    // b.v = _mm_sub_ps( _mm_setzero_ps(), a.v );

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      b.f[j] = -a.f[j];

    return b;
  }

  inline v4int operator !( const v4float &a )
  {
    v4int b;

    // b.v = _mm_cmpeq_ps( _mm_setzero_ps(), a.v );

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      b.i[j] = a.i[j] ? 0 : -1;

    return b;
  }

  // v4float prefix increment / decrement operators

  inline v4float operator ++( v4float &a )
  {
    v4float b;

    // __m128 t = _mm_add_ps( a.v, _mm_set1_ps( 1 ) );

    // a.v = t;
    // b.v = t;

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      b.f[j] = ++a.f[j];

    return b;
  }

  inline v4float operator --( v4float &a )
  {
    v4float b;

    // __m128 t = _mm_sub_ps( a.v, _mm_set1_ps( 1 ) );

    // a.v = t;
    // b.v = t;

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      b.f[j] = --a.f[j];

    return b;
  }

  // v4float postfix increment / decrement operators

  inline v4float operator ++( v4float &a, int )
  {
    v4float b;

    // __m128 a_v = a.v;

    // a.v = _mm_add_ps( a_v, _mm_set1_ps( 1 ) );
    // b.v = a_v;

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      b.f[j] = a.f[j]++;

    return b;
  }

  inline v4float operator --( v4float &a, int )
  {
    v4float b;

    // __m128 a_v = a.v;

    // a.v = _mm_sub_ps( a_v, _mm_set1_ps( 1 ) );
    // b.v = a_v;

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
      b.f[j] = a.f[j]--;

    return b;
  }

  // v4float binary operators

  #define BINARY(op,intrin)                                          \
  inline v4float operator op( const v4float &a, const v4float &b )   \
  {                                                                  \
    v4float c;                                                       \
    c.v = intrin( a.v, b.v );                                        \
    return c;                                                        \
  }

  BINARY( +, vaddq_f32 )
  BINARY( -, vsubq_f32 )
  BINARY( *, vmulq_f32 )
  BINARY( /, vdivq_f32 )

  #undef BINARY

  // #define BINARY(op)                                                 \
  // inline v4float operator op( const v4float &a, const v4float &b )   \
  // {                                                               \
  //   v4float c;                                                       \
  //   ALWAYS_VECTORIZE                                                 \
  //   for( int j = 0; j < 4; j++ )                                     \
  //     c.f[j] = a.f[j] op b.f[j];                                     \
  //   return c;                                                        \
  // }

  // BINARY(+)
  // BINARY(-)
  // BINARY(*)
  // BINARY(/)

  // #undef BINARY

  // v4float logical operators

  #define LOGICAL(op,intrin)                                       \
  inline v4int operator op( const v4float &a, const v4float &b )   \
  {                                                                \
    v4int c;                                                       \
    c.v = intrin( a.v, b.v );                                      \
    return c;                                                      \
  }

  LOGICAL(  <, vcltq_f32 )
  LOGICAL(  >, vcgtq_f32 )
  LOGICAL( ==, vceqq_f32 )
  LOGICAL( <=, vcleq_f32 )
  LOGICAL( >=, vcgeq_f32 )
  // LOGICAL( !=, _mm_cmpneq_ps )

  inline v4int operator !=( const v4float &a, const v4float &b )
  {
    v4int c;

    // r.neon_u32 = vmvnq_u32(vceqq_f32(a.neon_f32, b.neon_f32));
    // return type looks wrong here. try adding uint32x4_t vi to
    // the union. may need to do a cast.

    c.vui = vmvnq_u32( vceqq_f32( a.v, b.v ) );

    return c;
  }

  inline v4int operator &&( const v4float &a, const v4float &b )
  {
    v4int c;

    float32x4_t vzero = vdupq_n_f32(0.0f);

    // __m128 vzero = _mm_setzero_ps();

    // Is there a better way to do this than the SSE way?
    c.vsi = vandq_s32( vmvnq_u32( vceqq_f32( a.v,
					     vzero ) ),
		       vmvnq_u32( vceqq_f32( b.v,
					     vzero ) ) );

    // c.v = _mm_and_ps( _mm_cmpneq_ps( a.v, vzero ),
    //                   _mm_cmpneq_ps( b.v, vzero ) );

    return c;
  }

  inline v4int operator ||( const v4float &a, const v4float &b )
  {
    v4int c;

    float32x4_t vzero = vdupq_n_f32(0.0f);

    // __m128 vzero = _mm_setzero_ps();

    // Is there a better way to do this than the SSE way?
    c.vsi = vorrq_s32( vmvnq_u32( vceqq_f32( a.v,
					     vzero ) ),
		       vmvnq_u32( vceqq_f32( b.v,
					     vzero ) ) );

    // c.v = _mm_or_ps( _mm_cmpneq_ps( a.v, vzero ),
    //                  _mm_cmpneq_ps( b.v, vzero ) );

    return c;
  }

  #undef LOGICAL

  // #define LOGICAL(op)                                              \
  // inline v4int operator op( const v4float &a, const v4float &b )   \
  // {                                                                \
  //   v4int c;                                                       \
  //   ALWAYS_VECTORIZE                                               \
  //   for( int j = 0; j < 4; j++ )                                   \
  //     c.i[j] = - ( a.f[j] op b.f[j] );                             \
  //   return c;                                                      \
  // }

  // LOGICAL(< )
  // LOGICAL(> )
  // LOGICAL(==)
  // LOGICAL(!=)
  // LOGICAL(<=)
  // LOGICAL(>=)
  // LOGICAL(&&)
  // LOGICAL(||)

  // #undef LOGICAL

  // v4float math library functions

  #define CMATH_FR1(fn)                         \
  inline v4float fn( const v4float &a )         \
  {                                             \
    v4float b;                                  \
    ALWAYS_VECTORIZE                            \
    for( int j = 0; j < 4; j++ )                \
      b.f[j] = ::fn( a.f[j] );                  \
    return b;                                   \
  }

  #define CMATH_FR2(fn)                                         \
  inline v4float fn( const v4float &a, const v4float &b )       \
  {                                                             \
    v4float c;                                                  \
    ALWAYS_VECTORIZE                                            \
    for( int j = 0; j < 4; j++ )                                \
      c.f[j] = ::fn( a.f[j], b.f[j] );                          \
    return c;                                                   \
  }

  CMATH_FR1(acos)     CMATH_FR1(asin)  CMATH_FR1(atan) CMATH_FR2(atan2)
  CMATH_FR1(ceil)     CMATH_FR1(cos)   CMATH_FR1(cosh) CMATH_FR1(exp)
  CMATH_FR1(fabs)     CMATH_FR1(floor) CMATH_FR2(fmod) CMATH_FR1(log)
  CMATH_FR1(log10)    CMATH_FR2(pow)   CMATH_FR1(sin)  CMATH_FR1(sinh)
  CMATH_FR1(sqrt)     CMATH_FR1(tan)   CMATH_FR1(tanh)

  inline v4float copysign( const v4float &a, const v4float &b )
  {
    v4float c;
    float t;

    ALWAYS_VECTORIZE
    for( int j = 0; j < 4; j++ )
    {
      t = ::fabs( a.f[j] );
      if( b.f[j] < 0 ) t = -t;
      c.f[j] = t;
    }

    return c;
  }

  #undef CMATH_FR1
  #undef CMATH_FR2

  // v4float miscellaneous functions

  inline v4float rsqrt_approx( const v4float &a )
  {
    v4float b;

    b.v = vrsqrteq_f32( a.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.f[j] = ::sqrt( 1.0f / a.f[j] );

    return b;
  }

  inline v4float rsqrt( const v4float &a )
  {
    v4float b;

    float32x4_t a_v = a.v, b_v;

    b_v = vrsqrteq_f32( a_v );

    // Note: It is quicker to just call div_ps and sqrt_ps if more
    // refinement desired!
    b.v = vaddq_f32( b_v, vmulq_f32( vdupq_n_f32( 0.5f ),
                                     vsubq_f32( b_v,
                                                vmulq_f32( a_v,
                                                           vmulq_f32( b_v,
                                                                      vmulq_f32( b_v, b_v )
                                                                    )
                                                         )
                                              )
                                   )
                    );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.f[j] = ::sqrt( 1.0f / a.f[j] );

    return b;
  }

  inline v4float rcp_approx( const v4float &a )
  {
    v4float b;

    b.v = vrecpeq_f32( a.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.f[j] = 1.0f / a.f[j];

    return b;
  }

  inline v4float rcp( const v4float &a )
  {
    v4float b;

    float32x4_t a_v = a.v, b_v;

    b_v = vrecpeq_f32( a_v );

    b.v = vsubq_f32( vaddq_f32( b_v, b_v ),
                     vmulq_f32( a_v,
                                vmulq_f32( b_v, b_v )
                              )
                   );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.f[j] = 1.0f / a.f[j];

    return b;
  }

  inline v4float fma( const v4float &a, const v4float &b, const v4float &c )
  {
    v4float d;

    d.v = vfmaq_f32( a.v, b.v, c.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   d.f[j] = a.f[j] * b.f[j] + c.f[j];

    return d;
  }

  inline v4float fms( const v4float &a, const v4float &b, const v4float &c )
  {
    v4float d;

    d.v = vfmsq_f32( a.v, b.v, c.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   d.f[j] = a.f[j] * b.f[j] - c.f[j];

    return d;
  }

  inline v4float fnms( const v4float &a, const v4float &b, const v4float &c )
  {
    v4float d;

    d.v = vsubq_f32( vdupq_n_f32( 0.0f ), vfmsq_f32( a.v, b.v, c.v ) );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   d.f[j] = c.f[j] - a.f[j] * b.f[j];

    return d;
  }

  inline v4float clear_bits( const v4int &m, const v4float &a )
  {
    v4float b;

    b.vsi = vbicq_s32( m.vsi, a.vsi );

    // b.v = _mm_andnot_ps( m.v, a.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.i[j] = ( ~m.i[j] ) & a.i[j];

    return b;
  }

  inline v4float set_bits( const v4int &m, const v4float &a )
  {
    v4float b;

    b.vsi = vorrq_s32( m.vsi, a.vsi );

    // b.v = _mm_or_ps( m.v, a.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.i[j] = m.i[j] | a.i[j];

    return b;
  }

  inline v4float toggle_bits( const v4int &m, const v4float &a )
  {
    v4float b;

    b.vsi = veorq_s32( m.vsi, a.vsi );

    // b.v = _mm_xor_ps( m.v, a.v );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   b.i[j] = m.i[j] ^ a.i[j];

    return b;
  }

  inline void increment_4x1( float * ALIGNED(16) p, const v4float &a )
  {
    vst1q_f32( p, vaddq_f32( vld1q_f32( p ), a.v ) );

    // _mm_store_ps( p, _mm_add_ps( _mm_load_ps( p ), a.v ) );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   p[j] += a.f[j];
  }

  inline void decrement_4x1( float * ALIGNED(16) p, const v4float &a )
  {
    vst1q_f32( p, vsubq_f32( vld1q_f32( p ), a.v ) );

    // _mm_store_ps( p, _mm_sub_ps( _mm_load_ps( p ), a.v ) );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   p[j] -= a.f[j];
  }

  inline void scale_4x1( float * ALIGNED(16) p, const v4float &a )
  {
    vst1q_f32( p, vmulq_f32( vld1q_f32( p ), a.v ) );

    // _mm_store_ps( p, _mm_mul_ps( _mm_load_ps( p ), a.v ) );

    // ALWAYS_VECTORIZE
    // for( int j = 0; j < 4; j++ )
    //   p[j] *= a.f[j];
  }

  inline void trilinear( v4float & wl, v4float & wh )
  {
    float x = wl.f[0], y = wl.f[1], z = wl.f[2];

    wl.f[0] = ( ( 1.0f - x ) * ( 1.0f - y ) ) * ( 1.0f - z );
    wl.f[1] = ( ( 1.0f + x ) * ( 1.0f - y ) ) * ( 1.0f - z );
    wl.f[2] = ( ( 1.0f - x ) * ( 1.0f + y ) ) * ( 1.0f - z );
    wl.f[3] = ( ( 1.0f + x ) * ( 1.0f + y ) ) * ( 1.0f - z );

    wh.f[0] = ( ( 1.0f - x ) * ( 1.0f - y ) ) * ( 1.0f + z );
    wh.f[1] = ( ( 1.0f + x ) * ( 1.0f - y ) ) * ( 1.0f + z );
    wh.f[2] = ( ( 1.0f - x ) * ( 1.0f + y ) ) * ( 1.0f + z );
    wh.f[3] = ( ( 1.0f + x ) * ( 1.0f + y ) ) * ( 1.0f + z );
  }

} // namespace v4

#endif // _v4_neon_h_
