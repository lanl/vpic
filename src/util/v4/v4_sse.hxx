#ifndef _v4_sse_hxx_
#define _v4_sse_hxx_

#ifndef IN_v4_h
#error "Do not include v4_sse.hxx directly; use v4.h"
#endif

#define V4_ACCELERATION
#define V4_SSE_ACCELERATION

#include <xmmintrin.h>
#include <math.h>

#ifndef ALIGNED
#define ALIGNED(n)
#endif

// FIXME: IN PORTABLE, ALTIVEC, SPU
// - UPDATE V4INT, V4FLOAT

// This requires gcc-3.3 and up
// Also, Bug 12902 has not been resolved on gcc-3.x.x. See README.patches for
// details.  gcc-4.x.x does not seem to have this bug but may suffer from
// other problems (use "-fno-strict-aliasing" on these platforms)

namespace v4 {

  class v4;
  class v4int;
  class v4float;
  
  ////////////////
  // v4 base class
  
  class v4 {
    
    friend class v4int;
    friend class v4float;
      
    // v4 miscellenous friends

    friend inline int any( const v4 &a );
    friend inline int all( const v4 &a );
    friend inline v4 splat( const v4 &a, int n );
    friend inline v4 shuffle( const v4 &a,
                              int i0, int i1, int i2, int i3 );
    friend inline void swap( v4 &a, v4 &b );
    friend inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 );

    // v4int miscellaneous friends

    friend inline v4 czero(    const v4int &c, const v4 &a );
    friend inline v4 notczero( const v4int &c, const v4 &a );
    friend inline v4 merge(    const v4int &c, const v4 &a, const v4 &b );

    // v4 memory manipulation friends
        
    friend inline void load_4x1( const void * ALIGNED(16) p, v4 &a );
    friend inline void store_4x1( const v4 &a, void * ALIGNED(16) p );
    friend inline void stream_4x1( const v4 &a, void * ALIGNED(16) p );
    friend inline void clear_4x1( void * ALIGNED(16) dst );
    friend inline void copy_4x1( void * ALIGNED(16) dst,
                                 const void * ALIGNED(16) src );
    friend inline void swap_4x1( void * ALIGNED(16) a, void * ALIGNED(16) b );

    // v4 transposed memory manipulation friends

    friend inline void load_4x1_tr( const void *a0, const void *a1,
                                    const void *a2, const void *a3,
                                    v4 &a );
    friend inline void load_4x2_tr( const void * ALIGNED(8) a0,
                                    const void * ALIGNED(8) a1,
                                    const void * ALIGNED(8) a2,
                                    const void * ALIGNED(8) a3,
                                    v4 &a, v4 &b );
    friend inline void load_4x3_tr( const void * ALIGNED(16) a0,
                                    const void * ALIGNED(16) a1,
                                    const void * ALIGNED(16) a2,
                                    const void * ALIGNED(16) a3,
                                    v4 &a, v4 &b, v4 &c );
    friend inline void load_4x4_tr( const void * ALIGNED(16) a0,
                                    const void * ALIGNED(16) a1,
                                    const void * ALIGNED(16) a2,
                                    const void * ALIGNED(16) a3,
                                    v4 &a, v4 &b, v4 &c, v4 &d );
    
    friend inline void store_4x1_tr( const v4 &a,
                                     void *a0, void *a1, void *a2, void *a3 );
    friend inline void store_4x2_tr( const v4 &a, const v4 &b,
                                     void * ALIGNED(8) a0,
                                     void * ALIGNED(8) a1,
                                     void * ALIGNED(8) a2,
                                     void * ALIGNED(8) a3 );
    friend inline void store_4x3_tr( const v4 &a, const v4 &b, const v4 &c,
                                     void * ALIGNED(16) a0,
                                     void * ALIGNED(16) a1,
                                     void * ALIGNED(16) a2,
                                     void * ALIGNED(16) a3 );
    friend inline void store_4x4_tr( const v4 &a, const v4 &b,
                                     const v4 &c, const v4 &d,
                                     void * ALIGNED(16) a0,
                                     void * ALIGNED(16) a1,
                                     void * ALIGNED(16) a2,
                                     void * ALIGNED(16) a3 );

  protected:

    union {
      int i[4];
      float f[4];
      __m128 v;
    };
    
  public:

    v4() {}                    // Default constructor
    v4(const v4 &a) { v=a.v; } // Copy constructor
    ~v4() {}                   // Default destructor

  };
  
  // v4 miscellaneous functions

  inline int any( const v4 &a ) {
    return a.i[0] || a.i[1] || a.i[2] || a.i[3];
  }
  
  inline int all( const v4 &a ) {
    return a.i[0] && a.i[1] && a.i[2] && a.i[3];
  }
  
  // Note: n MUST BE AN IMMEDIATE!
  inline v4 splat( const v4 & a, int n ) {
    __m128 a_v = a.v;
    v4 b;
    b.v = _mm_shuffle_ps( a_v, a_v, n*0x55 );
    return b;
  }

  // Note: i0:3 MUST BE IMMEDIATES! */
  inline v4 shuffle( const v4 & a,
                     int i0, int i1, int i2, int i3 ) {
    __m128 a_v = a.v;
    v4 b;
    b.v = _mm_shuffle_ps( a_v, a_v, i0 + i1*4 + i2*16 + i3*64 );
    return b;
  }

  inline void swap( v4 &a, v4 &b ) { 
    __m128 a_v = a.v; a.v = b.v; b.v = a_v;
  }

  inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 ) {
    __m128 a0_v = a0.v, a1_v = a1.v, a2_v = a2.v, a3_v = a3.v, t, u;
    t    = _mm_unpackhi_ps( a0_v, a1_v );
    a0_v = _mm_unpacklo_ps( a0_v, a1_v );
    u    = _mm_unpackhi_ps( a2_v, a3_v );
    a2_v = _mm_unpacklo_ps( a2_v, a3_v );
    a1_v = _mm_movehl_ps( a2_v, a0_v );
    a0_v = _mm_movelh_ps( a0_v, a2_v );
    a2_v = _mm_movelh_ps( t, u );
    a3_v = _mm_movehl_ps( u, t );
    a0.v = a0_v; a1.v = a1_v; a2.v = a2_v; a3.v = a3_v;
  }

  // v4 memory manipulation functions
  
  inline void load_4x1( const void * ALIGNED(16) p, v4 &a ) {
    a.v = _mm_load_ps((float *)p);
  }

  inline void store_4x1( const v4 &a, void * ALIGNED(16) p ) {
    _mm_store_ps((float *)p,a.v);
  }

  inline void stream_4x1( const v4 &a, void * ALIGNED(16) p ) {
    _mm_stream_ps((float *)p,a.v);
  }

  inline void clear_4x1( void * ALIGNED(16) p ) {
    _mm_store_ps( (float *)p, _mm_setzero_ps() );
  }

  inline void copy_4x1( void * ALIGNED(16) dst,
                        const void * ALIGNED(16) src ) {
    _mm_store_ps( (float *)dst, _mm_load_ps( (const float *)src ) );
  }

  /* FIXME: MAKE ROBUST AGAINST ALIASING ISSUES */
  inline void swap_4x1( void * ALIGNED(16) a, void * ALIGNED(16) b ) {
    __m128 t = _mm_load_ps((float *)a);
    _mm_store_ps( (float *)a, _mm_load_ps( (float *)b ) );
    _mm_store_ps( (float *)b, t );
  }

  // v4 transposed memory manipulation functions

  inline void load_4x1_tr( const void *a0, const void *a1,
                           const void *a2, const void *a3, v4 &a ) {
    a.f[0] = ((const float *)a0)[0];
    a.f[1] = ((const float *)a1)[0];
    a.f[2] = ((const float *)a2)[0];
    a.f[3] = ((const float *)a3)[0];
  }

  inline void load_4x2_tr( const void * ALIGNED(8) a0,
                           const void * ALIGNED(8) a1,
                           const void * ALIGNED(8) a2,
                           const void * ALIGNED(8) a3,
                           v4 &a, v4 &b ) {
    __m128 a_v, b_v, t;
    b_v = _mm_setzero_ps();
    t   = _mm_loadh_pi( _mm_loadl_pi( b_v, (__m64 *)a0 ), (__m64 *)a1 );
    b_v = _mm_loadh_pi( _mm_loadl_pi( b_v, (__m64 *)a2 ), (__m64 *)a3 );
    a_v = _mm_shuffle_ps( t, b_v, 0x88 );
    b_v = _mm_shuffle_ps( t, b_v, 0xdd );
    a.v = a_v; b.v = b_v;
  }

  inline void load_4x3_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
                           v4 &a, v4 &b, v4 &c ) {
    __m128 a_v, b_v, c_v, t, u;
    t   = _mm_load_ps( (const float *)a0 );
    b_v = _mm_load_ps( (const float *)a1 );
    c_v = _mm_load_ps( (const float *)a2 );
    u   = _mm_load_ps( (const float *)a3 );
    a_v = _mm_unpacklo_ps( t, b_v );
    b_v = _mm_unpackhi_ps( t, b_v );
    t   = _mm_unpacklo_ps( c_v, u );
    u   = _mm_unpackhi_ps( c_v, u );
    c_v = _mm_movelh_ps( b_v, u );
    b_v = _mm_movehl_ps( t, a_v );
    a_v = _mm_movelh_ps( a_v, t );
    a.v = a_v; b.v = b_v; c.v = c_v;
  }

  inline void load_4x4_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
                           v4 &a, v4 &b, v4 &c, v4 &d ) {
    __m128 a_v, b_v, c_v, d_v, t, u;
    a_v = _mm_load_ps( (const float *)a0 );
    b_v = _mm_load_ps( (const float *)a1 );
    c_v = _mm_load_ps( (const float *)a2 );
    d_v = _mm_load_ps( (const float *)a3 );
    t   = _mm_unpackhi_ps( a_v, b_v );
    a_v = _mm_unpacklo_ps( a_v, b_v );
    u   = _mm_unpackhi_ps( c_v, d_v );
    c_v = _mm_unpacklo_ps( c_v, d_v );
    b_v = _mm_movehl_ps( c_v, a_v );
    a_v = _mm_movelh_ps( a_v, c_v );
    c_v = _mm_movelh_ps( t, u );
    d_v = _mm_movehl_ps( u, t );
    a.v = a_v; b.v = b_v; c.v = c_v; d.v = d_v;
  }

  inline void store_4x1_tr( const v4 &a,
                            void *a0, void *a1, void *a2, void *a3 ) {
    ((float *)a0)[0] = a.f[0];
    ((float *)a1)[0] = a.f[1];
    ((float *)a2)[0] = a.f[2];
    ((float *)a3)[0] = a.f[3];
  }

  inline void store_4x2_tr( const v4 &a, const v4 &b,
                            void * ALIGNED(8) a0, void * ALIGNED(8) a1,
                            void * ALIGNED(8) a2, void * ALIGNED(8) a3 ) {
    __m128 a_v = a.v, b_v = b.v, t;
    t = _mm_unpacklo_ps(a_v,b_v); // a0 b0 a1 b1 -> t
    _mm_storel_pi((__m64 *)a0,t); // a0 b0       -> a0
    _mm_storeh_pi((__m64 *)a1,t); // a1 b1       -> a1
    t = _mm_unpackhi_ps(a_v,b_v); // a2 b2 a3 b3 -> t
    _mm_storel_pi((__m64 *)a2,t); // a2 b2       -> a2
    _mm_storeh_pi((__m64 *)a3,t); // a3 b3       -> a3
  }

  inline void store_4x3_tr( const v4 &a, const v4 &b, const v4 &c,
                            void * ALIGNED(16) a0, void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2, void * ALIGNED(16) a3 ) {
    __m128 a_v = a.v, b_v = b.v, t;
    t = _mm_unpacklo_ps(a_v,b_v); // a0 b0 a1 b1 -> t
    _mm_storel_pi((__m64 *)a0,t); // a0 b0       -> a0
    _mm_storeh_pi((__m64 *)a1,t); // a1 b1       -> a1
    t = _mm_unpackhi_ps(a_v,b_v); // a2 b2 a3 b3 -> t
    _mm_storel_pi((__m64 *)a2,t); // a2 b2       -> a2
    _mm_storeh_pi((__m64 *)a3,t); // a3 b3       -> a3
    ((float *)a0)[2] = c.f[0];
    ((float *)a1)[2] = c.f[1];
    ((float *)a2)[2] = c.f[2];
    ((float *)a3)[2] = c.f[3];
  }

  /* FIXME: IS THIS FASTER THAN THE OLD WAY (HAD MORE STORE INSTR) */
  inline void store_4x4_tr( const v4 &a, const v4 &b, const v4 &c, const v4 &d,
                            void * ALIGNED(16) a0, void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2, void * ALIGNED(16) a3 ) {
    __m128 a_v = a.v, b_v = b.v, c_v = c.v, d_v = d.v, t, u;
    t   = _mm_unpackhi_ps( a_v, b_v );
    a_v = _mm_unpacklo_ps( a_v, b_v );
    u   = _mm_unpackhi_ps( c_v, d_v );
    c_v = _mm_unpacklo_ps( c_v, d_v );
    b_v = _mm_movehl_ps( c_v, a_v );
    a_v = _mm_movelh_ps( a_v, c_v );
    c_v = _mm_movelh_ps( t, u );
    d_v = _mm_movehl_ps( u, t );
    _mm_store_ps( (float *)a0, a_v );
    _mm_store_ps( (float *)a1, b_v );
    _mm_store_ps( (float *)a2, c_v );
    _mm_store_ps( (float *)a3, d_v );
  }

  //////////////
  // v4int class

  class v4int : public v4 {

    // v4int prefix unary operator friends

    friend inline v4int operator  +( const v4int & a );
    friend inline v4int operator  -( const v4int & a );
    friend inline v4int operator  ~( const v4int & a );
    friend inline v4int operator  !( const v4int & a );
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v4int prefix increment / decrement operator friends

    friend inline v4int operator ++( v4int & a );
    friend inline v4int operator --( v4int & a );

    // v4int postfix increment / decrement operator friends

    friend inline v4int operator ++( v4int & a, int );
    friend inline v4int operator --( v4int & a, int );

    // v4int binary operator friends

    friend inline v4int operator  +( const v4int &a, const v4int &b );
    friend inline v4int operator  -( const v4int &a, const v4int &b );
    friend inline v4int operator  *( const v4int &a, const v4int &b );
    friend inline v4int operator  /( const v4int &a, const v4int &b );
    friend inline v4int operator  %( const v4int &a, const v4int &b );
    friend inline v4int operator  ^( const v4int &a, const v4int &b );
    friend inline v4int operator  &( const v4int &a, const v4int &b );
    friend inline v4int operator  |( const v4int &a, const v4int &b );
    friend inline v4int operator <<( const v4int &a, const v4int &b );
    friend inline v4int operator >>( const v4int &a, const v4int &b );

    // v4int logical operator friends

    friend inline v4int operator  <( const v4int &a, const v4int &b );
    friend inline v4int operator  >( const v4int &a, const v4int &b );
    friend inline v4int operator ==( const v4int &a, const v4int &b );
    friend inline v4int operator !=( const v4int &a, const v4int &b );
    friend inline v4int operator <=( const v4int &a, const v4int &b );
    friend inline v4int operator >=( const v4int &a, const v4int &b );
    friend inline v4int operator &&( const v4int &a, const v4int &b );
    friend inline v4int operator ||( const v4int &a, const v4int &b );

    // v4int miscellaneous friends

    friend inline v4int abs( const v4int &a );
    friend inline v4    czero( const v4int &c, const v4 &a );
    friend inline v4 notczero( const v4int &c, const v4 &a );
    // FIXME: cswap, notcswap!
    friend inline v4 merge( const v4int &c, const v4 &t, const v4 &f );

    // v4float unary operator friends

    friend inline v4int operator  !( const v4float & a ); 

    // v4float logical operator friends

    friend inline v4int operator  <( const v4float &a, const v4float &b );
    friend inline v4int operator  >( const v4float &a, const v4float &b );
    friend inline v4int operator ==( const v4float &a, const v4float &b );
    friend inline v4int operator !=( const v4float &a, const v4float &b );
    friend inline v4int operator <=( const v4float &a, const v4float &b );
    friend inline v4int operator >=( const v4float &a, const v4float &b );
    friend inline v4int operator &&( const v4float &a, const v4float &b );
    friend inline v4int operator ||( const v4float &a, const v4float &b );

    // v4float miscellaneous friends

    friend inline v4float clear_bits(  const v4int &m, const v4float &a );
    friend inline v4float set_bits(    const v4int &m, const v4float &a );
    friend inline v4float toggle_bits( const v4int &m, const v4float &a );

  public:

    // v4int constructors / destructors
    
    v4int() {}                                // Default constructor
    v4int( const v4int &a ) { v = a.v; }      // Copy constructor
    v4int( const v4 &a ) { v = a.v; }         // Init from mixed
    v4int( int a ) {                          // Init from scalar
      union { int i; float f; } u;
      u.i = a;
      v = _mm_set1_ps( u.f );
    }
    v4int( int i0, int i1, int i2, int i3 ) { // Init from scalars
      union { int i; float f; } u0, u1, u2, u3;
      u0.i = i0; u1.i = i1; u2.i = i2; u3.i = i3;
      v = _mm_setr_ps( u0.f, u1.f, u2.f, u3.f );
    }
    ~v4int() {};                              // Destructor
    
    // v4int assignment operators
  
#   define ASSIGN(op)			          \
    inline v4int &operator op( const v4int &b ) { \
      i[0] op b.i[0];                             \
      i[1] op b.i[1];                             \
      i[2] op b.i[2];                             \
      i[3] op b.i[3];                             \
      return *this;                               \
    }

    inline v4int &operator =(const v4int &b) {
      v = b.v;
      return *this;
    }

    ASSIGN(+=)
    ASSIGN(-=)
    ASSIGN(*=)
    ASSIGN(/=)
    ASSIGN(%=)

    inline v4int &operator ^=(const v4int &b) {
      v = _mm_xor_ps( v, b.v );
      return *this;
    }

    inline v4int &operator &=(const v4int &b) {
      v = _mm_and_ps( v, b.v );
      return *this;
    }

    inline v4int &operator |=(const v4int &b) {
      v = _mm_or_ps( v, b.v );
      return *this;
    }

    ASSIGN(<<=)
    ASSIGN(>>=)

#   undef ASSIGN

    // v4int member access operator
    
    inline int &operator []( int n ) { return i[n]; }
    inline int  operator ()( int n ) { return i[n]; }

  };

  // v4int prefix unary operators

# define PREFIX_UNARY(op)                       \
  inline v4int operator op( const v4int & a ) { \
    v4int b;                                    \
    b.i[0] = (op a.i[0]);                       \
    b.i[1] = (op a.i[1]);                       \
    b.i[2] = (op a.i[2]);                       \
    b.i[3] = (op a.i[3]);                       \
    return b;                                   \
  }

  inline v4int operator +( const v4int & a ) {
    v4int b;
    b.v = a.v;
    return b;
  }

  PREFIX_UNARY(-)

  inline v4int operator !( const v4int & a ) {
    v4int b;
    b.i[0] = -(!a.i[0]);
    b.i[1] = -(!a.i[1]);
    b.i[2] = -(!a.i[2]);
    b.i[3] = -(!a.i[3]);
    return b;
  }

  inline v4int operator ~( const v4int & a ) {
    v4int b;
    union { int i; float f; } u;
    u.i = -1;
    b.v = _mm_xor_ps( a.v, _mm_set1_ps( u.f ) );
    return b;
  }
  
# undef PREFIX_UNARY

  // v4int prefix increment / decrement

# define PREFIX_INCDEC(op)                      \
  inline v4int operator op( v4int & a ) {       \
    v4int b;                                    \
    b.i[0] = (op a.i[0]);                       \
    b.i[1] = (op a.i[1]);                       \
    b.i[2] = (op a.i[2]);                       \
    b.i[3] = (op a.i[3]);                       \
    return b;                                   \
  }

  PREFIX_INCDEC(++)
  PREFIX_INCDEC(--)

# undef PREFIX_INCDEC

  // v4int postfix increment / decrement

# define POSTFIX_INCDEC(op)                    \
  inline v4int operator op( v4int & a, int ) { \
    v4int b;                                   \
    b.i[0] = (a.i[0] op);                      \
    b.i[1] = (a.i[1] op);                      \
    b.i[2] = (a.i[2] op);                      \
    b.i[3] = (a.i[3] op);                      \
    return b;                                  \
  }

  POSTFIX_INCDEC(++)
  POSTFIX_INCDEC(--)

# undef POSTFIX_INCDEC

  // v4int binary operators
  
# define BINARY(op)                                             \
  inline v4int operator op( const v4int &a, const v4int &b ) {	\
    v4int c;                                                    \
    c.i[0] = a.i[0] op b.i[0];                                  \
    c.i[1] = a.i[1] op b.i[1];                                  \
    c.i[2] = a.i[2] op b.i[2];                                  \
    c.i[3] = a.i[3] op b.i[3];                                  \
    return c;                                                   \
  }

  BINARY(+)
  BINARY(-)
  BINARY(*)
  BINARY(/)
  BINARY(%)

  inline v4int operator ^( const v4int &a, const v4int &b ) {
    v4int c;
    c.v = _mm_xor_ps( a.v, b.v );
    return c;
  }

  inline v4int operator &( const v4int &a, const v4int &b ) {
    v4int c;
    c.v = _mm_and_ps( a.v, b.v );
    return c;
  }

  inline v4int operator |( const v4int &a, const v4int &b ) {
    v4int c;
    c.v = _mm_or_ps( a.v, b.v );
    return c;
  }

  BINARY(<<)
  BINARY(>>)

# undef BINARY

  // v4int logical operators

# define LOGICAL(op)                                           \
  inline v4int operator op( const v4int &a, const v4int &b ) { \
    v4int c;                                                   \
    c.i[0] = -(a.i[0] op b.i[0]);                              \
    c.i[1] = -(a.i[1] op b.i[1]);                              \
    c.i[2] = -(a.i[2] op b.i[2]);                              \
    c.i[3] = -(a.i[3] op b.i[3]);                              \
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
  
# undef LOGICAL

  // v4int miscellaneous functions

  inline v4int abs( const v4int &a ) {
    v4int b;
    b.i[0] = (a.i[0]>=0) ? a.i[0] : -a.i[0];
    b.i[1] = (a.i[1]>=0) ? a.i[1] : -a.i[1];
    b.i[2] = (a.i[2]>=0) ? a.i[2] : -a.i[2];
    b.i[3] = (a.i[3]>=0) ? a.i[3] : -a.i[3];
    return b;
  }

  inline v4 czero( const v4int &c, const v4 &a ) {
    v4 b;
    b.v = _mm_andnot_ps(c.v,a.v);
    return b;
  }

  inline v4 notczero( const v4int &c, const v4 &a ) {
    v4 b;
    b.v = _mm_and_ps(c.v,a.v);
    return b;
  }
  
  inline v4 merge( const v4int &c, const v4 &t, const v4 &f ) {
    __m128 c_v = c.v;
    v4 tf;
    tf.v = _mm_or_ps(_mm_andnot_ps(c_v,f.v),_mm_and_ps(c_v,t.v));
    return tf;
  }

  ////////////////
  // v4float class

  class v4float : public v4 {

    // v4float prefix unary operator friends

    friend inline v4float operator  +( const v4float &a );
    friend inline v4float operator  -( const v4float &a );
    friend inline v4float operator  ~( const v4float &a );
    friend inline v4int   operator  !( const v4float &a );
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v4float prefix increment / decrement operator friends

    friend inline v4float operator ++( v4float &a );
    friend inline v4float operator --( v4float &a );

    // v4float postfix increment / decrement operator friends

    friend inline v4float operator ++( v4float &a, int );
    friend inline v4float operator --( v4float &a, int );

    // v4float binary operator friends

    friend inline v4float operator  +( const v4float &a, const v4float &b );
    friend inline v4float operator  -( const v4float &a, const v4float &b );
    friend inline v4float operator  *( const v4float &a, const v4float &b );
    friend inline v4float operator  /( const v4float &a, const v4float &b );

    // v4float logical operator friends

    friend inline v4int operator  <( const v4float &a, const v4float &b );
    friend inline v4int operator  >( const v4float &a, const v4float &b );
    friend inline v4int operator ==( const v4float &a, const v4float &b );
    friend inline v4int operator !=( const v4float &a, const v4float &b );
    friend inline v4int operator <=( const v4float &a, const v4float &b );
    friend inline v4int operator >=( const v4float &a, const v4float &b );
    friend inline v4int operator &&( const v4float &a, const v4float &b );
    friend inline v4int operator ||( const v4float &a, const v4float &b );

    // v4float math library friends

#   define CMATH_FR1(fn) friend inline v4float fn( const v4float &a )
#   define CMATH_FR2(fn) friend inline v4float fn( const v4float &a,  \
                                                   const v4float &b )

    CMATH_FR1(acos);  CMATH_FR1(asin);  CMATH_FR1(atan); CMATH_FR2(atan2);
    CMATH_FR1(ceil);  CMATH_FR1(cos);   CMATH_FR1(cosh); CMATH_FR1(exp);
    CMATH_FR1(fabs);  CMATH_FR1(floor); CMATH_FR2(fmod); CMATH_FR1(log);
    CMATH_FR1(log10); CMATH_FR2(pow);   CMATH_FR1(sin);  CMATH_FR1(sinh);
    CMATH_FR1(sqrt);  CMATH_FR1(tan);   CMATH_FR1(tanh);

    CMATH_FR2(copysign);

#   undef CMATH_FR1
#   undef CMATH_FR2

    // v4float miscellaneous friends

    friend inline v4float rsqrt_approx( const v4float &a );
    friend inline v4float rsqrt( const v4float &a );
    friend inline v4float rcp_approx( const v4float &a );
    friend inline v4float rcp( const v4float &a );
    friend inline v4float fma(  const v4float &a, const v4float &b, const v4float &c );
    friend inline v4float fms(  const v4float &a, const v4float &b, const v4float &c );
    friend inline v4float fnms( const v4float &a, const v4float &b, const v4float &c );
    friend inline v4float clear_bits(  const v4int &m, const v4float &a );
    friend inline v4float set_bits(    const v4int &m, const v4float &a );
    friend inline v4float toggle_bits( const v4int &m, const v4float &a );
    friend inline void increment_4x1( float * ALIGNED(16) p, const v4float &a );
    friend inline void decrement_4x1( float * ALIGNED(16) p, const v4float &a );
    friend inline void scale_4x1(     float * ALIGNED(16) p, const v4float &a );
    // FIXME: crack
    
  public:

    // v4float constructors / destructors
    
    v4float() {}                                  // Default constructor
    v4float( const v4float &a ) { v = a.v; }      // Copy constructor
    v4float( const v4 &a ) { v = a.v; }           // Initialize from mixed
    v4float( float a ) {                          // Initialize from scalar
      v = _mm_set1_ps( a );
    }
    v4float( float f0, float f1,
             float f2, float f3 ) {               // Initalize from scalars
      v = _mm_setr_ps( f0, f1, f2, f3 );
    }
    ~v4float() {}                                 // Destructor

    // v4float assignment operators

#   define ASSIGN(op,intrin)				\
    inline v4float &operator op(const v4float &b) {	\
      v = intrin(v,b.v);				\
      return *this;					\
    }

    inline v4float &operator =(const v4float &b) {
      v = b.v;
      return *this;
    }

    ASSIGN(+=,_mm_add_ps)
    ASSIGN(-=,_mm_sub_ps)
    ASSIGN(*=,_mm_mul_ps)
    ASSIGN(/=,_mm_div_ps)

#   undef ASSIGN

    // v4float member access operator

    inline float &operator []( int n ) { return f[n]; }
    inline float  operator ()( int n ) { return f[n]; }

  };

  // v4float prefix unary operators

  inline v4float operator +( const v4float &a ) {
    v4float b;
    b.v = a.v;
    return b;
  }

  inline v4float operator -( const v4float &a ) {
    v4float b;
    b.v = _mm_sub_ps(_mm_setzero_ps(),a.v);
    return b;
  }

  inline v4int operator !( const v4float &a ) {
    v4int b;
    b.v = _mm_cmpeq_ps(_mm_setzero_ps(),a.v);
    return b;
  }

  // v4float prefix increment / decrement operators

  inline v4float operator ++( v4float &a ) {
    v4float b;
    __m128 t = _mm_add_ps( a.v, _mm_set1_ps( 1 ) );
    a.v = t;
    b.v = t;
    return b;
  }

  inline v4float operator --( v4float &a ) {
    v4float b;
    __m128 t = _mm_sub_ps( a.v, _mm_set1_ps( 1 ) );
    a.v = t;
    b.v = t;
    return b;
  }

  // v4float postfix increment / decrement operators

  inline v4float operator ++( v4float &a, int ) {
    v4float b;
    __m128 a_v = a.v;
    a.v = _mm_add_ps( a_v, _mm_set1_ps( 1 ) );
    b.v = a_v;
    return b;
  }

  inline v4float operator --( v4float &a, int ) {
    v4float b;
    __m128 a_v = a.v;
    a.v = _mm_sub_ps(a_v, _mm_set1_ps( 1 ) );
    b.v = a_v;
    return b;
  }

  // v4float binary operators
    
# define BINARY(op,intrin)                                           \
  inline v4float operator op( const v4float &a, const v4float &b ) { \
    v4float c;                                                       \
    c.v = intrin(a.v,b.v);                                           \
    return c;                                                        \
  }

  BINARY(+,_mm_add_ps)
  BINARY(-,_mm_sub_ps)
  BINARY(*,_mm_mul_ps)
  BINARY(/,_mm_div_ps)

# undef BINARY

  // v4float logical operators

# define LOGICAL(op,intrin)                                        \
  inline v4int operator op( const v4float &a, const v4float &b ) { \
    v4int c;                                                       \
    c.v = intrin(a.v,b.v);                                         \
    return c;                                                      \
  }

  LOGICAL(<, _mm_cmplt_ps )
  LOGICAL(>, _mm_cmpgt_ps )
  LOGICAL(==,_mm_cmpeq_ps )
  LOGICAL(!=,_mm_cmpneq_ps)
  LOGICAL(<=,_mm_cmple_ps )
  LOGICAL(>=,_mm_cmpge_ps )

  inline v4int operator &&( const v4float &a, const v4float &b ) {
    v4int c;
    __m128 vzero = _mm_setzero_ps();
    c.v = _mm_and_ps(_mm_cmpneq_ps(a.v,vzero),_mm_cmpneq_ps(b.v,vzero));
    return c;
  }

  inline v4int operator ||( const v4float &a, const v4float &b ) {
    v4int c;
    __m128 vzero = _mm_setzero_ps();
    c.v = _mm_or_ps(_mm_cmpneq_ps(a.v,vzero),_mm_cmpneq_ps(b.v,vzero));
    return c;
  }

# undef LOGICAL

  // v4float math library functions

# define CMATH_FR1(fn)                          \
  inline v4float fn( const v4float &a ) {       \
    v4float b;                                  \
    b.f[0] = ::fn(a.f[0]);                      \
    b.f[1] = ::fn(a.f[1]);                      \
    b.f[2] = ::fn(a.f[2]);                      \
    b.f[3] = ::fn(a.f[3]);                      \
    return b;                                   \
  }

# define CMATH_FR2(fn)                                          \
  inline v4float fn( const v4float &a, const v4float &b ) {     \
    v4float c;                                                  \
    c.f[0] = ::fn(a.f[0],b.f[0]);                               \
    c.f[1] = ::fn(a.f[1],b.f[1]);                               \
    c.f[2] = ::fn(a.f[2],b.f[2]);                               \
    c.f[3] = ::fn(a.f[3],b.f[3]);                               \
    return c;                                                   \
  }

  CMATH_FR1(acos)     CMATH_FR1(asin)  CMATH_FR1(atan) CMATH_FR2(atan2)
  CMATH_FR1(ceil)     CMATH_FR1(cos)   CMATH_FR1(cosh) CMATH_FR1(exp)
  /*CMATH_FR1(fabs)*/ CMATH_FR1(floor) CMATH_FR2(fmod) CMATH_FR1(log)
  CMATH_FR1(log10)    CMATH_FR2(pow)   CMATH_FR1(sin)  CMATH_FR1(sinh)
  /*CMATH_FR1(sqrt)*/ CMATH_FR1(tan)   CMATH_FR1(tanh)

  inline v4float fabs( const v4float &a ) {
    v4float b;
    b.v = _mm_andnot_ps( _mm_set1_ps( -0.f ), a.v );
    return b;
  }

  inline v4float sqrt( const v4float &a ) {
    v4float b;
    b.v = _mm_sqrt_ps(a.v);
    return b;
  }

  inline v4float copysign( const v4float &a, const v4float &b ) {
    v4float c;
    __m128 t = _mm_set1_ps( -0.f );
    c.v = _mm_or_ps( _mm_and_ps( t, b.v ), _mm_andnot_ps( t, a.v ) );
    return c;
  }

# undef CMATH_FR1
# undef CMATH_FR2

  // v4float miscelleanous functions
  
  inline v4float rsqrt_approx( const v4float &a ) {
    v4float b;
    b.v = _mm_rsqrt_ps(a.v);
    return b;
  }
  
  inline v4float rsqrt( const v4float &a ) {
    v4float b;
    __m128 a_v = a.v, b_v;
    b_v = _mm_rsqrt_ps(a_v);
    // Note: It is quicker to just call div_ps and sqrt_ps if more
    // refinement desired!
    b.v = _mm_add_ps(b_v,_mm_mul_ps(_mm_set1_ps(0.5f),
                                    _mm_sub_ps(b_v,_mm_mul_ps(a_v,
                                                   _mm_mul_ps(b_v,
                                                   _mm_mul_ps(b_v,b_v))))));
    return b;
  }

  inline v4float rcp_approx( const v4float &a ) {
    v4float b;
    b.v = _mm_rcp_ps(a.v);
    return b;
  }
  
  inline v4float rcp( const v4float &a ) {
    v4float b;
    __m128 a_v = a.v, b_v;
    b_v = _mm_rcp_ps(a_v);
    b.v = _mm_sub_ps(_mm_add_ps(b_v,b_v),_mm_mul_ps(a_v,_mm_mul_ps(b_v,b_v)));
    return b;
  }

  inline v4float fma(  const v4float &a, const v4float &b, const v4float &c ) {
    v4float d;
    d.v = _mm_add_ps( _mm_mul_ps( a.v, b.v ), c.v );
    return d;
  }

  inline v4float fms(  const v4float &a, const v4float &b, const v4float &c ) {
    v4float d;
    d.v = _mm_sub_ps( _mm_mul_ps( a.v, b.v ), c.v );
    return d;
  }

  inline v4float fnms( const v4float &a, const v4float &b, const v4float &c ) {
    v4float d;
    d.v = _mm_sub_ps( c.v, _mm_mul_ps( a.v, b.v ) );
    return d;
  }

  inline v4float clear_bits( const v4int &m, const v4float &a ) {
    v4float b;
    b.v = _mm_andnot_ps( m.v, a.v );
    return b;
  }

  inline v4float set_bits( const v4int &m, const v4float &a ) {
    v4float b;
    b.v = _mm_or_ps( m.v, a.v );
    return b;
  }

  inline v4float toggle_bits( const v4int &m, const v4float &a ) {
    v4float b;
    b.v = _mm_xor_ps( m.v, a.v );
    return b;
  }

  inline void increment_4x1( float * ALIGNED(16) p, const v4float &a ) {
    _mm_store_ps( p, _mm_add_ps( _mm_load_ps( p ), a.v ) );
  }

  inline void decrement_4x1( float * ALIGNED(16) p, const v4float &a ) {
    _mm_store_ps( p, _mm_sub_ps( _mm_load_ps( p ), a.v ) );
  }

  inline void scale_4x1( float * ALIGNED(16) p, const v4float &a ) {
    _mm_store_ps( p, _mm_mul_ps( _mm_load_ps( p ), a.v ) );
  }

} // namespace v4

#endif // _v4_sse_hxx_
