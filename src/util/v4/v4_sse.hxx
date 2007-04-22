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
    friend inline v4 splat( const v4 &a, const int n );
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
  
  inline v4 splat( const v4 & a, const int n ) {
    __m128 a_v = a.v;
    v4 b;
    // Note: _mm_shuffle_ps has uses an imm8. Compiler can usually
    //       optimize out this switch at compile time.
    switch(n) {
    case 0: b.v = _mm_shuffle_ps( a_v, a_v, 0x00 ); break;
    case 1: b.v = _mm_shuffle_ps( a_v, a_v, 0x55 ); break;
    case 2: b.v = _mm_shuffle_ps( a_v, a_v, 0xaa ); break;
    case 3: b.v = _mm_shuffle_ps( a_v, a_v, 0xff ); break;
    default: break;
    }
    return b;
  }

  inline void swap( v4 &a, v4 &b ) { 
    __m128 t = a.v, u = b.v;
    a.v = u;
    b.v = t;
  }

  inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 ) {
    __m128 a0_v = a0.v, a1_v = a1.v, a2_v = a2.v, a3_v = a3.v, t;

    t    = a0_v;                               // t  <-  1  2  3  4
    a0_v = _mm_movelh_ps( a0_v, a1_v );        // a0 <-  1  2  5  6
    a1_v = _mm_movehl_ps( a1_v,    t );        // a1 <-  3  4  7  8
    t    = a2_v;                               // t  <-  9 10 11 12
    a2_v = _mm_movelh_ps( a2_v, a3_v );        // a2 <-  9 10 13 14
    a3_v = _mm_movehl_ps( a3_v,    t );        // a3 <- 11 12 15 16
    t    = a0_v;                               // t  <-  1  2  5  6
    a0_v = _mm_shuffle_ps( a0_v, a2_v, 0x88 ); // a0 <-  1  5  9 13
    t    = _mm_shuffle_ps(    t, a2_v, 0xdd ); // t  <-  2  6 10 13
    a2_v = a1_v;                               // a2 <-  3  4  7  8
    a1_v = t;                                  // a1 <-  1  5  9 13
    t    = a2_v;                               // t  <-  3  4  7  8
    a2_v = _mm_shuffle_ps( a2_v, a3_v, 0x88 ); // a2 <-  3  7 11 15
    a3_v = _mm_shuffle_ps(    t, a3_v, 0xdd ); // a3 <-  4  8 12 16

    a0.v = a0_v; a1.v = a1_v; a2.v = a2_v; a3.v = a3_v;
  }

  // v4 memory manipulation functions
  
  inline void load_4x1( const void * ALIGNED(16) p, v4 &a ) {
    a.v = _mm_load_ps((float * ALIGNED(16))p);
  }

  inline void store_4x1( const v4 &a, void * ALIGNED(16) p ) {
    _mm_store_ps((float * ALIGNED(16))p,a.v);
  }

  inline void stream_4x1( const v4 &a, void * ALIGNED(16) p ) {
    _mm_stream_ps((float * ALIGNED(16))p,a.v);
  }

  // FIXME: Ordering semantics
  inline void copy_4x1( void * ALIGNED(16) dst, const void * ALIGNED(16) src ) {
    _mm_store_ps( (float * ALIGNED(16))dst,
                  _mm_load_ps( (float * ALIGNED(16))src ) );
  }

  inline void swap_4x1( void * ALIGNED(16) a, void * ALIGNED(16) b ) {
    __m128 t = _mm_load_ps((float * ALIGNED(16))a);
    _mm_store_ps( (float * ALIGNED(16))a,
                  _mm_load_ps( (float * ALIGNED(16))b ) );
    _mm_store_ps( (float * ALIGNED(16))b, t );
  }

  // v4 transposed memory manipulation functions

  inline void load_4x1_tr( const void *a0, const void *a1,
                           const void *a2, const void *a3, v4 &a ) {
    a.i[0] = ((const int *)a0)[0];
    a.i[1] = ((const int *)a1)[0];
    a.i[2] = ((const int *)a2)[0];
    a.i[3] = ((const int *)a3)[0];
  }
  
  inline void load_4x2_tr( const void * ALIGNED(8) a0,
                           const void * ALIGNED(8) a1,
                           const void * ALIGNED(8) a2,
                           const void * ALIGNED(8) a3,
                           v4 &a, v4 &b ) {
    __m128 t = a.v, a_v = t, b_v;
    a_v = _mm_loadl_pi(a_v,(__m64 * ALIGNED(8))a0); // a0 b0 .. .. -> a
    a_v = _mm_loadh_pi(a_v,(__m64 * ALIGNED(8))a1); // a0 b0 a1 b1 -> a
    b_v = a_v;                                      // a0 b0 a1 b1 -> b
    t =   _mm_loadl_pi(t,  (__m64 * ALIGNED(8))a2); // a2 b2 .. .. -> t
    t =   _mm_loadh_pi(t,  (__m64 * ALIGNED(8))a3); // a2 b2 a3 b3 -> t
    a_v = _mm_shuffle_ps(a_v,t,0x88);               // a0 a1 a2 a3 -> a
    b_v = _mm_shuffle_ps(b_v,t,0xdd);               // b0 b1 b2 b3 -> b
    a.v = a_v;
    b.v = b_v;
  }
  
  inline void load_4x3_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
                           v4 &a, v4 &b, v4 &c ) {
    __m128 t = a.v, u = t, a_v = t, b_v, c_v = t;
    a_v = _mm_loadl_pi(a_v, (__m64 * ALIGNED(16))a0);    // a0 b0 .. .. -> a
    c_v = _mm_loadl_pi(c_v,((__m64 * ALIGNED(16))a0)+1); // c0 .. .. .. -> c
    a_v = _mm_loadh_pi(a_v, (__m64 * ALIGNED(16))a1);    // a0 b0 a1 b1 -> a
    c_v = _mm_loadh_pi(c_v,((__m64 * ALIGNED(16))a1)+1); // c0 .. c1 .. -> c
    b_v = a_v;                                           // a0 b0 a1 b1 -> b
    t =   _mm_loadl_pi(t,   (__m64 * ALIGNED(16))a2);    // a2 b2 .. .. -> t 
    u =   _mm_loadl_pi(u,  ((__m64 * ALIGNED(16))a2)+1); // c2 .. .. .. -> u 
    t =   _mm_loadh_pi(t,   (__m64 * ALIGNED(16))a3);    // a2 b2 a3 b3 -> t
    u =   _mm_loadh_pi(u,  ((__m64 * ALIGNED(16))a3)+1); // c2 .. c3 .. -> u
    a.v = _mm_shuffle_ps(a_v,t,0x88);                    // a0 a1 a2 a3 -> a
    b.v = _mm_shuffle_ps(b_v,t,0xdd);                    // b0 b1 b2 b3 -> b
    c.v = _mm_shuffle_ps(c_v,u,0x88);                    // c0 c1 c2 c3 -> c
  }

  inline void load_4x4_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
                           v4 &a, v4 &b, v4 &c, v4 &d ) {
    __m128 t = a.v, u = t, a_v = t, b_v, c_v = t, d_v;
    a_v = _mm_loadl_pi(a_v, (__m64 * ALIGNED(16))a0);    // a0 b0 .. .. -> a
    c_v = _mm_loadl_pi(c_v,((__m64 * ALIGNED(16))a0)+1); // c0 d0 .. .. -> c
    a_v = _mm_loadh_pi(a_v, (__m64 * ALIGNED(16))a1);    // a0 b0 a1 b1 -> a
    c_v = _mm_loadh_pi(c_v,((__m64 * ALIGNED(16))a1)+1); // c0 d0 c1 d1 -> c
    b_v = a_v;                                           // a0 b0 a1 b1 -> b
    d_v = c_v;                                           // c0 d0 c1 d1 -> d
    t =   _mm_loadl_pi(t,   (__m64 * ALIGNED(16))a2);    // a2 b2 .. .. -> t 
    u =   _mm_loadl_pi(u,  ((__m64 * ALIGNED(16))a2)+1); // c2 d2 .. .. -> u 
    t =   _mm_loadh_pi(t,   (__m64 * ALIGNED(16))a3);    // a2 b2 a3 b3 -> t
    u =   _mm_loadh_pi(u,  ((__m64 * ALIGNED(16))a3)+1); // c2 d2 c3 d3 -> u
    a.v = _mm_shuffle_ps(a_v,t,0x88);                    // a0 a1 a2 a3 -> a
    b.v = _mm_shuffle_ps(b_v,t,0xdd);                    // b0 b1 b2 b3 -> b
    c.v = _mm_shuffle_ps(c_v,u,0x88);                    // c0 c1 c2 c3 -> c
    d.v = _mm_shuffle_ps(d_v,u,0xdd);                    // d0 d1 d2 d3 -> d
  }

  inline void store_4x1_tr( const v4 &a,
                            void *a0, void *a1, void *a2, void *a3 ) {
    ((int *)a0)[0] = a.i[0];
    ((int *)a1)[0] = a.i[1];
    ((int *)a2)[0] = a.i[2];
    ((int *)a3)[0] = a.i[3];
  }

  inline void store_4x2_tr( const v4 &a, const v4 &b,
                            void * ALIGNED(8) a0, void * ALIGNED(8) a1,
                            void * ALIGNED(8) a2, void * ALIGNED(8) a3 ) {
    __m128 a_v = a.v, b_v = b.v, t;
    t = _mm_unpacklo_ps(a_v,b_v);            // a0 b0 a1 b1 -> t
    _mm_storel_pi((__m64 * ALIGNED(8))a0,t); // a0 b0       -> a0
    _mm_storeh_pi((__m64 * ALIGNED(8))a1,t); // a1 b1       -> a1
    t = _mm_unpackhi_ps(a_v,b_v);            // a2 b2 a3 b3 -> t
    _mm_storel_pi((__m64 * ALIGNED(8))a2,t); // a2 b2       -> a2
    _mm_storeh_pi((__m64 * ALIGNED(8))a3,t); // a3 b3       -> a3
  }

  inline void store_4x3_tr( const v4 &a, const v4 &b, const v4 &c,
                            void * ALIGNED(16) a0, void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2, void * ALIGNED(16) a3 ) {
    __m128 a_v = a.v, b_v = b.v, t;
    t = _mm_unpacklo_ps(a_v,b_v);             // a0 b0 a1 b1 -> t
    _mm_storel_pi((__m64 * ALIGNED(16))a0,t); // a0 b0       -> a0
    _mm_storeh_pi((__m64 * ALIGNED(16))a1,t); // a1 b1       -> a1
    t = _mm_unpackhi_ps(a_v,b_v);             // a2 b2 a3 b3 -> t
    _mm_storel_pi((__m64 * ALIGNED(16))a2,t); // a2 b2       -> a2
    _mm_storeh_pi((__m64 * ALIGNED(16))a3,t); // a3 b3       -> a3
    ((int * ALIGNED(16))a0)[2] = c.i[0];
    ((int * ALIGNED(16))a1)[2] = c.i[1];
    ((int * ALIGNED(16))a2)[2] = c.i[2];
    ((int * ALIGNED(16))a3)[2] = c.i[3];
  }
  
  inline void store_4x4_tr( const v4 &a, const v4 &b, const v4 &c, const v4 &d,
                            void * ALIGNED(16) a0, void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2, void * ALIGNED(16) a3 ) {
    __m128 a_v = a.v, b_v = b.v, c_v = c.v, d_v = d.v, t, u;
    t = _mm_unpacklo_ps(a_v,b_v);                 // a0 b0 a1 b1 -> t
    u = _mm_unpacklo_ps(c_v,d_v);                 // c0 d0 c1 d1 -> u
    _mm_storel_pi( (__m64 * ALIGNED(16))a0,   t); // a0 b0 .. .. -> a0
    _mm_storel_pi(((__m64 * ALIGNED(16))a0)+1,u); // a0 b0 c0 d0 -> a0
    _mm_storeh_pi( (__m64 * ALIGNED(16))a1,   t); // a1 b1 .. .. -> a1
    _mm_storeh_pi(((__m64 * ALIGNED(16))a1)+1,u); // a1 b1 a2 b2 -> a1
    t = _mm_unpackhi_ps(a_v,b_v);                 // a0 b0 a1 b1 -> t
    u = _mm_unpackhi_ps(c_v,d_v);                 // c0 d0 c1 d1 -> u
    _mm_storel_pi( (__m64 * ALIGNED(16))a2,   t); // a2 b2 .. .. -> a2
    _mm_storel_pi(((__m64 * ALIGNED(16))a2)+1,u); // a2 b2 c3 d3 -> a2
    _mm_storeh_pi( (__m64 * ALIGNED(16))a3,   t); // a3 b3 .. .. -> a3
    _mm_storeh_pi(((__m64 * ALIGNED(16))a3)+1,u); // a3 b3 c3 d3 -> a3
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
    
    v4int() {}                              // Default constructor
    v4int( const v4int &a ) { v = a.v; }    // Copy constructor
    v4int( const v4 &a ) { v = a.v; }       // Initialize from mixed
    v4int( const int &a ) {                 // Initialize from scalar
      union { int i; float f; } u;
      u.i = a;
      v = _mm_load_ss(&(u.f));
      v = _mm_shuffle_ps(v,v,0);
    }
    v4int( const int &i0, const int &i1,
           const int &i2, const int &i3 ) { // Initialize from scalars
      i[0] = i0; i[1] = i1; i[2] = i2; i[3] = i3;
    }
    ~v4int() {};                            // Destructor
    
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
    
    int &operator()(const int n) {
      return i[n];
    }

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
    b.i[0] = a.i[0] ? 0 : -1;
    b.i[1] = a.i[1] ? 0 : -1;
    b.i[2] = a.i[2] ? 0 : -1;
    b.i[3] = a.i[3] ? 0 : -1;
    return b;
  }

  inline v4int operator ~( const v4int & a ) {
    v4int b;
    __m128 si;
    union { int i; float f; } u;
    u.i = -1;
    si  = _mm_load_ss(&(u.f));
    b.v = _mm_xor_ps(a.v,_mm_shuffle_ps(si,si,0));
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
    c.i[0] = (a.i[0] op b.i[0]) ? -1 : 0;                      \
    c.i[1] = (a.i[1] op b.i[1]) ? -1 : 0;                      \
    c.i[2] = (a.i[2] op b.i[2]) ? -1 : 0;                      \
    c.i[3] = (a.i[3] op b.i[3]) ? -1 : 0;                      \
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
    v4float( const float &a ) {                   // Initialize from scalar
      v = _mm_load_ss((float *)&a);
      v = _mm_shuffle_ps(v,v,0);
    }
    v4float( const float &f0, const float &f1,
             const float &f2, const float &f3 ) { // Initalize from scalars
      f[0] = f0; f[1] = f1; f[2] = f2; f[3] = f3;
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

    float &operator()(const int n) { return f[n]; }

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
    __m128 t;
    float one = 1.;
    t   = _mm_load_ss(&one);
    t   = _mm_add_ps(a.v,_mm_shuffle_ps(t,t,0));
    a.v = t;
    b.v = t;
    return b;
  }

  inline v4float operator --( v4float &a ) {
    v4float b;
    __m128 t;
    float one = 1.;
    t    = _mm_load_ss(&one);
    t    = _mm_sub_ps(a.v,_mm_shuffle_ps(t,t,0));
    a.v  = t;
    b.v  = t;
    return b;
  }

  // v4float postfix increment / decrement operators

  inline v4float operator ++( v4float &a, int ) {
    v4float b;
    __m128 a_v = a.v, t;
    float one = 1.;
    t   = _mm_load_ss(&one);
    a.v = _mm_add_ps(a_v,_mm_shuffle_ps(t,t,0));
    b.v = a_v;
    return b;
  }

  inline v4float operator --( v4float &a, int ) {
    v4float b;
    __m128 a_v = a.v, t;
    float one = 1.;
    t   = _mm_load_ss(&one);
    a.v = _mm_sub_ps(a_v,_mm_shuffle_ps(t,t,0));
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
    union { int i; float f; } u;
    __m128 t;
    u.i = 1<<31;
    t   = _mm_load_ss( &u.f );
    b.v = _mm_andnot_ps( _mm_shuffle_ps( t, t, 0 ), a.v );
    return b;
  }

  inline v4float sqrt( const v4float &a ) {
    v4float b;
    b.v = _mm_sqrt_ps(a.v);
    return b;
  }

  inline v4float copysign( const v4float &a, const v4float &b ) {
    v4float c;
    union { int i; float f; } u;
    __m128 t;
    u.i = 1<<31;
    t   = _mm_load_ss( &u.f );
    t   = _mm_shuffle_ps( t, t, 0 );
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
    __m128 a_v = a.v, b_v, t;
    float half = 0.5;
    b_v = _mm_rsqrt_ps(a_v);
    t   = _mm_load_ss(&half);
    // Note: It is quicker to just call div_ps and sqrt_ps if more refinement
    // desired!
    b.v = _mm_add_ps(b_v,_mm_mul_ps(_mm_shuffle_ps(t,t,0),
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
