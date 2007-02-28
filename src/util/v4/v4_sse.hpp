/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Imported from earlier V4PIC versions
 *
 */

#ifndef _v4_sse_hpp_
#define _v4_sse_hpp_

#include <ostream>
#include <cmath>

// This requires gcc-3.3 and up
// Also, Bug 12902 has not been resolved on gcc. See README.patches for
// details.

namespace v4 {
  
# include <xmmintrin.h>

  class v4;
  class v4int;
  class v4float;
  
  /////////////////////////////////////////////////////////////////////////////
    
  ///////////////////
  // v4 base class //
  ///////////////////
  
  class v4 {
    
    /////////////////////////////
    // v4 manipulation friends //
    /////////////////////////////
      
    friend inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 );
    friend inline void cmov( const v4 &c, const v4 &a, v4 &b );
    friend inline void czero( const v4 &c, v4 &a );
    friend inline void notcmov( const v4 &c, const v4 &a, v4 &b );
    friend inline void notczero( const v4 &c, v4 &a );
  
    /////////////////////////////////
    // Memory manipulation friends //
    /////////////////////////////////
        
    friend inline void load( const void *p, v4 &a );
    friend inline void half_swizzle( const void *a0, const void *a1,
                                     const void *a2, const void *a3,
                                     v4 &a, v4 &b );
    friend inline void swizzle( const void *a0, const void *a1,
                                const void *a2, const void *a3,
                                v4 &a, v4 &b, v4 &c, v4 &d );
    
    friend inline void store( const v4 &a, void *p );
    friend inline void stream( const v4 &a, void *p );
    friend inline void half_deswizzle( const v4 &a, const v4 &b,
                                       void *a0, void *a1,
                                       void *a2, void *a3 );
    friend inline void deswizzle( const v4 &a, const v4 &b,
                                  const v4 &c, const v4 &d,
                                  void *a0, void *a1, void *a2, void *a3 );
    
  protected:
    union {
      int i[4];
      float f[4];
      __m128 v;
    };
    
  public:
    v4()            {};
    v4(const v4 &a) {
      v=a.v;
    }
    ~v4()           {};
  };
  
  ///////////////////////////////
  // v4 manipulation functions //
  ///////////////////////////////

  inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 ) {
    __m128 t = a0.v;                           // t  <-  1  2  3  4
    a0.v = _mm_movelh_ps( a0.v, a1.v );        // a0 <-  1  2  5  6
    a1.v = _mm_movehl_ps( a1.v,    t );        // a1 <-  3  4  7  8
    t    = a2.v;                               // t  <-  9 10 11 12
    a2.v = _mm_movelh_ps( a2.v, a3.v );        // a2 <-  9 10 13 14
    a3.v = _mm_movehl_ps( a3.v,    t );        // a3 <- 11 12 15 16
    t    = a0.v;                               // t  <-  1  2  5  6
    a0.v = _mm_shuffle_ps( a0.v, a2.v, 0x88 ); // a0 <-  1  5  9 13
    t    = _mm_shuffle_ps(    t, a2.v, 0xdd ); // t  <-  2  6 10 13
    a2.v = a1.v;                               // a2 <-  3  4  7  8
    a1.v = t;                                  // a1 <-  1  5  9 13
    t    = a2.v;                               // t  <-  3  4  7  8
    a2.v = _mm_shuffle_ps( a2.v, a3.v, 0x88 ); // a2 <-  3  7 11 15
    a3.v = _mm_shuffle_ps(    t, a3.v, 0xdd ); // a3 <-  4  8 12 16
  }

  inline void cmov( const v4 &c, const v4 &a, v4 &b ) {
    b.v = _mm_or_ps(_mm_andnot_ps(c.v,b.v),_mm_and_ps(c.v,a.v));
  }

  inline void czero( const v4 &c, v4 &a ) {
    a.v = _mm_andnot_ps(c.v,a.v);
  }

  inline void notcmov( const v4 &c, const v4 &a, v4 &b ) {
    b.v = _mm_or_ps(_mm_and_ps(c.v,b.v),_mm_andnot_ps(c.v,a.v));
  }

  inline void notczero( const v4 &c, v4 &a ) {
    a.v = _mm_and_ps(c.v,a.v);
  }

  ///////////////////////////////////
  // Memory manipulation functions //
  ///////////////////////////////////
  
  inline void load( const void *p, v4 &a ) {
    a.v = _mm_load_ps((float *)p);
  }
  
  inline void half_swizzle( const void *a0, const void *a1,
                            const void *a2, const void *a3, v4 &a, v4 &b ) {
    __m128 t;
    a.v = _mm_loadl_pi(a.v,(__m64 *)a0); // a0 b0 .. .. -> a
    a.v = _mm_loadh_pi(a.v,(__m64 *)a1); // a0 b0 a1 b1 -> a
    b.v = a.v;                           // a0 b0 a1 b1 -> b
    t =   _mm_loadl_pi(b.v,(__m64 *)a2); // a2 b2 .. .. -> t ... avoid warn
    t =   _mm_loadh_pi(t,  (__m64 *)a3); // a2 b2 a3 b3 -> t
    a.v = _mm_shuffle_ps(a.v,t,0x88);    // a0 a1 a2 a3 -> a
    b.v = _mm_shuffle_ps(b.v,t,0xdd);    // b0 b1 b2 b3 -> b
  }
  
  inline void swizzle( const void *a0, const void *a1,
                       const void *a2, const void *a3,
                       v4 &a, v4 &b, v4 &c, v4 &d ) {
    __m128 t, u;
    a.v = _mm_loadl_pi(a.v, (__m64 *)a0);    // a0 b0 .. .. -> a
    c.v = _mm_loadl_pi(c.v,((__m64 *)a0)+1); // c0 d0 .. .. -> c
    a.v = _mm_loadh_pi(a.v, (__m64 *)a1);    // a0 b0 a1 b1 -> a
    c.v = _mm_loadh_pi(c.v,((__m64 *)a1)+1); // c0 d0 c1 d1 -> c
    b.v = a.v;                               // a0 b0 a1 b1 -> b
    d.v = c.v;                               // c0 d0 c1 d1 -> d
    t =   _mm_loadl_pi(b.v, (__m64 *)a2);    // a2 b2 .. .. -> t ... avoid warn
    u =   _mm_loadl_pi(d.v,((__m64 *)a2)+1); // c2 d2 .. .. -> u ... avoid warn
    t =   _mm_loadh_pi(t,   (__m64 *)a3);    // a2 b2 a3 b3 -> t
    u =   _mm_loadh_pi(u,  ((__m64 *)a3)+1); // c2 d2 c3 d3 -> u
    a.v = _mm_shuffle_ps(a.v,t,0x88);        // a0 a1 a2 a3 -> a
    b.v = _mm_shuffle_ps(b.v,t,0xdd);        // b0 b1 b2 b3 -> b
    c.v = _mm_shuffle_ps(c.v,u,0x88);        // c0 c1 c2 c3 -> c
    d.v = _mm_shuffle_ps(d.v,u,0xdd);        // d0 d1 d2 d3 -> d
  }

  inline void store( const v4 &a, void *p ) {
    _mm_store_ps((float *)p,a.v);
  }

  inline void stream( const v4 &a, void *p ) {
    _mm_stream_ps((float *)p,a.v);
  }

  inline void half_deswizzle( const v4 &a, const v4 &b,
                              void *a0, void *a1, void *a2, void *a3 ) {
    __m128 t;
    t = _mm_unpacklo_ps(a.v,b.v); // a0 b0 a1 b1 -> t
    _mm_storel_pi((__m64 *)a0,t); // a0 b0       -> a0
    _mm_storeh_pi((__m64 *)a1,t); // a1 b1       -> a1
    t = _mm_unpackhi_ps(a.v,b.v); // a2 b2 a3 b3 -> t
    _mm_storel_pi((__m64 *)a2,t); // a2 b2       -> a2
    _mm_storeh_pi((__m64 *)a3,t); // a3 b3       -> a3
  }
  
  inline void deswizzle( const v4 &a, const v4 &b, const v4 &c, const v4 &d,
                         void *a0, void *a1, void *a2, void *a3 ) {
    __m128 t, u;
    t = _mm_unpacklo_ps(a.v,b.v);     // a0 b0 a1 b1 -> t
    u = _mm_unpacklo_ps(c.v,d.v);     // c0 d0 c1 d1 -> u
    _mm_storel_pi( (__m64 *)a0,   t); // a0 b0 .. .. -> a0
    _mm_storel_pi(((__m64 *)a0)+1,u); // a0 b0 c0 d0 -> a0
    _mm_storeh_pi( (__m64 *)a1,   t); // a1 b1 .. .. -> a1
    _mm_storeh_pi(((__m64 *)a1)+1,u); // a1 b1 a2 b2 -> a1
    t = _mm_unpackhi_ps(a.v,b.v);     // a0 b0 a1 b1 -> t
    u = _mm_unpackhi_ps(c.v,d.v);     // c0 d0 c1 d1 -> u
    _mm_storel_pi( (__m64 *)a2,   t); // a2 b2 .. .. -> a2
    _mm_storel_pi(((__m64 *)a2)+1,u); // a2 b2 c3 d3 -> a2
    _mm_storeh_pi( (__m64 *)a3,   t); // a3 b3 .. .. -> a3
    _mm_storeh_pi(((__m64 *)a3)+1,u); // a3 b3 c3 d3 -> a3
  }

  /////////////////////////////////////////////////////////////////////////////

  /////////////////
  // v4int class //
  /////////////////

  class v4int : public v4 {

    ////////////////////
    // Friend classes //
    ////////////////////

    friend class v4float;

    //////////////////////////////////
    // Integer manipulation friends //
    //////////////////////////////////
  
    friend inline v4int abs( const v4int &a );

    //////////////////////////////////
    // Logical manipulation friends //
    //////////////////////////////////
    
    friend inline int any( const v4int &a );
    friend inline int all( const v4int &a );
    
    //////////////////////////////////
    // ostream manipulation friends //
    //////////////////////////////////
    
    friend inline std::ostream &operator <<( std::ostream &s, const v4int &a );

  public:

    ////////////////////////////////
    // Constructors / destructors //
    ////////////////////////////////
    
    v4int() {};                           // Default constructor
    v4int( const v4int &a ) { v = a.v; }  // Copy constructor
    v4int( const int a ) {                // Initialize from scalar
      v = _mm_load_ss((float *)&a);
      v = _mm_shuffle_ps(v,v,0);
    };
    v4int( const int i0, const int i1,
           const int i2, const int i3 ) { // Initialize from scalars
      i[0] = i0;
      i[1] = i1;
      i[2] = i2;
      i[3] = i3;
    };
    ~v4int() {};                          // Destructor
    
    ////////////////////////////
    // Member access operator //
    ////////////////////////////
    
    int &operator()(const int n) {
      return i[n];
    };
    
    /////////////////////////////////////
    // Overload prefix unary operators //
    /////////////////////////////////////

#   define prefix_unary(op)                     \
    inline v4int operator op() {                \
      v4int b;                                  \
      b.i[0] = (op i[0]);                       \
      b.i[1] = (op i[1]);                       \
      b.i[2] = (op i[2]);                       \
      b.i[3] = (op i[3]);                       \
      return b;                                 \
    }
    prefix_unary(+);
    prefix_unary(-);
    inline v4int operator !() {
      v4int b;
      b.i[0] = i[0] ? 0 : -1;
      b.i[1] = i[1] ? 0 : -1;
      b.i[2] = i[2] ? 0 : -1;
      b.i[3] = i[3] ? 0 : -1;
      return b;
    }
    inline v4int operator ~() {
      v4int b;
      __m128 si;
      int i = -1;
      si = _mm_load_ss((float *)&i);
      b.v = _mm_xor_ps(v,_mm_shuffle_ps(si,si,0));
      return b;
    }
    prefix_unary(++);
    prefix_unary(--);
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector
#   undef prefix_unary

    //////////////////////////////////////
    // Overload postfix unary operators //
    //////////////////////////////////////

#   define postfix_unary(op)                    \
    inline v4int operator op(int) {             \
      v4int b;                                  \
      b.i[0] = (i[0] op);                       \
      b.i[1] = (i[1] op);                       \
      b.i[2] = (i[2] op);                       \
      b.i[3] = (i[3] op);                       \
      return b;                                 \
    }
    postfix_unary(++);
    postfix_unary(--);
#   undef postfix_unary
  
    ///////////////////////////////
    // Overload binary operators //
    ///////////////////////////////
  
#   define binary(op)				\
    inline v4int operator op(const v4int &b) {	\
      v4int c;					\
      c.i[0] = i[0] op b.i[0];			\
      c.i[1] = i[1] op b.i[1];			\
      c.i[2] = i[2] op b.i[2];			\
      c.i[3] = i[3] op b.i[3];			\
      return c;					\
    }
    binary(+);
    binary(-);
    binary(*);
    binary(/);
    binary(%);
    inline v4int operator ^(const v4int &b) {
      v4int c;
      c.v = _mm_xor_ps( v, b.v );
      return c;
    }
    inline v4int operator &(const v4int &b) {
      v4int c;
      c.v = _mm_and_ps( v, b.v );
      return c;
    }
    inline v4int operator |(const v4int &b) {
      v4int c;
      c.v = _mm_or_ps( v, b.v );
      return c;
    }
    binary(<<);
    binary(>>);
#   undef binary
#   define binary_logical(op)			\
    inline v4int operator op(const v4int &b) {	\
      v4int c;					\
      c.i[0] = (i[0] op b.i[0]) ? -1 : 0;	\
      c.i[1] = (i[1] op b.i[1]) ? -1 : 0;	\
      c.i[2] = (i[2] op b.i[2]) ? -1 : 0;	\
      c.i[3] = (i[3] op b.i[3]) ? -1 : 0;	\
      return c;					\
    }
    binary_logical(<);
    binary_logical(>);
    binary_logical(==);
    binary_logical(!=);
    binary_logical(<=);
    binary_logical(>=);
    binary_logical(&&);
    binary_logical(||);
#   undef binary_logical
  
    ///////////////////////////////////
    // Overload assignment operators //
    ///////////////////////////////////
  
#   define assign(op)			        \
    inline v4int &operator op(const v4int &b) {	\
      i[0] op b.i[0];				\
      i[1] op b.i[1];				\
      i[2] op b.i[2];				\
      i[3] op b.i[3];				\
      return *this;				\
    }
    inline v4int &operator =(const v4int &b) {
      v = b.v;
      return *this;
    }
    assign(+=);
    assign(-=);
    assign(*=);
    assign(/=);
    assign(%=);
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
    assign(<<=);
    assign(>>=);
#   undef assign

  };

  ////////////////////////////////////
  // Integer manipulation functions //
  ////////////////////////////////////

  inline v4int abs( const v4int &a ) {
    v4int b;
    b.i[0] = (a.i[0]>=0) ? a.i[0] : -a.i[0];
    b.i[1] = (a.i[1]>=0) ? a.i[1] : -a.i[1];
    b.i[2] = (a.i[2]>=0) ? a.i[2] : -a.i[2];
    b.i[3] = (a.i[3]>=0) ? a.i[3] : -a.i[3];
    return b;
  }

  ////////////////////////////////////
  // Logical manipulation functions //
  ////////////////////////////////////
  
  inline int any( const v4int &a ) {
    return a.i[0] || a.i[1] || a.i[2] || a.i[3];
  }
  
  inline int all( const v4int &a ) {
    return a.i[0] && a.i[1] && a.i[2] && a.i[3];
  }
  
  ////////////////////////////////////
  // ostream manipulation functions //
  ////////////////////////////////////
  
  inline std::ostream &operator <<( std::ostream &s, const v4int &a ) {
    s << a.i[0] << " " << a.i[1] << " " << a.i[2] << " " << a.i[3];
    return s;
  }

  /////////////////////////////////////////////////////////////////////////////

  ///////////////////
  // v4float class //
  ///////////////////

  class v4float : public v4 {

    /////////////////////////////////////////
    // Floating point manipulation friends //
    /////////////////////////////////////////

#   define cmath_fr1(fn) friend inline v4float fn( const v4float &a )
#   define cmath_fr2(fn) friend inline v4float fn( const v4float &a,  \
                                                   const v4float &b )
    cmath_fr1(acos);  cmath_fr1(asin);  cmath_fr1(atan); cmath_fr2(atan2);
    cmath_fr1(ceil);  cmath_fr1(cos);   cmath_fr1(cosh); cmath_fr1(exp);
    cmath_fr1(fabs);  cmath_fr1(floor); cmath_fr2(fmod); cmath_fr1(log);
    cmath_fr1(log10); cmath_fr2(pow);   cmath_fr1(sin);  cmath_fr1(sinh);
    cmath_fr1(sqrt);  cmath_fr1(tan);   cmath_fr1(tanh);
    // Not implemented: frexp, ldexp, modf
#   if 0 // C99 library functions
    cmath_fr1(acosh);     cmath_fr1(asinh);      cmath_fr1(atanh);
    cmath_fr1(expm1);     cmath_fr1(log1p);      cmath_fr1(logb);     
    cmath_fr1(exp2);      cmath_fr1(log2);       cmath_fr2(hypot);    
    cmath_fr1(cbrt);      cmath_fr2(copysign);   cmath_fr1(erf);
    cmath_fr1(erfc);      cmath_fr1(lgamma);     cmath_fr1(tgamma);  
    cmath_fr1(rint);      cmath_fr2(nextafter);  cmath_fr2(nexttoward); 
    cmath_fr2(remainder); cmath_fr1(nearbyint);  cmath_fr1(round);
    cmath_fr1(trunc);     cmath_fr2(fdim);       cmath_fr2(fmax);
    cmath_fr2(fmin);
    // Not implemented: nan scalbn ilogb scalbln remquo lrint lround fma scalb
#   endif
#   undef cmath_fr1
#   undef cmath_fr2

    //////////////////////////////////
    // Logical manipulation friends //
    //////////////////////////////////

    friend inline int any( const v4float &a );
    friend inline int all( const v4float &a );

    //////////////////////////////////
    // ostream manipulation friends //
    //////////////////////////////////

    friend inline std::ostream &operator <<( std::ostream &s,
                                             const v4float &a );

    //////////////////////////
    // Other useful friends //
    //////////////////////////

    friend inline v4float rsqrt_approx( const v4float &a );
    friend inline v4float rsqrt( const v4float &a );
    friend inline v4float rcp_approx( const v4float &a );
    friend inline v4float rcp( const v4float &a );
    
  public:

    //////////////////
    // Constructors //
    //////////////////
    
    v4float() {};                               // Default constructor
    v4float( const v4float &a ) {               // Copy constructor
      v = a.v;
    };
    v4float( const float a ) {                  // Initialize from scalar
      v = _mm_load_ss((float *)&a);
      v = _mm_shuffle_ps(v,v,0);
    };
    v4float( const float f0, const float f1,
             const float f2, const float f3 ) { // Initalize from scalars
      f[0] = f0; f[1] = f1; f[2] = f2; f[3] = f3;
    };
    v4float( const v4int &a ) {                 // int->float conversion
      f[0] = a.i[0]; f[1] = a.i[1]; f[2] = a.i[2]; f[3] = a.i[3];
    };
    ~v4float() {};                              // Destructor

    ////////////////////////////
    // Member access operator //
    ////////////////////////////

    float &operator()(const int n) { return f[n]; };

    /////////////////////////////////////
    // Overload prefix unary operators //
    /////////////////////////////////////

    inline v4float operator +() {
      v4float b;
      b.v = v;
      return b;
    }
    inline v4float operator -() {
      v4float b;
      b.v = _mm_sub_ps(_mm_setzero_ps(),v);
      return b;
    }
    inline v4int operator !() {
      v4int b;
      b.v = _mm_cmpeq_ps(_mm_setzero_ps(),v);
      return b;
    }
    inline v4float operator ++() {
      v4float b;
      __m128 sone;
      float one = 1.;
      sone = _mm_load_ss(&one);
      v = _mm_add_ps(v,_mm_shuffle_ps(sone,sone,0));
      b.v = v;
      return b;
    }
    inline v4float operator --() {
      v4float b;
      __m128 sone;
      float one = 1.;
      sone = _mm_load_ss(&one);
      v = _mm_sub_ps(v,_mm_shuffle_ps(sone,sone,0));
      b.v = v;
      return b;
    }
    // Note: Referencing (*) and deferencing (&) apply to the whole vector

    //////////////////////////////////////
    // Overload postfix unary operators //
    //////////////////////////////////////

    inline v4float operator ++(int) {
      v4float b;
      __m128 sone;
      float one = 1.;
      sone = _mm_load_ss(&one);
      b.v = v;
      v = _mm_add_ps(v,_mm_shuffle_ps(sone,sone,0));
      return b;
    }
    inline v4float operator --(int) {
      v4float b;
      __m128 sone;
      float one = 1.;
      sone = _mm_load_ss(&one);
      b.v = v;
      v = _mm_sub_ps(v,_mm_shuffle_ps(sone,sone,0));
      return b;
    }
    
    ///////////////////////////////
    // Overload binary operators //
    ///////////////////////////////
    
#   define binary(op,intrin,ret)		\
    inline ret operator op(const v4float &b) {	\
      ret c;					\
      c.v = intrin(v,b.v);			\
      return c;					\
    }
    binary(+,_mm_add_ps,v4float);
    binary(-,_mm_sub_ps,v4float);
    binary(*,_mm_mul_ps,v4float);
    binary(/,_mm_div_ps,v4float);
    binary(<,_mm_cmplt_ps,v4int);
    binary(>,_mm_cmpgt_ps,v4int);
    binary(==,_mm_cmpeq_ps,v4int);
    binary(!=,_mm_cmpneq_ps,v4int);
    binary(<=,_mm_cmple_ps,v4int);
    binary(>=,_mm_cmpge_ps,v4int);
    inline v4int operator &&(const v4float &b) {
      v4int c;
      __m128 vzero = _mm_setzero_ps();
      c.v = _mm_and_ps(_mm_cmpneq_ps(v,vzero),_mm_cmpneq_ps(b.v,vzero));
      return c;
    }
    inline v4int operator ||(const v4float &b) {
      v4int c;
      __m128 vzero = _mm_setzero_ps();
      c.v = _mm_or_ps(_mm_cmpneq_ps(v,vzero),_mm_cmpneq_ps(b.v,vzero));
      return c;
    }
#   undef binary

    ///////////////////////////////////
    // Overload assignment operators //
    ///////////////////////////////////
    
    inline v4float &operator =(const v4float &b) {
      v = b.v;
      return *this;
    }
#   define assign(op,intrin)				\
    inline v4float &operator op(const v4float &b) {	\
      v = intrin(v,b.v);				\
      return *this;					\
    }
    assign(+=,_mm_add_ps);
    assign(-=,_mm_sub_ps);
    assign(*=,_mm_mul_ps);
    assign(/=,_mm_div_ps);
#   undef assign

  };

  ///////////////////////////////////////////
  // Floating point manipulation functions //
  ///////////////////////////////////////////

# define cmath_fr1(fn)                          \
  inline v4float fn( const v4float &a ) {       \
    v4float b;                                  \
    b.f[0] = std::fn(a.f[0]);                   \
    b.f[1] = std::fn(a.f[1]);                   \
    b.f[2] = std::fn(a.f[2]);                   \
    b.f[3] = std::fn(a.f[3]);                   \
    return b;                                   \
  }
# define cmath_fr2(fn)                                          \
  inline v4float fn( const v4float &a, const v4float &b ) {     \
    v4float c;                                                  \
    c.f[0] = std::fn(a.f[0],b.f[0]);                            \
    c.f[1] = std::fn(a.f[1],b.f[1]);                            \
    c.f[2] = std::fn(a.f[2],b.f[2]);                            \
    c.f[3] = std::fn(a.f[3],b.f[3]);                            \
    return c;                                                   \
  }
  cmath_fr1(acos)     cmath_fr1(asin)  cmath_fr1(atan) cmath_fr2(atan2)
  cmath_fr1(ceil)     cmath_fr1(cos)   cmath_fr1(cosh) cmath_fr1(exp)
  cmath_fr1(fabs)     cmath_fr1(floor) cmath_fr2(fmod) cmath_fr1(log)
  cmath_fr1(log10)    cmath_fr2(pow)   cmath_fr1(sin)  cmath_fr1(sinh)
  /*cmath_fr1(sqrt)*/ cmath_fr1(tan)   cmath_fr1(tanh)
  inline v4float sqrt( const v4float &a ) {
    v4float b;
    b.v = _mm_sqrt_ps(a.v);
    return b;
  }
  // Not implemented: frexp, ldexp, modf
# if 0 // C99 library functions
  cmath_fr1(acosh)     cmath_fr1(asinh)      cmath_fr1(atanh)
  cmath_fr1(expm1)     cmath_fr1(log1p)      cmath_fr1(logb)
  cmath_fr1(exp2)      cmath_fr1(log2)       cmath_fr2(hypot)
  cmath_fr1(cbrt)      cmath_fr2(copysign)   cmath_fr1(erf)
  cmath_fr1(erfc)      cmath_fr1(lgamma)     cmath_fr1(tgamma)
  cmath_fr1(rint)      cmath_fr2(nextafter)  cmath_fr2(nexttoward)
  cmath_fr2(remainder) cmath_fr1(nearbyint)  cmath_fr1(round)
  cmath_fr1(trunc)     cmath_fr2(fdim)       cmath_fr2(fmax)
  cmath_fr2(fmin)
  // Not implemented: nan scalbn ilogb scalbln remquo lrint lround fma scalb
# endif
# undef cmath_fr1
# undef cmath_fr2

  ////////////////////////////////////
  // Logical manipulation functions //
  ////////////////////////////////////
  
  inline int any( const v4float &a ) {
    return a.f[0] || a.f[1] || a.f[2] || a.f[3];
  }

  inline int all( const v4float &a ) {
    return a.f[0] && a.f[1] && a.f[2] && a.f[3];
  }

  ////////////////////////////////////
  // ostream manipulation functions //
  ////////////////////////////////////

  inline std::ostream &operator <<( std::ostream &s, const v4float &a ) {
    s << a.f[0] << " " << a.f[1] << " " << a.f[2] << " " << a.f[3];
    return s;
  }

  ////////////////////////////
  // Other useful functions //
  ////////////////////////////
  
  inline v4float rsqrt_approx( const v4float &a ) {
    v4float b;
    b.v = _mm_rsqrt_ps(a.v);
    return b;
  }
  
  inline v4float rsqrt( const v4float &a ) {
    v4float b;
    __m128 t;
    float half = 0.5;
    b.v = _mm_rsqrt_ps(a.v);
    t = _mm_load_ss(&half);
    b.v = _mm_add_ps(b.v,_mm_mul_ps(_mm_shuffle_ps(t,t,0),
                                    _mm_sub_ps(b.v,_mm_mul_ps(a.v,
                                                   _mm_mul_ps(b.v,
                                                   _mm_mul_ps(b.v,b.v))))));
    return b;
  }

  inline v4float rcp_approx( const v4float &a ) {
    v4float b;
    b.v = _mm_rcp_ps(a.v);
    return b;
  }
  
  inline v4float rcp( const v4float &a ) {
    v4float b;
    b.v = _mm_rcp_ps(a.v);
    b.v = _mm_sub_ps(_mm_add_ps(b.v,b.v),_mm_mul_ps(a.v,_mm_mul_ps(b.v,b.v)));
    return b;
  }

} // namespace v4

#endif
