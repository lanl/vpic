#ifndef _v4_spu_hxx_
#define _v4_spu_hxx_

#ifndef IN_v4_h
#error "Do not include v4_spu.hxx directly; use v4.h"
#endif

#define V4_ACCELERATION
#define V4_SPU_ACCELERATION

#include <spu_intrinsics.h> // Requires -Wno-long-long to include
#include <stdint.h>
#include <math.h>

#ifndef ALIGNED
#define ALIGNED(n)
#endif

// FIXME: CHECK SPU_CMPEQ WITH VEC_FLOAT4 ZEROS!

namespace v4 {

  const vec_uchar16 _packe   = {  0, 1, 2, 3,    8, 9,10,11,
                                 16,17,18,19,   24,25,26,27 };
  const vec_uchar16 _packo   = {  4, 5, 6, 7,   12,13,14,15,
                                 20,21,22,23,   28,29,30,31 };
  const vec_uchar16 _unpackl = {  0, 1, 2, 3,   16,17,18,19,
                                  4, 5, 6, 7,   20,21,22,23 };
  const vec_uchar16 _unpackh = {  8, 9,10,11,   24,25,26,27,
                                 12,13,14,15,   28,29,30,31 };


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
    // Note: Half aligned values are permissible in the 4x2_tr variants!

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
      vec_float4 v;
    };
    
  public:

    v4() {}                      // Default constructor
    v4(const v4 &a) { v = a.v; } // Copy constructor
    ~v4() {}                     // Default destructor

  };
  
  // v4 miscellaneous functions

  // FIXME: Castless variant?
  inline int any( const v4 &a ) {
    return spu_extract( spu_gather( spu_cmpeq( (vec_int4)a.v,
                                               spu_splats( 0 ) ) ),
                        0 )==15 ? 0 : -1; 
  }
  
  // FIXME: Castless variant?
  inline int all( const v4 &a ) {
    return spu_extract( spu_gather( spu_cmpeq( (vec_int4)a.v,
                                               spu_splats( 0 ) ) ),
                        0 )==0  ? -1 : 0;
  }
  
  inline v4 splat( const v4 & a, const int n ) {
    v4 b;
    b.v = spu_splats( spu_extract( a.v, n ) );
    return b;
  }

  inline void swap( v4 &a, v4 &b ) { 
    vec_float4 t = a.v; a.v = b.v; b.v = t;
  }

  inline void transpose( v4 &a, v4 &b, v4 &c, v4 &d ) {
    vec_uchar16 unpackl = _unpackl;
    vec_uchar16 unpackh = _unpackh;
    vec_float4 a0 = a.v;                            // a0 =  0  1  2  3
    vec_float4 b0 = b.v;                            // b0 =  4  5  6  7
    vec_float4 c1 = c.v;                            // c1 =  8  9 10 11
    vec_float4 d1 = d.v;                            // d1 = 12 13 14 15
  
    // Step 1: Interleave top and bottom half

    vec_float4 a1 = spu_shuffle( a0, c1, unpackl ); // a1 =  0  8  1  9
    vec_float4 b1 = spu_shuffle( b0, d1, unpackl ); // b1 =  4 12  5 13
    c1            = spu_shuffle( a0, c1, unpackh ); // c1 =  2 10  3 11
    d1            = spu_shuffle( b0, d1, unpackh ); // d1 =  6 14  7 15

    // Step 2: Interleave even and odd rows

    a.v           = spu_shuffle( a1, b1, unpackl ); // a  =  0  4  8 12
    b.v           = spu_shuffle( a1, b1, unpackh ); // b  =  1  5  9 13
    c.v           = spu_shuffle( c1, d1, unpackl ); // c  =  2  6 10 14
    d.v           = spu_shuffle( c1, d1, unpackh ); // d  =  3  7 11 15
  }

  // v4 memory manipulation functions
  
  inline void load_4x1( const void * ALIGNED(16) p, v4 &a ) {
    a.v = *((const vec_float4 * ALIGNED(16))p);
  }

  inline void store_4x1( const v4 &a, void * ALIGNED(16) p ) {
    *((vec_float4 * ALIGNED(16))p) = a.v;
  }

  inline void stream_4x1( const v4 &a, void * ALIGNED(16) p ) {
    *((vec_float4 * ALIGNED(16))p) = a.v;
  }

  // FIXME: Ordering semantics
  inline void copy_4x1(       void * ALIGNED(16) dst,
                        const void * ALIGNED(16) src ) {
    *((vec_float4 * ALIGNED(16))dst) = *((const vec_float4 * ALIGNED(16))src);
  }

  inline void swap_4x1( void * ALIGNED(16) a, void * ALIGNED(16) b ) {
    vec_float4 t                   = *((vec_float4 * ALIGNED(16))a);
    *((vec_float4 * ALIGNED(16))a) = *((vec_float4 * ALIGNED(16))b);
    *((vec_float4 * ALIGNED(16))b) = t;
  }

  // v4 transposed memory manipulation functions

  inline void load_4x1_tr( const void *a0, const void *a1,
                           const void *a2, const void *a3, v4 &a ) {
    a.f[0] = ((const float *)a0)[0];
    a.f[1] = ((const float *)a1)[0];
    a.f[2] = ((const float *)a2)[0];
    a.f[3] = ((const float *)a3)[0];
  }

  // Note: load_4x4_tr with last two shuffles discared may be
  // preferable for 128-byte aligned results.  (Compiler optimization
  // may actually do that for free!)

  inline void load_4x2_tr( const void * ALIGNED(8) pa,
                           const void * ALIGNED(8) pb,
                           const void * ALIGNED(8) pc,
                           const void * ALIGNED(8) pd,
                           v4 &a, v4 &b ) {
    // FIXME: a castless variant?
    vec_llong2 a_v = { *(const int64_t * ALIGNED(8))pa,
                       *(const int64_t * ALIGNED(8))pb }; // 0 4 1 5
    vec_llong2 b_v = { *(const int64_t * ALIGNED(8))pc,
                       *(const int64_t * ALIGNED(8))pd }; // 2 6 3 7
    a.v = (vec_float4)spu_shuffle( a_v, b_v, _packe ); // 0 1 2 3
    b.v = (vec_float4)spu_shuffle( a_v, b_v, _packo ); // 4 5 6 7
  }
  
  inline void load_4x3_tr( const void * ALIGNED(16) pa,
                           const void * ALIGNED(16) pb,
                           const void * ALIGNED(16) pc,
                           const void * ALIGNED(16) pd,
                           v4 &a, v4 &b, v4 &c ) {
    vec_uchar16 unpackl = _unpackl;
    vec_uchar16 unpackh = _unpackh;
    vec_float4 a0 = *((const vec_float4 * ALIGNED(16))pa); // a0 =  0  1  2 x
    vec_float4 b0 = *((const vec_float4 * ALIGNED(16))pb); // b0 =  4  5  6 x
    vec_float4 c1 = *((const vec_float4 * ALIGNED(16))pc); // c1 =  8  9 10 x
    vec_float4 d1 = *((const vec_float4 * ALIGNED(16))pd); // d1 = 12 13 14 x
  
    // Step 1: Interleave top and bottom half

    vec_float4 a1 = spu_shuffle( a0, c1, unpackl ); // a1 =  0  8  1  9
    vec_float4 b1 = spu_shuffle( b0, d1, unpackl ); // b1 =  4 12  5 13
    c1            = spu_shuffle( a0, c1, unpackh ); // c1 =  2 10  x  x
    d1            = spu_shuffle( b0, d1, unpackh ); // d1 =  6 14  x  x

    // Step 2: Interleave even and odd rows

    a.v           = spu_shuffle( a1, b1, unpackl ); // a  =  0  4  8 12
    b.v           = spu_shuffle( a1, b1, unpackh ); // b  =  1  5  9 13
    c.v           = spu_shuffle( c1, d1, unpackl ); // c  =  2  6 10 14
  }

  inline void load_4x4_tr( const void * ALIGNED(16) pa,
                           const void * ALIGNED(16) pb,
                           const void * ALIGNED(16) pc,
                           const void * ALIGNED(16) pd,
                           v4 &a, v4 &b, v4 &c, v4 &d ) {
    vec_uchar16 unpackl = _unpackl;
    vec_uchar16 unpackh = _unpackh;
    vec_float4 a0 = *((const vec_float4 * ALIGNED(16))pa); // a0 =  0  1  2  3
    vec_float4 b0 = *((const vec_float4 * ALIGNED(16))pb); // b0 =  4  5  6  7
    vec_float4 c1 = *((const vec_float4 * ALIGNED(16))pc); // c1 =  8  9 10 11
    vec_float4 d1 = *((const vec_float4 * ALIGNED(16))pd); // d1 = 12 13 14 15
  
    // Step 1: Interleave top and bottom half

    vec_float4 a1 = spu_shuffle( a0, c1, unpackl ); // a1 =  0  8  1  9
    vec_float4 b1 = spu_shuffle( b0, d1, unpackl ); // b1 =  4 12  5 13
    c1            = spu_shuffle( a0, c1, unpackh ); // c1 =  2 10  3 11
    d1            = spu_shuffle( b0, d1, unpackh ); // d1 =  6 14  7 15

    // Step 2: Interleave even and odd rows

    a.v           = spu_shuffle( a1, b1, unpackl ); // a  =  0  4  8 12
    b.v           = spu_shuffle( a1, b1, unpackh ); // b  =  1  5  9 13
    c.v           = spu_shuffle( c1, d1, unpackl ); // c  =  2  6 10 14
    d.v           = spu_shuffle( c1, d1, unpackh ); // c  =  3  7 11 15

  }

  inline void store_4x1_tr( const v4 &a,
                            void *a0, void *a1, void *a2, void *a3 ) {
    ((float *)a0)[0] = a.f[0];
    ((float *)a1)[0] = a.f[1];
    ((float *)a2)[0] = a.f[2];
    ((float *)a3)[0] = a.f[3];
  }

  inline void store_4x2_tr( const v4 &a, const v4 &b,
                            void * ALIGNED(8) pa,
                            void * ALIGNED(8) pb,
                            void * ALIGNED(8) pc,
                            void * ALIGNED(8) pd ) {
    // FIXME: A castless variant??
    vec_llong2 t  = (vec_llong2)a.v;                // t  =  0  1  2  3
    vec_llong2 b1 = (vec_llong2)b.v;                // b1 =  4  5  6  7

    vec_llong2 a1 = spu_shuffle( t, b1, _unpackl ); // a1 =  0  4  1  5
    b1            = spu_shuffle( t, b1, _unpackh ); // b1 =  2  6  3  7
    ((int64_t * ALIGNED(8))pa)[0] = spu_extract( a1, 0 );
    ((int64_t * ALIGNED(8))pb)[0] = spu_extract( a1, 1 );
    ((int64_t * ALIGNED(8))pc)[0] = spu_extract( b1, 0 );
    ((int64_t * ALIGNED(8))pd)[0] = spu_extract( b1, 1 );
  }

  inline void store_4x3_tr( const v4 &a, const v4 &b, const v4 &c,
                            void * ALIGNED(16) pa, void * ALIGNED(16) pb,
                            void * ALIGNED(16) pc, void * ALIGNED(16) pd ) {
    // FIXME: A castless variant??
    vec_llong2 t  = (vec_llong2)a.v;                // t  =  0  1  2  3
    vec_llong2 b1 = (vec_llong2)b.v;                // b1 =  4  5  6  7
    vec_float4 c1 = c.v;                            // c1 =  8  9 10 12

    vec_llong2 a1 = spu_shuffle( t, b1, _unpackl ); // a1 =  0  4  1  5
    b1            = spu_shuffle( t, b1, _unpackh ); // b1 =  2  6  3  7

    ((int64_t * ALIGNED(8) )pa)[0] = spu_extract( a1, 0 );
    ((float   * ALIGNED(16))pa)[2] = spu_extract( c1, 0 );
    ((int64_t * ALIGNED(8) )pb)[0] = spu_extract( a1, 1 );
    ((float   * ALIGNED(16))pb)[2] = spu_extract( c1, 1 );
    ((int64_t * ALIGNED(8) )pc)[0] = spu_extract( b1, 0 );
    ((float   * ALIGNED(16))pc)[2] = spu_extract( c1, 2 );
    ((int64_t * ALIGNED(8) )pd)[0] = spu_extract( b1, 1 );
    ((float   * ALIGNED(16))pd)[2] = spu_extract( c1, 3 );
  }
  
  inline void store_4x4_tr( const v4 &a, const v4 &b, const v4 &c, const v4 &d,
                            void * ALIGNED(16) pa, void * ALIGNED(16) pb,
                            void * ALIGNED(16) pc, void * ALIGNED(16) pd ) {
    vec_uchar16 unpackl = _unpackl;
    vec_uchar16 unpackh = _unpackh;
    vec_float4 a0 = a.v;                            // a0 =  0  1  2  3
    vec_float4 b0 = b.v;                            // b0 =  4  5  6  7
    vec_float4 c1 = c.v;                            // c1 =  8  9 10 11
    vec_float4 d1 = d.v;                            // d1 = 12 13 14 15
  
    // Step 1: Interleave top and bottom half

    vec_float4 a1 = spu_shuffle( a0, c1, unpackl ); // a1 =  0  8  1  9
    vec_float4 b1 = spu_shuffle( b0, d1, unpackl ); // b1 =  4 12  5 13
    c1            = spu_shuffle( a0, c1, unpackh ); // c1 =  2 10  3 11
    d1            = spu_shuffle( b0, d1, unpackh ); // d1 =  6 14  7 15

    // Step 2: Interleave even and odd rows

    *((vec_float4 * ALIGNED(16))pa) = spu_shuffle( a1, b1, unpackl );
    *((vec_float4 * ALIGNED(16))pb) = spu_shuffle( a1, b1, unpackh );
    *((vec_float4 * ALIGNED(16))pc) = spu_shuffle( c1, d1, unpackl );
    *((vec_float4 * ALIGNED(16))pd) = spu_shuffle( c1, d1, unpackh );
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

    // v4float miscellaenous friends

    friend inline v4float clear_bits(  const v4int &m, const v4float &a );
    friend inline v4float set_bits(    const v4int &m, const v4float &a );
    friend inline v4float toggle_bits( const v4int &m, const v4float &a );

  public:

    // v4int constructors / destructors
    
    v4int() {}                              // Default constructor
    v4int( const v4int &a ) {               // Copy constructor
      v = a.v;
    }
    v4int( const v4 &a ) {                  // Initialize from mixed
      v = a.v;
    }
    v4int( const int &a ) {                 // Initialize from scalar
      v = (vec_float4)spu_splats( a );
    }
    v4int( const int &i0, const int &i1,
           const int &i2, const int &i3 ) { // Initialize from scalars
      // FIXME: vec_int4 t = { i0, i1, i2, i3 }; segfaults gcc
      i[0] = i0; i[1] = i1; i[2] = i2; i[3] = i3;
    }
    ~v4int() {}                             // Destructor
    
    // v4int assignment operators
  
#   define ASSIGN(op,instructions)                \
    inline v4int &operator op( const v4int &b ) { \
      instructions;                               \
      return *this;                               \
    }

    ASSIGN(=,   v  = b.v )
    ASSIGN(+=,  v  = (vec_float4)spu_add( (vec_int4)v, (vec_int4)b.v ) )
    ASSIGN(-=,  v  = (vec_float4)spu_sub( (vec_int4)v, (vec_int4)b.v ) )
    ASSIGN(*=,  i[0] *= b.i[0]; i[1] *= b.i[1];
                i[2] *= b.i[2]; i[3] *= b.i[3] ) // FIXME: Sigh ...
    ASSIGN(/=,  i[0] /= b.i[0]; i[1] /= b.i[1];
                i[2] /= b.i[2]; i[3] /= b.i[3] ) // FIXME: Sigh ...
    ASSIGN(%=,  i[0] %= b.i[0]; i[1] %= b.i[1];
                i[2] %= b.i[2]; i[3] %= b.i[3] ) // FIXME: Sigh ...
    ASSIGN(^=,  v = (vec_float4)spu_xor( (vec_int4)v, (vec_int4)b.v ) )
    ASSIGN(&=,  v = (vec_float4)spu_and( (vec_int4)v, (vec_int4)b.v ) )
    ASSIGN(|=,  v = (vec_float4)spu_or(  (vec_int4)v, (vec_int4)b.v ) )
    ASSIGN(<<=, i[0] <<= b.i[0]; i[1] <<= b.i[1];
                i[2] <<= b.i[2]; i[3] <<= b.i[3] ) // FIXME: Sigh ...
    ASSIGN(>>=, i[0] >>= b.i[0]; i[1] >>= b.i[1];
                i[2] >>= b.i[2]; i[3] >>= b.i[3] ) // FIXME: Sigh ...

#   undef ASSIGN

    // v4int member access operator
    
    inline int &operator []( const int n ) { return i[n]; }
    inline int  operator ()( const int n ) {
      return spu_extract( (vec_int4)v, n );
    }

  };

  // v4int prefix unary operators

# define PREFIX_UNARY(op,instructions)          \
  inline v4int operator op( const v4int & a ) { \
    v4int b;                                    \
    instructions;                               \
    return b;                                   \
  }

  PREFIX_UNARY(+, b.v = a.v )
  PREFIX_UNARY(-, b.v = (vec_float4)spu_sub( spu_splats( 0 ),
                                             (vec_int4)a.v ) )
  PREFIX_UNARY(!, b.v = (vec_float4)spu_cmpeq( spu_splats( 0 ),
                                               (vec_int4)a.v ) )
  PREFIX_UNARY(~, b.v = (vec_float4)spu_xor( spu_splats( -1 ),
                                             (vec_int4)a.v ) )
  // FIXME: Sigh
  
# undef PREFIX_UNARY

  // v4int prefix increment / decrement

# define PREFIX_INCDEC(op,intrinsic)                                    \
  inline v4int operator op( v4int & a ) {                               \
    vec_float4 a_v = a.v;                                               \
    v4int b;                                                            \
    a_v = (vec_float4)intrinsic( (vec_int4)a_v, spu_splats( 1 ) );      \
    a.v = a_v;                                                          \
    b.v = a_v;                                                          \
    return b;                                                           \
  }

  PREFIX_INCDEC(++,spu_add)
  PREFIX_INCDEC(--,spu_sub)

# undef PREFIX_INCDEC

  // v4int postfix increment / decrement

# define POSTFIX_INCDEC(op,intrinsic)                                   \
  inline v4int operator op( v4int & a, int ) {                          \
    vec_float4 a_v = a.v;                                               \
    v4int b;                                                            \
    b.v = a_v;                                                          \
    a.v = (vec_float4)intrinsic( (vec_int4)a_v, spu_splats( 1 ) );      \
    return b;                                                           \
  }

  POSTFIX_INCDEC(++,spu_add)
  POSTFIX_INCDEC(--,spu_sub)

# undef POSTFIX_INCDEC

  // v4int binary operators
  
# define BINARY(op,instructions)                                \
  inline v4int operator op( const v4int &a, const v4int &b ) {	\
    v4int c;                                                    \
    instructions;                                               \
    return c;                                                   \
  }

  BINARY(+,  c.v = (vec_float4)spu_add( (vec_int4)a.v, (vec_int4)b.v ) )
  BINARY(-,  c.v = (vec_float4)spu_sub( (vec_int4)a.v, (vec_int4)b.v ) )
  BINARY(*,  c.i[0] = a.i[0]*b.i[0]; c.i[1] = a.i[1]*b.i[1];
             c.i[2] = a.i[2]*b.i[2]; c.i[3] = a.i[3]*b.i[3] ) // FIXME: Sigh
  BINARY(/,  c.i[0] = a.i[0]/b.i[0]; c.i[1] = a.i[1]/b.i[1];
             c.i[2] = a.i[2]/b.i[2]; c.i[3] = a.i[3]/b.i[3] ) // FIXME: Sigh
  BINARY(%,  c.i[0] = a.i[0]%b.i[0]; c.i[1] = a.i[1]%b.i[1];
             c.i[2] = a.i[2]%b.i[2]; c.i[3] = a.i[3]%b.i[3] ) // FIXME: Sigh
  BINARY(^,  c.v = (vec_float4)spu_xor( (vec_int4)a.v, (vec_int4)b.v ) )
  BINARY(&,  c.v = (vec_float4)spu_and( (vec_int4)a.v, (vec_int4)b.v ) )
  BINARY(|,  c.v = (vec_float4)spu_or(  (vec_int4)a.v, (vec_int4)b.v ) )
  BINARY(<<, c.i[0] = a.i[0]<<b.i[0]; c.i[1] = a.i[1]<<b.i[1];
             c.i[2] = a.i[2]<<b.i[2]; c.i[3] = a.i[3]<<b.i[3] ) // FIXME: Sigh
  BINARY(>>, c.i[0] = a.i[0]>>b.i[0]; c.i[1] = a.i[1]>>b.i[1];
             c.i[2] = a.i[2]>>b.i[2]; c.i[3] = a.i[3]>>b.i[3] ) // FIXME: Sigh

# undef BINARY

  // v4int logical operators

# define LOGICAL(op,instructions)                              \
  inline v4int operator op( const v4int &a, const v4int &b ) { \
    v4int c;                                                   \
    instructions;                                              \
    return c;                                                  \
  }

  LOGICAL(<,  c.v = (vec_float4)spu_cmpgt( (vec_int4)b.v, (vec_int4)a.v ) )
  LOGICAL(>,  c.v = (vec_float4)spu_cmpgt( (vec_int4)a.v, (vec_int4)b.v ) )
  LOGICAL(==, c.v = (vec_float4)spu_cmpeq( (vec_int4)a.v, (vec_int4)b.v ) )
  LOGICAL(!=, c.v = (vec_float4)spu_xor( spu_splats( 0xffffffff ),
                                         spu_cmpeq( (vec_int4)a.v,
                                                    (vec_int4)b.v ) ) )
  /**/                                                           // FIXME: Sigh
  LOGICAL(<=, c.v = (vec_float4)spu_xor( spu_splats( 0xffffffff ),
                                         spu_cmpgt( (vec_int4)a.v,
                                                    (vec_int4)b.v ) ) )
  /**/                                                           // FIXME: Sigh
  LOGICAL(>=, c.v = (vec_float4)spu_xor( spu_splats( 0xffffffff ),
                                         spu_cmpgt( (vec_int4)b.v,
                                                    (vec_int4)a.v ) ) )
  /**/                                                           // FIXME: Sigh
  LOGICAL(&&, c.v = (vec_float4)spu_nor( spu_cmpeq( (vec_int4)a.v,
                                                    spu_splats( 0 ) ),
                                         spu_cmpeq( (vec_int4)b.v,
                                                    spu_splats( 0 ) ) ) )
  LOGICAL(||, c.v = (vec_float4)spu_nand( spu_cmpeq( (vec_int4)a.v,
                                                     spu_splats( 0 ) ),
                                          spu_cmpeq( (vec_int4)b.v,
                                                     spu_splats( 0 ) ) ) )

# undef LOGICAL

  // v4int miscellaneous functions

  inline v4int abs( const v4int &a ) { // FIXME: This could be done better
    v4int b;
    b.i[0] = (a.i[0]>=0) ? a.i[0] : -a.i[0];
    b.i[1] = (a.i[1]>=0) ? a.i[1] : -a.i[1];
    b.i[2] = (a.i[2]>=0) ? a.i[2] : -a.i[2];
    b.i[3] = (a.i[3]>=0) ? a.i[3] : -a.i[3];
    return b;
  }

  inline v4 czero( const v4int &c, const v4 &a ) {
    v4 b;
    b.v = spu_andc( a.v, c.v );
    return b;
  }

  inline v4 notczero( const v4int &c, const v4 &a ) {
    v4 b;
    b.v = spu_and( a.v, c.v );
    return b;
  }
  
  inline v4 merge( const v4int &c, const v4 &t, const v4 &f ) {
    v4 m;
    m.v = spu_sel( f.v, t.v, (vec_uint4)c.v );
    return m;
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
    friend inline v4float fma(  const v4float &a,
                                const v4float &b,
                                const v4float &c );
    friend inline v4float fms(  const v4float &a,
                                const v4float &b,
                                const v4float &c );
    friend inline v4float fnms( const v4float &a,
                                const v4float &b,
                                const v4float &c );
    friend inline v4float clear_bits(  const v4int &m, const v4float &a );
    friend inline v4float set_bits(    const v4int &m, const v4float &a );
    friend inline v4float toggle_bits( const v4int &m, const v4float &a );
    friend inline void increment_4x1( float * ALIGNED(16) p,
                                      const v4float &a );
    friend inline void decrement_4x1( float * ALIGNED(16) p,
                                      const v4float &a );
    friend inline void scale_4x1(     float * ALIGNED(16) p,
                                      const v4float &a );
    // FIXME: crack
    
  public:

    // v4float constructors / destructors
    
    v4float() {}                                  // Default constructor
    v4float( const v4float &a ) {                 // Copy constructor
      v = a.v;
    }
    v4float( const v4 &a ) {                      // Initialize from mixed
      v = a.v;
    }
    v4float( const float &a ) {                   // Initialize from scalar
      v = spu_splats(a);
    }
    v4float( const float &f0, const float &f1,
             const float &f2, const float &f3 ) { // Initalize from scalars
      // FIXME: vec_float4 t = { f0, f1, f2, f3 };
      // ... seg faults gcc .. MAYBE NOT ANYMORE
      f[0] = f0; f[1] = f1; f[2] = f2; f[3] = f3;
    }
    ~v4float() {}                                 // Destructor

    // v4float assignment operators

#   define ASSIGN(op,intrinsic)                         \
    inline v4float &operator op( const v4float &b ) {	\
      v = intrinsic( v, b.v );                          \
      return *this;                                     \
    }

    inline v4float &operator =( const v4float &b ) {
      v = b.v;
      return *this;
    }

    ASSIGN(+=,spu_add)
    ASSIGN(-=,spu_sub)
    ASSIGN(*=,spu_mul)

    inline v4float &operator /=( const v4float &a ) {
      vec_float4 a_v = a.v, b_v;

      // Compute an estimate of the reciprocal of a (12-bit accurate)

      b_v = spu_re( a_v );
     
      // FIXME: CHECK NUMERICS ... MAY WANT ADDITIONAL N-R STEPS:
      // FIXME: THE TWO MUL-ADD BASED REFINMENT SEEMS LESS ACCURATE WHEN
      // USED IN THIS CONTEXT!
      // b_v = spu_nmsub( a_v, spu_mul( b_v, b_v ), spu_add( b_v, b_v ) );

      // Compute a * refined( (1/b)_estimate ) to get result a/b

      v = spu_mul( v, spu_nmsub( a_v, spu_mul( b_v, b_v ),
                                 spu_add( b_v, b_v ) ) );

      return *this;
    }

#   undef ASSIGN

    // v4float member access operator

    inline float &operator []( const int n ) { return f[n]; }
    inline float  operator ()( const int n ) { return spu_extract( v, n ); }

  };

  // v4float prefix unary operators

  inline v4float operator +( const v4float &a ) {
    v4float b;
    b.v = a.v;
    return b;
  }

  inline v4float operator -( const v4float &a ) {
    v4float b;
    b.v = spu_sub( spu_splats( 0.f ), a.v );
    return b;
  }

  inline v4int operator !( const v4float &a ) {
    v4int b;
    b.v = (vec_float4)spu_cmpeq( a.v, spu_splats( 0.f ) );
    return b;
  }

  // v4float prefix increment / decrement operators

  inline v4float operator ++( v4float &a ) {
    vec_float4 a_v = a.v;
    v4float b;
    a_v = spu_add( a_v, spu_splats( 1.f ) );
    a.v = a_v;
    b.v = a_v;
    return b;
  }

  inline v4float operator --( v4float &a ) {
    vec_float4 a_v = a.v;
    v4float b;
    a_v = spu_sub( a_v, spu_splats( 1.f ) );
    a.v = a_v;
    b.v = a_v;
    return b;
  }

  // v4float postfix increment / decrement operators

  inline v4float operator ++( v4float &a, int ) {
    vec_float4 a_v = a.v;
    v4float b;
    b.v = a_v;
    a_v = spu_add( a_v, spu_splats( 1.f ) );
    a.v = a_v;
    return b;
  }

  inline v4float operator --( v4float &a, int ) {
    vec_float4 a_v = a.v;
    v4float b;
    b.v = a_v;
    a_v = spu_sub( a_v, spu_splats( 1.f ) );
    a.v = a_v;
    return b;
  }

  // v4float binary operators
    
# define BINARY(op,intrinsic)                                        \
  inline v4float operator op( const v4float &a, const v4float &b ) { \
    v4float c;                                                       \
    c.v = intrinsic( a.v, b.v );                                     \
    return c;                                                        \
  }

  BINARY(+,spu_add)
  BINARY(-,spu_sub)
  BINARY(*,spu_mul)

  inline v4float operator /( const v4float &n, const v4float &a ) {
    vec_float4 a_v = a.v, b_v;
    v4float c;
    
    // Compute an estimate of the reciprocal of a (12-bit accurate)
    
    b_v = spu_re( a_v );
    
    // FIXME: CHECK NUMERICS ... MAY WANT ADDITIONAL N-R STEPS:
    // FIXME: THE TWO MUL-ADD BASED REFINMENT SEEMS LESS ACCURATE WHEN
    // USED IN THIS CONTEXT!
    // b_v = spu_nmsub( a_v, spu_mul( b_v, b_v ), spu_add( b_v, b_v ) );

    // Compute a * refined( (1/b)_estimate ) to get result a/b
    
    c.v = spu_mul( n.v, spu_nmsub( a_v, spu_mul( b_v, b_v ),
                                   spu_add( b_v, b_v ) ) );

    return c;
  }

# undef BINARY

  // v4float logical operators

# define LOGICAL(op,assembly)                                      \
  inline v4int operator op( const v4float &a, const v4float &b ) { \
    v4int c;                                                       \
    c.v = assembly;                                                \
    return c;                                                      \
  }

  LOGICAL(< ,(vec_float4)spu_cmpgt( b.v, a.v ))
  LOGICAL(> ,(vec_float4)spu_cmpgt( a.v, b.v ))
  LOGICAL(==,(vec_float4)spu_cmpeq( a.v, b.v ))
  LOGICAL(!=,(vec_float4)spu_xor( spu_splats( 0xffffffff ),
                                  spu_cmpeq( a.v, b.v ) ))
  // FIXME: ALMOST CERTAINLY NOT IEEE.  MAY BE A BETTER WAY.
  LOGICAL(<=,(vec_float4)spu_xor( spu_splats( 0xffffffff ),
                                  spu_cmpgt( a.v, b.v ) )) // FIXME: SEE ABOVE
  LOGICAL(>=,(vec_float4)spu_xor( spu_splats( 0xffffffff ),
                                  spu_cmpgt( b.v, a.v ) )) // FIXME: SEE ABOVE
  // FIXME: Does logical true for a floating point just mean any bits high
  // or is it more subtle (i.e. logical true for a floating point number
  // could mean IEEE compare to unsigned zero is false).
  LOGICAL(&&,(vec_float4)spu_nor(  spu_cmpeq( a.v, spu_splats( 0.f ) ),
                                   spu_cmpeq( b.v, spu_splats( 0.f ) ) ))
  LOGICAL(||,(vec_float4)spu_nand( spu_cmpeq( a.v, spu_splats( 0.f ) ),
                                   spu_cmpeq( b.v, spu_splats( 0.f ) ) ))

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
    b.v = spu_andc( a.v, (vec_float4)spu_splats( 1<<31 ) );
    return b;            
  }

  inline v4float sqrt( const v4float &a ) {
    vec_float4 a_v = a.v, b_v;
    vec_float4 half = spu_splats( 0.5f ), one = spu_splats( 1.f );
    v4float b;

    // Compute an estimate of the rsqrt (12-bit accurate)

    b_v = spu_rsqrte( a_v );

    // FIXME: CHECK NUMERICS.  MIGHT WANT TO USE ADDITIONAL REFINEMENT:
    // b_v = spu_madd( spu_nmsub( spu_mul( b_v, b_v ), a_v, one ),
    //                            spu_mul( b_v, half ),
    //                            b_v );

    // Compute the sqrt(a) via a*refined_rsqrt_estimate(a) ~ sqrt(a)

    b.v = spu_mul( a_v, spu_madd( spu_nmsub( spu_mul( b_v, b_v ), a_v, one ),
                                  spu_mul( b_v, half ),
                                  b_v ) );

    return b;
  }

  inline v4float copysign( const v4float &a, const v4float &b ) {
    v4float c;
    c.v = spu_or( spu_andc( a.v, (vec_float4)spu_splats( 1<<31 ) ),
                  spu_and(  b.v, (vec_float4)spu_splats( 1<<31 ) ) );
    return c;
  }

# undef CMATH_FR1
# undef CMATH_FR2

  // v4float miscelleanous functions
  
  inline v4float rsqrt_approx( const v4float &a ) {
    vec_float4 a_v = a.v;
    v4float b;
    b.v = spu_rsqrte( a_v );
    return b;
  }
  
  inline v4float rsqrt( const v4float &a ) {
    vec_float4 a_v = a.v, b_v;
    vec_float4 half = spu_splats( 0.5f ), one = spu_splats( 1.f );
    v4float b;

    // Compute an estimate of the rsqrt (12-bit accurate)

    b_v = spu_rsqrte( a_v );

    // Refine the estimate with N-R.
    // FIXME: CHECK NUMERICS.  MIGHT WANT TO USE ADDITIONAL REFINEMENT.

    b.v = spu_madd( spu_nmsub( spu_mul( b_v, b_v ), a_v, one ),
                    spu_mul( b_v, half ),
                    b_v );
    
    return b;
  }

  inline v4float rcp_approx( const v4float &a ) {
    vec_float4 a_v = a.v;
    v4float b;
    b.v = spu_re( a_v );
    return b;
  }
  
  inline v4float rcp( const v4float &a ) {
    vec_float4 a_v = a.v, b_v, one = spu_splats( 1.f );
    v4float b;

    // Compute an estimate of the reciprocal of a (12-bit accurate)

    b_v = spu_re( a_v );
    
    // Perform Newton-Raphson refinement (one step should give 24-bit)
    // FIXME: CHECK NUMERICS ... MAY WANT TO USE ADDITIONAL REFINEMNT.
    
    b.v = spu_madd( spu_nmsub( b_v, a_v, one ), b_v, b_v );

    return b;
  }

  inline v4float fma(  const v4float &a, const v4float &b, const v4float &c ) {
    v4float d;
    d.v = spu_madd( a.v, b.v, c.v );
    return d;
  }

  inline v4float fms(  const v4float &a, const v4float &b, const v4float &c ) {
    v4float d;
    d.v = spu_msub( a.v, b.v, c.v );
    return d;
  }

  inline v4float fnms( const v4float &a, const v4float &b, const v4float &c ) {
    v4float d;
    d.v = spu_nmsub( a.v, b.v, c.v );
    return d;
  }

  inline v4float clear_bits( const v4int &m, const v4float &a ) {
    v4float b;
    b.v = spu_andc( a.v, m.v );
    return b;
  }

  inline v4float set_bits( const v4int &m, const v4float &a ) {
    v4float b;
    b.v = spu_or( a.v, m.v );
    return b;
  }

  inline v4float toggle_bits( const v4int &m, const v4float &a ) {
    v4float b;
    b.v = spu_xor( a.v, m.v );
    return b;
  }

  inline void increment_4x1( float * ALIGNED(16) p, const v4float &a ) {
    *((vec_float4 * ALIGNED(16))p) =
      spu_add( *((vec_float4 * ALIGNED(16))p), a.v );
  }

  inline void decrement_4x1( float * ALIGNED(16) p, const v4float &a ) {
    *((vec_float4 * ALIGNED(16))p) =
      spu_sub( *((vec_float4 * ALIGNED(16))p), a.v );
  }

  inline void scale_4x1( float * ALIGNED(16) p, const v4float &a ) {
    *((vec_float4 * ALIGNED(16))p) =
      spu_mul( *((vec_float4 * ALIGNED(16))p), a.v );
  }

} // namespace v4

#endif // _v4_spu_hxx_
