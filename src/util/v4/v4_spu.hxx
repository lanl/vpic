#ifndef _v4_spu_hxx_
#define _v4_spu_hxx_

#ifndef IN_v4_h
#error "Do not include v4_spu.hxx directly; use v4.h"
#endif

#define V4_ACCELERATION
#define V4_SPU_ACCELERATION

#include <spu_intrinsics.h> /* Requires -Wno-long-long to include */
#include <stdint.h>
#include <math.h>

#ifndef ALIGNED
#define ALIGNED
#endif

namespace v4 {

  // FIXME: CHECK VECTOR CASTING BETWEEN TYPES IS FREE IN COMPILER GENERATED
  // ASSEMBLY

  // FIXME: IT MAY BE FASTER WAY TO ASSEMBLY THESE VECTOR CONSTANTS ON THE FLY
  // (e.g spu_splats)

  const vec_uint4  _vufalse = {          0,          0,
                                         0,          0 };
  const vec_uint4  _vutrue  = { 0xffffffff, 0xffffffff,
                                0xffffffff, 0xffffffff };

  const vec_int4   _vizero  = {  0,  0,  0,  0 };
  const vec_int4   _vione   = {  1,  1,  1,  1 };

  const vec_float4 _vfzero  = {   0,   0,   0,   0 };
  const vec_float4 _vfhalf  = { 0.5, 0.5, 0.5, 0.5 };
  const vec_float4 _vfone   = {   1,   1,   1,   1 };


  const vec_uchar16 _tr0 = {  0, 1, 2, 3,    4, 5, 6, 7,
                             16,17,18,19,   20,21,22,23 };
  const vec_uchar16 _tr1 = {  8, 9,10,11,   12,13,14,15,
                             24,25,26,27,   28,29,30,31 };
  const vec_uchar16 _tr2 = {  0, 1, 2, 3,   16,17,18,19,
                              8, 9,10,11,   24,25,26,27 };
  const vec_uchar16 _tr3 = {  4, 5, 6, 7,   20,21,22,23,
                             12,13,14,15,   28,29,30,31 };
  const vec_uchar16 _tr4 = {  0, 1, 2, 3,    8, 9,10,11,
                             16,17,18,19,   24,25,26,27 };
  const vec_uchar16 _tr5 = {  4, 5, 6, 7,   12,13,14,15,
                             20,21,22,23,   28,29,30,31 };

  const vec_uint4   _visignbits = { 1<<31, 1<<31, 1<<31, 1<<31 };

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
        
    friend inline void load_4x1( const void * ALIGNED p, v4 &a );
    friend inline void store_4x1( const v4 &a, void * ALIGNED p );
    friend inline void stream_4x1( const v4 &a, void * ALIGNED p );
    friend inline void copy_4x1( void * ALIGNED dst,
                                 const void * ALIGNED src );
    friend inline void swap_4x1( void * ALIGNED a, void * ALIGNED b );

    // v4 transposed memory manipulation friends
    // Note: Half aligned values are permissible in the 4x2_tr variants!

    friend inline void load_4x1_tr( const void *a0, const void *a1,
                                    const void *a2, const void *a3,
                                    v4 &a );
    friend inline void load_4x2_tr( const void * ALIGNED a0,
                                    const void * ALIGNED a1,
                                    const void * ALIGNED a2,
                                    const void * ALIGNED a3,
                                    v4 &a, v4 &b );
    friend inline void load_4x3_tr( const void * ALIGNED a0,
                                    const void * ALIGNED a1,
                                    const void * ALIGNED a2,
                                    const void * ALIGNED a3,
                                    v4 &a, v4 &b, v4 &c );
    friend inline void load_4x4_tr( const void * ALIGNED a0,
                                    const void * ALIGNED a1,
                                    const void * ALIGNED a2,
                                    const void * ALIGNED a3,
                                    v4 &a, v4 &b, v4 &c, v4 &d );
    
    friend inline void store_4x1_tr( const v4 &a,
                                     void *a0, void *a1, void *a2, void *a3 );
    friend inline void store_4x2_tr( const v4 &a, const v4 &b,
                                     void * ALIGNED a0, void * ALIGNED a1,
                                     void * ALIGNED a2, void * ALIGNED a3 );
    friend inline void store_4x3_tr( const v4 &a, const v4 &b, const v4 &c,
                                     void * ALIGNED a0, void * ALIGNED a1,
                                     void * ALIGNED a2, void * ALIGNED a3 );
    friend inline void store_4x4_tr( const v4 &a, const v4 &b,
                                     const v4 &c, const v4 &d,
                                     void * ALIGNED a0, void * ALIGNED a1,
                                     void * ALIGNED a2, void * ALIGNED a3 );

  protected:

    union {
      int i[4];
      float f[4];
      vec_float4 vf;
      vec_int4   vi; // FIXME: ALWAYS USE vf WITH CASTING??
    };
    
  public:

    v4() {}                        // Default constructor
    v4(const v4 &a) { vi = a.vi; } // Copy constructor
    ~v4() {}                       // Default destructor

  };
  
  // v4 miscellaneous functions

  inline int any( const v4 &a ) {
    return
      spu_extract( spu_gather( spu_cmpeq( a.vi, _vizero ) ), 0 )==15 ? 0 : -1; 
  }
  
  inline int all( const v4 &a ) {
    return
      spu_extract( spu_gather( spu_cmpeq( a.vi, _vizero ) ), 0 )==0 ? -1 : 0;
  }
  
  inline v4 splat( const v4 & a, const int n ) {
    v4 b;
    b.vi = spu_splats( spu_extract( a.vi, n ) ); // FIXME: A better way?
    return b;
  }

  inline void swap( v4 &a, v4 &b ) { 
    vec_int4 t = a.vi; a.vi = b.vi; b.vi = t;
  }

  inline void transpose( v4 &a, v4 &b, v4 &c, v4 &d ) {
    vec_int4 a0 = a.vi;                        // a0 =  0  1  2  3
    vec_int4 b0 = b.vi;                        // b0 =  4  5  6  7
    vec_int4 c0 = c.vi;                        // c0 =  8  9 10 11
    vec_int4 d0 = d.vi;                        // d0 = 12 13 14 15

    // Step 1: Transpose the block matrix

    vec_int4 a1 = spu_shuffle( a0, c0, _tr0 ); // a1 =  0  1  8  9
    vec_int4 b1 = spu_shuffle( b0, d0, _tr0 ); // b1 =  4  5 12 13
    vec_int4 c1 = spu_shuffle( a0, c0, _tr1 ); // c1 =  2  3 10 11
    vec_int4 d1 = spu_shuffle( b0, d0, _tr1 ); // d1 =  6  7 14 15

    // Step 2: Transpose 2x2 subblocks of matrix

    a.vi = spu_shuffle( a1, b1, _tr2 );        // a  =  0  4  8 12
    b.vi = spu_shuffle( a1, b1, _tr3 );        // b  =  1  5  9 13
    c.vi = spu_shuffle( c1, d1, _tr2 );        // c  =  2  6 10 14
    d.vi = spu_shuffle( c1, d1, _tr3 );        // d  =  3  7 11 15
  }

  // v4 memory manipulation functions
  
  inline void load_4x1( const void * ALIGNED p, v4 &a ) {
    a.vi = *((const vec_int4 * ALIGNED)p);
  }

  inline void store_4x1( const v4 &a, void * ALIGNED p ) {
    *((vec_int4 * ALIGNED)p) = a.vi;
  }

  inline void stream_4x1( const v4 &a, void * ALIGNED p ) {
    *((vec_int4 * ALIGNED)p) = a.vi;
  }

  // FIXME: Ordering semantics
  inline void copy_4x1( void * ALIGNED dst, const void * ALIGNED src ) {
    *((vec_int4 * ALIGNED)dst) = *((const vec_int4 * ALIGNED)src);
  }

  inline void swap_4x1( void * ALIGNED a, void * ALIGNED b ) {
    vec_int4 t               = *((vec_int4 * ALIGNED)a);
    *((vec_int4 * ALIGNED)a) = *((vec_int4 * ALIGNED)b);
    *((vec_int4 * ALIGNED)b) = t;
  }

  // v4 transposed memory manipulation functions

  inline void load_4x1_tr( const void *a0, const void *a1,
                           const void *a2, const void *a3, v4 &a ) {
    a.i[0] = ((const int *)a0)[0];
    a.i[1] = ((const int *)a1)[0];
    a.i[2] = ((const int *)a2)[0];
    a.i[3] = ((const int *)a3)[0];
  }

  // Note: load_4x4_tr with last two shuffles discared may be
  // preferable for 128-byte aligned results.  (Compiler optimization
  // may actually do that for free!)

  inline void load_4x2_tr( const void * ALIGNED pa, const void * ALIGNED pb,
                           const void * ALIGNED pc, const void * ALIGNED pd,
                           v4 &a, v4 &b ) {
    vec_llong2 a_v = { *(const int64_t * ALIGNED)pa,
                       *(const int64_t * ALIGNED)pb }; // 0 4 1 5
    vec_llong2 b_v = { *(const int64_t * ALIGNED)pc,
                       *(const int64_t * ALIGNED)pd }; // 2 6 3 7
    a.vi = (vec_int4)spu_shuffle( a_v, b_v, _tr4 ); // 0 1 2 3
    b.vi = (vec_int4)spu_shuffle( a_v, b_v, _tr5 ); // 4 5 6 7
  }
  
  inline void load_4x3_tr( const void * ALIGNED pa, const void * ALIGNED pb,
                           const void * ALIGNED pc, const void * ALIGNED pd,
                           v4 &a, v4 &b, v4 &c ) {
    vec_int4 a0 = *((const vec_int4 * ALIGNED)pa);
    vec_int4 b0 = *((const vec_int4 * ALIGNED)pb);
    vec_int4 c0 = *((const vec_int4 * ALIGNED)pc);
    vec_int4 d0 = *((const vec_int4 * ALIGNED)pd);

    // Step 1: Transpose the block matrix

    vec_int4 a1 = spu_shuffle( a0, c0, _tr0 ); // a1 =  0  1  8  9
    vec_int4 b1 = spu_shuffle( b0, d0, _tr0 ); // b1 =  4  5 12 13
    vec_int4 c1 = spu_shuffle( a0, c0, _tr1 ); // c1 =  2  3 10 11
    vec_int4 d1 = spu_shuffle( b0, d0, _tr1 ); // d1 =  6  7 14 15

    // Step 2: Transpose 2x2 subblocks of matrix

    a.vi = spu_shuffle( a1, b1, _tr2 );        // a  =  0  4  8 12
    b.vi = spu_shuffle( a1, b1, _tr3 );        // b  =  1  5  9 13
    c.vi = spu_shuffle( c1, d1, _tr2 );        // c  =  2  6 10 14
  }

  inline void load_4x4_tr( const void * ALIGNED pa, const void * ALIGNED pb,
                           const void * ALIGNED pc, const void * ALIGNED pd,
                           v4 &a, v4 &b, v4 &c, v4 &d ) {
    vec_int4 a0 = *((const vec_int4 * ALIGNED)pa);
    vec_int4 b0 = *((const vec_int4 * ALIGNED)pb);
    vec_int4 c0 = *((const vec_int4 * ALIGNED)pc);
    vec_int4 d0 = *((const vec_int4 * ALIGNED)pd);

    // Step 1: Transpose the block matrix

    vec_int4 a1 = spu_shuffle( a0, c0, _tr0 ); // a1 =  0  1  8  9
    vec_int4 b1 = spu_shuffle( b0, d0, _tr0 ); // b1 =  4  5 12 13
    vec_int4 c1 = spu_shuffle( a0, c0, _tr1 ); // c1 =  2  3 10 11
    vec_int4 d1 = spu_shuffle( b0, d0, _tr1 ); // d1 =  6  7 14 15

    // Step 2: Transpose 2x2 subblocks of matrix

    a.vi = spu_shuffle( a1, b1, _tr2 );        // a  =  0  4  8 12
    b.vi = spu_shuffle( a1, b1, _tr3 );        // b  =  1  5  9 13
    c.vi = spu_shuffle( c1, d1, _tr2 );        // c  =  2  6 10 14
    d.vi = spu_shuffle( c1, d1, _tr3 );        // d  =  3  7 11 15
  }

  inline void store_4x1_tr( const v4 &a,
                            void *a0, void *a1, void *a2, void *a3 ) {
    ((int *)a0)[0] = a.i[0];
    ((int *)a1)[0] = a.i[1];
    ((int *)a2)[0] = a.i[2];
    ((int *)a3)[0] = a.i[3];
  }

  inline void store_4x2_tr( const v4 &a, const v4 &b,
                            void * ALIGNED pa, void * ALIGNED pb,
                            void * ALIGNED pc, void * ALIGNED pd ) {
    vec_int4 a1 = a.vi;                        // a =  0  1  2  3
    vec_int4 b1 = b.vi;                        // b =  4  5  6  7
    vec_int4 c1 = spu_shuffle( a1, a1, _tr1 ); // c =  2  3  x  x
    vec_int4 d1 = spu_shuffle( b1, b1, _tr1 ); // d =  6  7  x  x

    // Step 2: Transpose 2x2 subblocks of matrix

    vec_int4 a0 = spu_shuffle( a1, b1, _tr2 ); // a  =  0  4  x  x
    vec_int4 b0 = spu_shuffle( a1, b1, _tr3 ); // b  =  1  5  x  x
    vec_int4 c0 = spu_shuffle( c1, d1, _tr2 ); // c  =  2  6  x  x
    vec_int4 d0 = spu_shuffle( c1, d1, _tr3 ); // d  =  3  7  x  x

    // Store the 2 columns of the matrix

    ((int64_t * ALIGNED)pa)[0] = spu_extract( (vec_llong2)a0, 0 );
    ((int64_t * ALIGNED)pb)[0] = spu_extract( (vec_llong2)b0, 0 );
    ((int64_t * ALIGNED)pc)[0] = spu_extract( (vec_llong2)c0, 0 );
    ((int64_t * ALIGNED)pd)[0] = spu_extract( (vec_llong2)d0, 0 );
  }

  inline void store_4x3_tr( const v4 &a, const v4 &b, const v4 &c,
                            void * ALIGNED pa, void * ALIGNED pb,
                            void * ALIGNED pc, void * ALIGNED pd ) {
    vec_int4 a0 = a.vi;                        // a =  0  1  2  3
    vec_int4 b0 = b.vi;                        // b =  4  5  6  7
    vec_int4 c0 = c.vi;                        // c =  8  9 10 11
    vec_int4 d0;                               // d =  x  x  x  x ... no warn
    
    // Step 1: Transpose the block matrix

    vec_int4 a1 = spu_shuffle( a0, c0, _tr0 ); // a =  0  1  8  9
    vec_int4 b1 = spu_shuffle( b0, b0, _tr0 ); // b =  4  5  x  x
    vec_int4 c1 = spu_shuffle( a0, c0, _tr1 ); // c =  2  3 10 11
    vec_int4 d1 = spu_shuffle( b0, b0, _tr1 ); // d =  6  7  x  x

    // Step 2: Transpose 2x2 subblocks of matrix

    a0 = spu_shuffle( a1, b1, _tr2 );          // a  =  0  4  8  x
    b0 = spu_shuffle( a1, b1, _tr3 );          // b  =  1  5  9  x
    c0 = spu_shuffle( c1, d1, _tr2 );          // c  =  2  6 10  x
    d0 = spu_shuffle( c1, d1, _tr3 );          // d  =  3  7 11  x

    // Store the 3 columns of the matrix

    ((int64_t * ALIGNED)pa)[0] = spu_extract( (vec_llong2)a0, 0 );
    ((int64_t * ALIGNED)pb)[0] = spu_extract( (vec_llong2)b0, 0 );
    ((int64_t * ALIGNED)pc)[0] = spu_extract( (vec_llong2)c0, 0 );
    ((int64_t * ALIGNED)pd)[0] = spu_extract( (vec_llong2)d0, 0 );

    ((int * ALIGNED)pa)[2]       = spu_extract( a0, 2 );
    ((int * ALIGNED)pb)[2]       = spu_extract( b0, 2 );
    ((int * ALIGNED)pc)[2]       = spu_extract( c0, 2 );
    ((int * ALIGNED)pd)[2]       = spu_extract( d0, 2 );
  }
  
  inline void store_4x4_tr( const v4 &a, const v4 &b, const v4 &c, const v4 &d,
                            void * ALIGNED pa, void * ALIGNED pb,
                            void * ALIGNED pc, void * ALIGNED pd ) {
    vec_int4 a0 = a.vi;
    vec_int4 b0 = b.vi;
    vec_int4 c0 = c.vi;
    vec_int4 d0 = d.vi;

    // Step 1: Transpose the block matrix

    vec_int4 a1 = spu_shuffle( a0, c0, _tr0 ); // a1 =  0  1  8  9
    vec_int4 b1 = spu_shuffle( b0, d0, _tr0 ); // b1 =  4  5 12 13
    vec_int4 c1 = spu_shuffle( a0, c0, _tr1 ); // c1 =  2  3 10 11
    vec_int4 d1 = spu_shuffle( b0, d0, _tr1 ); // d1 =  6  7 14 15

    // Step 2: Transpose 2x2 subblocks of matrix

    *((vec_int4 * ALIGNED)pa) = spu_shuffle( a1, b1, _tr2 );
    *((vec_int4 * ALIGNED)pb) = spu_shuffle( a1, b1, _tr3 );
    *((vec_int4 * ALIGNED)pc) = spu_shuffle( c1, d1, _tr2 );
    *((vec_int4 * ALIGNED)pd) = spu_shuffle( c1, d1, _tr3 );
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
      vi = a.vi;
    }
    v4int( const v4 &a ) {                  // Initialize from mixed
      vi = a.vi;
    }
    v4int( const int &a ) {                 // Initialize from scalar
      vi = spu_splats( a );
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

    ASSIGN(=,   vi = b.vi )
    ASSIGN(+=,  vi = spu_add( vi, b.vi ) )
    ASSIGN(-=,  vi = spu_sub( vi, b.vi ) )
    ASSIGN(*=,  i[0] *= b.i[0]; i[1] *= b.i[1];
                i[2] *= b.i[2]; i[3] *= b.i[3] ) // FIXME: Sigh ...
    ASSIGN(/=,  i[0] /= b.i[0]; i[1] /= b.i[1];
                i[2] /= b.i[2]; i[3] /= b.i[3] ) // FIXME: Sigh ...
    ASSIGN(%=,  i[0] %= b.i[0]; i[1] %= b.i[1];
                i[2] %= b.i[2]; i[3] %= b.i[3] ) // FIXME: Sigh ...
    ASSIGN(^=,  vi = spu_xor( vi, b.vi ) )
    ASSIGN(&=,  vi = spu_and( vi, b.vi ) )
    ASSIGN(|=,  vi = spu_or(  vi, b.vi ) )
    ASSIGN(<<=, i[0] <<= b.i[0]; i[1] <<= b.i[1];
                i[2] <<= b.i[2]; i[3] <<= b.i[3] ) // FIXME: Sigh ...
    ASSIGN(>>=, i[0] >>= b.i[0]; i[1] >>= b.i[1];
                i[2] >>= b.i[2]; i[3] >>= b.i[3] ) // FIXME: Sigh ...

#   undef ASSIGN

    // v4int member access operator
    
    int &operator()(const int n) {
      return i[n];
    }

  };

  // v4int prefix unary operators

# define PREFIX_UNARY(op,instructions)          \
  inline v4int operator op( const v4int & a ) { \
    v4int b;                                    \
    instructions;                               \
    return b;                                   \
  }

  PREFIX_UNARY(+, b.vi = a.vi )
  PREFIX_UNARY(-, b.vi = spu_sub( _vizero, a.vi ) )
  PREFIX_UNARY(!, b.vi = (vec_int4)spu_cmpeq( _vizero, a.vi ) )
  PREFIX_UNARY(~, b.vi = spu_xor( (vec_int4)_vutrue, a.vi ) ) // FIXME: Sigh
  
# undef PREFIX_UNARY

  // v4int prefix increment / decrement

# define PREFIX_INCDEC(op,intrinsic)            \
  inline v4int operator op( v4int & a ) {       \
    vec_int4 a_vi = a.vi;                       \
    v4int b;                                    \
    a_vi = intrinsic( a_vi, _vione );           \
    a.vi = a_vi;                                \
    b.vi = a_vi;                                \
    return b;                                   \
  }

  PREFIX_INCDEC(++,spu_add)
  PREFIX_INCDEC(--,spu_sub)

# undef PREFIX_INCDEC

  // v4int postfix increment / decrement

# define POSTFIX_INCDEC(op,intrinsic)           \
  inline v4int operator op( v4int & a, int ) {  \
    vec_int4 a_vi = a.vi;                     \
    v4int b;                                    \
    b.vi = a_vi;                                \
    a.vi = intrinsic( a_vi, _vione );           \
    return b;                                   \
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

  BINARY(+,  c.vi = spu_add( a.vi, b.vi ) )
  BINARY(-,  c.vi = spu_sub( a.vi, b.vi ) )
  BINARY(*,  c.i[0] = a.i[0]*b.i[0]; c.i[1] = a.i[1]*b.i[1];
             c.i[2] = a.i[2]*b.i[2]; c.i[3] = a.i[3]*b.i[3] ) // FIXME: Sigh
  BINARY(/,  c.i[0] = a.i[0]/b.i[0]; c.i[1] = a.i[1]/b.i[1];
             c.i[2] = a.i[2]/b.i[2]; c.i[3] = a.i[3]/b.i[3] ) // FIXME: Sigh
  BINARY(%,  c.i[0] = a.i[0]%b.i[0]; c.i[1] = a.i[1]%b.i[1];
             c.i[2] = a.i[2]%b.i[2]; c.i[3] = a.i[3]%b.i[3] ) // FIXME: Sigh
  BINARY(^,  c.vi = spu_xor( a.vi, b.vi ) )
  BINARY(&,  c.vi = spu_and( a.vi, b.vi ) )
  BINARY(|,  c.vi = spu_or(  a.vi, b.vi ) )
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

  LOGICAL(<,  c.vi = (vec_int4)spu_cmpgt( b.vi, a.vi ) )
  LOGICAL(>,  c.vi = (vec_int4)spu_cmpgt( a.vi, b.vi ) )
  LOGICAL(==, c.vi = (vec_int4)spu_cmpeq( a.vi, b.vi ) )
  LOGICAL(!=, c.vi = (vec_int4)spu_xor( _vutrue, spu_cmpeq( a.vi, b.vi ) ) ) // FIXME: Sigh
  LOGICAL(<=, c.vi = (vec_int4)spu_xor( _vutrue, spu_cmpgt( a.vi, b.vi ) ) ) // FIXME: Sigh
  LOGICAL(>=, c.vi = (vec_int4)spu_xor( _vutrue, spu_cmpgt( b.vi, a.vi ) ) ) // FIXME: Sigh
  LOGICAL(&&, c.vi = (vec_int4)spu_nor(  spu_cmpeq( a.vi, _vizero ),
                                         spu_cmpeq( b.vi, _vizero ) ) )
  LOGICAL(||, c.vi = (vec_int4)spu_nand( spu_cmpeq( a.vi, _vizero ),
                                         spu_cmpeq( b.vi, _vizero ) ) )
  
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
    b.vi = spu_andc( a.vi, c.vi );
    return b;
  }

  inline v4 notczero( const v4int &c, const v4 &a ) {
    v4 b;
    b.vi = spu_and( a.vi, c.vi );
    return b;
  }
  
  inline v4 merge( const v4int &c, const v4 &t, const v4 &f ) {
    v4 m;
    m.vi = spu_sel( f.vi, t.vi, (vec_uint4)c.vi );
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
    friend inline v4float fma(  const v4float &a, const v4float &b, const v4float &c );
    friend inline v4float fms(  const v4float &a, const v4float &b, const v4float &c );
    friend inline v4float fnms( const v4float &a, const v4float &b, const v4float &c );
    friend inline v4float clear_bits(  const v4int &m, const v4float &a );
    friend inline v4float set_bits(    const v4int &m, const v4float &a );
    friend inline v4float toggle_bits( const v4int &m, const v4float &a );
    friend inline void increment_4x1( float * ALIGNED p, const v4float &a );
    friend inline void decrement_4x1( float * ALIGNED p, const v4float &a );
    friend inline void scale_4x1(     float * ALIGNED p, const v4float &a );
    // FIXME: crack
    
  public:

    // v4float constructors / destructors
    
    v4float() {}                                  // Default constructor
    v4float( const v4float &a ) {                 // Copy constructor
      vf = a.vf;
    }
    v4float( const v4 &a ) {                      // Initialize from mixed
      vf = a.vf;
    }
    v4float( const float &a ) {                   // Initialize from scalar
      vf = spu_splats(a);
    }
    v4float( const float &f0, const float &f1,
             const float &f2, const float &f3 ) { // Initalize from scalars
      // FIXME: vec_float4 t = { f0, f1, f2, f3 }; ... seg faults gcc .. MAYBE NOT ANYMORE
      f[0] = f0; f[1] = f1; f[2] = f2; f[3] = f3;
    }
    ~v4float() {}                                 // Destructor

    // v4float assignment operators

#   define ASSIGN(op,intrinsic)                         \
    inline v4float &operator op( const v4float &b ) {	\
      vf = intrinsic( vf, b.vf );                       \
      return *this;                                     \
    }

    inline v4float &operator =( const v4float &b ) {
      vf = b.vf;
      return *this;
    }

    ASSIGN(+=,spu_add)
    ASSIGN(-=,spu_sub)
    ASSIGN(*=,spu_mul)

    inline v4float &operator /=( const v4float &a ) {
      vec_float4 a_vf = a.vf, b_vf;

      // Compute an estimate of the reciprocal of a (12-bit accurate)

      b_vf = spu_re( a_vf );
     
      // FIXME: CHECK NUMERICS ... MAY WANT ADDITIONAL N-R STEPS:
      // b_vf = spu_nmsub(a_vf,spu_mul(b_vf,b_vf),spu_add(b_vf,b_vf));

      // Compute a * refined( (1/b)_estimate ) to get result a/b

      vf = spu_mul( vf,
                    spu_nmsub(a_vf,spu_mul(b_vf,b_vf),spu_add(b_vf,b_vf) ) );

      return *this;
    }

#   undef ASSIGN

    // v4float member access operator

    float &operator()(const int n) { return f[n]; }

  };

  // v4float prefix unary operators

  inline v4float operator +( const v4float &a ) {
    v4float b;
    b.vf = a.vf;
    return b;
  }

  inline v4float operator -( const v4float &a ) {
    v4float b;
    b.vf = spu_sub( _vfzero, a.vf );
    return b;
  }

  inline v4int operator !( const v4float &a ) {
    v4int b;
    b.vi = (vec_int4)spu_cmpeq( a.vf, _vfzero );
    return b;
  }

  // v4float prefix increment / decrement operators

  inline v4float operator ++( v4float &a ) {
    vec_float4 a_vf = a.vf;
    v4float b;
    a_vf = spu_add( a_vf, _vfone );
    a.vf = a_vf;
    b.vf = a_vf;
    return b;
  }

  inline v4float operator --( v4float &a ) {
    vec_float4 a_vf = a.vf;
    v4float b;
    a_vf = spu_sub( a_vf, _vfone );
    a.vf = a_vf;
    b.vf = a_vf;
    return b;
  }

  // v4float postfix increment / decrement operators

  inline v4float operator ++( v4float &a, int ) {
    vec_float4 a_vf = a.vf;
    v4float b;
    b.vf = a_vf;
    a_vf = spu_add( a_vf, _vfone );
    a.vf = a_vf;
    return b;
  }

  inline v4float operator --( v4float &a, int ) {
    vec_float4 a_vf = a.vf;
    v4float b;
    b.vf = a_vf;
    a_vf = spu_sub( a_vf, _vfone );
    a.vf = a_vf;
    return b;
  }

  // v4float binary operators
    
# define BINARY(op,intrinsic)                                        \
  inline v4float operator op( const v4float &a, const v4float &b ) { \
    v4float c;                                                       \
    c.vf = intrinsic( a.vf, b.vf );                                  \
    return c;                                                        \
  }

  BINARY(+,spu_add)
  BINARY(-,spu_sub)
  BINARY(*,spu_mul)

  inline v4float operator /( const v4float &n, const v4float &a ) {
    vec_float4 a_vf = a.vf, b_vf;
    v4float c;
    
    // Compute an estimate of the reciprocal of a (12-bit accurate)
    
    b_vf = spu_re( a_vf );
    
    // FIXME: CHECK NUMERICS ... MAY WANT ADDITIONAL N-R STEPS:
    // b_vf = spu_nmsub(a_vf,spu_mul(b_vf,b_vf),spu_add(b_vf,b_vf));
    
    // Compute a * refined( (1/b)_estimate ) to get result a/b
    
    c.vf = spu_mul( n.vf,
                    spu_nmsub(a_vf,spu_mul(b_vf,b_vf),spu_add(b_vf,b_vf) ) );

    return c;
  }

# undef BINARY

  // v4float logical operators

# define LOGICAL(op,assembly)                                      \
  inline v4int operator op( const v4float &a, const v4float &b ) { \
    v4int c;                                                       \
    c.vi = assembly;                                               \
    return c;                                                      \
  }

  LOGICAL(< ,(vec_int4)spu_cmpgt(b.vf,a.vf))
  LOGICAL(> ,(vec_int4)spu_cmpgt(a.vf,b.vf))
  LOGICAL(==,(vec_int4)spu_cmpeq(a.vf,b.vf))
  LOGICAL(!=,(vec_int4)spu_xor(_vutrue,spu_cmpeq(a.vf,b.vf))) // FIXME: ALMOST CERTAINLY NOT IEEE.  MAY BE A BETTER WAY.
  LOGICAL(<=,(vec_int4)spu_xor(_vutrue,spu_cmpgt(a.vf,b.vf))) // FIXME: SEE ABOVE
  LOGICAL(>=,(vec_int4)spu_xor(_vutrue,spu_cmpgt(b.vf,a.vf))) // FIXME: SEE ABOVE
  LOGICAL(&&,(vec_int4)spu_nor(  spu_cmpeq(a.vf,_vfzero),
                                 spu_cmpeq(b.vf,_vfzero) ))
  LOGICAL(||,(vec_int4)spu_nand( spu_cmpeq(a.vf,_vfzero),
                                 spu_cmpeq(b.vf,_vfzero) ))

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
    b.vf = (vec_float4)spu_andc( (vec_uint4)a.vf, (vec_uint4)_visignbits );
    return b;            
  }

  inline v4float sqrt( const v4float &a ) {
    vec_float4 a_vf = a.vf, b_vf;
    v4float b;

    // Compute an estimate of the rsqrt (12-bit accurate)

    b_vf = spu_rsqrte( a_vf );

    // FIXME: CHECK NUMERICS.  MIGHT WANT TO USE ADDITIONAL REFINEMENT:
    // b_vf = spu_madd( _vfhalf,
    //                  spu_nmsub( a_vf,
    //                             spu_mul( b_vf, spu_mul( b_vf, b_vf ) ),
    //                             b_vf ),
    //                  b_vf );

    // Compute the sqrt(a) via a*refined_rsqrt_estimate(a) ~ sqrt(a)

    b.vf = spu_mul( a_vf,
                    spu_madd( _vfhalf,
                              spu_nmsub( a_vf,
                                         spu_mul( b_vf, spu_mul( b_vf, b_vf ) ),
                                         b_vf ),
                              b_vf ) );

    return b;
  }

  inline v4float copysign( const v4float &a, const v4float &b ) {
    v4float c;
    c.vf = (vec_float4)spu_or( spu_andc( (vec_uint4)a.vf, (vec_uint4)_visignbits ),
                               spu_and(  (vec_uint4)b.vf, (vec_uint4)_visignbits ) );
    return c;
  }

# undef CMATH_FR1
# undef CMATH_FR2

  // v4float miscelleanous functions
  
  inline v4float rsqrt_approx( const v4float &a ) {
    vec_float4 a_vf = a.vf;
    v4float b;
    b.vf = spu_rsqrte( a_vf );
    return b;
  }
  
  inline v4float rsqrt( const v4float &a ) {
    vec_float4 a_vf = a.vf, b_vf;
    v4float b;

    // Compute an estimate of the rsqrt (12-bit accurate)

    b_vf = spu_rsqrte( a_vf );

    // Refine the estimate with N-R.
    // FIXME: CHECK NUMERICS.  MIGHT WANT TO USE ADDITIONAL REFINEMENT.

    b.vf = spu_madd( _vfhalf,
                     spu_nmsub( a_vf,
                                spu_mul( b_vf, spu_mul( b_vf, b_vf ) ),
                                b_vf ),
                     b_vf );
    
    return b;
  }

  inline v4float rcp_approx( const v4float &a ) {
    vec_float4 a_vf = a.vf;
    v4float b;
    b.vf = spu_re( a_vf );
    return b;
  }
  
  inline v4float rcp( const v4float &a ) {
    vec_float4 a_vf = a.vf, b_vf;
    v4float b;

    // Compute an estimate of the reciprocal of a (12-bit accurate)

    b_vf = spu_re( a_vf );
    
    // Perform Newton-Raphson refinement (one step should give 24-bit)
    // FIXME: CHECK NUMERICS ... MAY WANT TO USE ADDITIONAL REFINEMNT.
    
    b.vf = spu_nmsub(a_vf,spu_mul(b_vf,b_vf),spu_add(b_vf,b_vf));

    return b;
  }

  inline v4float fma(  const v4float &a, const v4float &b, const v4float &c ) {
    v4float d;
    d.vf = spu_madd( a.vf, b.vf, c.vf );
    return d;
  }

  inline v4float fms(  const v4float &a, const v4float &b, const v4float &c ) {
    v4float d;
    d.vf = spu_msub( a.vf, b.vf, c.vf );
    return d;
  }

  inline v4float fnms( const v4float &a, const v4float &b, const v4float &c ) {
    v4float d;
    d.vf = spu_nmsub( a.vf, b.vf, c.vf );
    return d;
  }

  inline v4float clear_bits( const v4int &m, const v4float &a ) {
    v4float b;
    b.vf = (vec_float4)spu_andc( (vec_uint4)a.vf, (vec_uint4)m.vi );
    return b;
  }

  inline v4float set_bits( const v4int &m, const v4float &a ) {
    v4float b;
    b.vf = (vec_float4)spu_or( (vec_uint4)a.vf, (vec_uint4)m.vi );
    return b;
  }

  inline v4float toggle_bits( const v4int &m, const v4float &a ) {
    v4float b;
    b.vf = (vec_float4)spu_xor( (vec_uint4)a.vf, (vec_uint4)m.vi );
    return b;
  }

  inline void increment_4x1( float * ALIGNED p, const v4float &a ) {
    *((vec_float4 * ALIGNED)p) = spu_add( *((vec_float4 * ALIGNED)p), a.vf );
  }

  inline void decrement_4x1( float * ALIGNED p, const v4float &a ) {
    *((vec_float4 * ALIGNED)p) = spu_sub( *((vec_float4 * ALIGNED)p), a.vf );
  }

  inline void scale_4x1( float * ALIGNED p, const v4float &a ) {
    *((vec_float4 * ALIGNED)p) = spu_mul( *((vec_float4 * ALIGNED)p), a.vf );
  }

} // namespace v4

#endif /* _v4_spu_hxx_ */
