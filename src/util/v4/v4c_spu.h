#ifndef _v4c_spu_h_
#define _v4c_spu_h_
#ifdef CELL_SPU_BUILD

/* Since XLC does not seem able to handle C++ of V4 for the SPU and
   generates terrible assembly, this file contains a "baby" version
   of V4 that is compatible with C such that XLC can compile it and
   not generate terrible assembly. */

/* FIXME: THIS HEADER IS NOT PRODUCTION QUALITY YET!! */

#include <spu_intrinsics.h>

/* All macros below are robust unless otherwise noted. */

#define USING_V4C                                               \
   const vec_uchar16 _unpackl = {  0, 1, 2, 3,   16,17,18,19,   \
                                   4, 5, 6, 7,   20,21,22,23 }; \
   const vec_uchar16 _unpackh = {  8, 9,10,11,   24,25,26,27,   \
                                  12,13,14,15,   28,29,30,31 }; \
   vec_float4 _a, _b, _c, _d

/* FIXME: THIS MACRO IS NOT ROBUST! 
   MULTIPLE SEMANTIC STATEMENTS
   MUST DECLARE _a, _b, _c, _d, unpackl = _unpackl, unpackh = _unpackh TO USE
   a, b, c and d ARE REFERENCED THREE TIMES */

#define TRANSPOSE( a, b, c, d )                                  \
  /* Step 1: Interleave top and bottom half */                   \
  _a  = spu_shuffle( (a), (c), _unpackl ); /* a =  0  8  1  9 */ \
  _b  = spu_shuffle( (b), (d), _unpackl ); /* b =  4 12  5 13 */ \
  _c  = spu_shuffle( (a), (c), _unpackh ); /* c =  2 10  3 11 */ \
  _d  = spu_shuffle( (b), (d), _unpackh ); /* d =  6 14  7 15 */ \
  /* Step 2: Interleave even and odd rows */                     \
  (a) = spu_shuffle( _a,  _b,  _unpackl ); /* a =  0  4  8 12 */ \
  (b) = spu_shuffle( _a,  _b,  _unpackh ); /* b =  1  5  9 13 */ \
  (c) = spu_shuffle( _c,  _d,  _unpackl ); /* c =  2  6 10 14 */ \
  (d) = spu_shuffle( _c,  _d,  _unpackh )  /* d =  3  7 11 15 */

/* FIXME: THIS MACRO IS NOT ROBUST!
   MULTIPLE SEMANTIC STATEMENTS
   MUST DECLARE _a, _b, unpackl = _unpackl, unpackh = _unpackh TO USE
   a AND b ARE REFERENCED TWICE */

#define HALF_TRANSPOSE( a, b, c, d )                             \
  /* Step 1: Interleave top and bottom half */                   \
  _a  = spu_shuffle( (a), (c), _unpackl ); /* a =  0  8  1  9 */ \
  _b  = spu_shuffle( (b), (d), _unpackl ); /* b =  4 12  5 13 */ \
  /* Step 2: Interleave even and odd rows */                     \
  (a) = spu_shuffle( _a,  _b,  _unpackl ); /* a =  0  4  8 12 */ \
  (b) = spu_shuffle( _a,  _b,  _unpackh )  /* b =  1  5  9 13 */ 

/* FIXME: THIS MACRO MAY HAVE A POTENTIAL STRICT ALIASING ISSUE (PROBABLY
   UNAVOIDABLE THE WAY SIMD WORKS ON GCC AND XLC) */
#define LOAD_4x1( p, v )       (v) = *((const vec_float4 *)(p))
#define LOAD_INT_4x1( p, v )   (v) = *((const vec_int4   *)(p))
#define LOAD_UINT_4x1( p, v )  (v) = *((const vec_uint4  *)(p))

/* FIXME: THIS MACRO MAY HAVE A POTENTIAL STRICT ALIASING ISSUE (PROBABLY
   UNAVOIDABLE THE WAY SIMD WORKS ON GCC AND XLC) */
#define STORE_4x1( v, p )     *((vec_float4 *)(p)) = (v)

/* FIXME: THIS MACRO IS TECHNICALLY NOT ROBUST! p IS REFERENCED TWICE! 
   THIS MACRO MAY HAVE A POTENTIAL STRICT ALISING ISSUE (PROBABLY
   UNAVOIDABLE THE WAY SIMD WORKS ON GCC AND XLC) */
#define INCREMENT_4x1( p, v ) \
  *((vec_float4 *)(p)) = spu_add( *((vec_float4 *)(p)), v )

/* NOTE: gcc allows for example a*b to automagically become the right kind
   of spu_mul but xlc does not.  Sigh ... */

#define VEC_FLOAT4( a )  spu_splats( (float)(a) )
#define VEC_UINT4( a )   spu_splats( (unsigned int)(a) )
#define VEC_INT4( a )    spu_splats( (int)(a) )

#define EXTRACT( v, n )  spu_extract( (v), (n) )
#define GATHER( a )      spu_extract( spu_gather( (a) ), 0 )

#define ADD( a, b )      spu_add(   (a), (b) )
#define SUB( a, b )      spu_sub(   (a), (b) )
#define MUL( a, b )      spu_mul(   (a), (b) )

#define FMA(  a, b, c )  spu_madd(  (a), (b), (c) )
#define FMS(  a, b, c )  spu_msub(  (a), (b), (c) )
#define FNMS( a, b, c )  spu_nmsub( (a), (b), (c) )

#define MERGE( c, t, f ) spu_sel(   (f), (t), (c) )
#define CZERO( c, a )    spu_andc(  (a), (vec_float4)(c) ) /* USE typeof?? */
#define NOTCZERO( c, a ) spu_and(   (a), (vec_float4)(c) ) /* USE typeof?? */

#define CMPGT( a, b )    spu_cmpgt( (a), (b) )
#define CMPLT( a, b )    spu_cmpgt( (b), (a) )
#define CMPEQ( a, b )    spu_cmpeq( (b), (a) )

#define AND( a, b )      spu_and( (a), (b) )
#define OR( a, b )       spu_or(  (a), (b) )
#define XOR( a, b )      spu_xor( (a), (b) )
#define NOT( a )         spu_xor( spu_splats( 0xffffffff ), (a) )

/* Compute an estimate of the rsqrt (12-bit accurate) */
/* Perform Newton-Raphson refinement (one step should give 24-bit) */
/* FIXME: CHECK NUMERICS.  MIGHT WANT TO USE ADDITIONAL REFINEMENT.
    _b = spu_madd( spu_nmsub( spu_mul( _b, _b ), _a, one ), \
                   spu_mul( _b, one_half ),                 \
                   _b ),                                    \
*/
/* FIXME: THIS MACRO IS NOT ROBUST!
   MUST DECLARE vec_float4 _a, _b and one_half = 0.5 to use.  */

#define RSQRT( a )                                          \
  ( _a = (a),                                               \
    _b = spu_rsqrte( _a ),                                  \
    spu_madd( spu_nmsub( spu_mul( _b, _b ), _a, one ),      \
              spu_mul( _b, one_half ),                      \
              _b ) )

/* Compute an estimate of the reciprocal of a (12-bit accurate) */
/* Perform Newton-Raphson refinement (one step should give 24-bit) */
/* FIXME: CHECK NUMERICS ... MAY WANT TO USE ADDITIONAL REFINEMNT.
    _b = spu_madd( spu_nmsub( _b, _a, one ), _b, _b ),    \
*/

/* FIXME: THIS MACRO IS NOT ROBUST! 
   MUST DECLARE vec_float4 _a, _b and one = 1 to use. */
#define RCP( a )                                          \
  ( _a = (a),                                             \
    _b = spu_re( _a ),                                    \
    spu_madd( spu_nmsub( _b, _a, one ), _b, _b ) )

#endif // CELL_SPU_BUILD 
#endif // _v4c_spu_h_ 

