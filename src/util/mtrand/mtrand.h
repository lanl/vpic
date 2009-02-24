/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Adapted from earlier V4PIC versions.
 *
 */

#ifndef _mtrand_h_
#define _mtrand_h_

#include "util_base.h"

// Setup the mt_rng opaque handle

struct mt_rng;
typedef struct mt_rng mt_rng_t;

// Note: All these trinary operators should be evaluated at compile
// time for any self-respecting compiler.

#define mt_CRAND_MAX  ((CHAR_BIT*sizeof(signed char)       <32) ? SCHAR_MAX : 0x7fffffff  )
#define mt_HRAND_MAX  ((CHAR_BIT*sizeof(short int)         <32) ? SHRT_MAX  : 0x7fffffff  )
#define mt_RAND_MAX   ((CHAR_BIT*sizeof(int)               <32) ? INT_MAX   : 0x7fffffff  )
#define mt_LRAND_MAX  ((CHAR_BIT*sizeof(long int)          <32) ? LONG_MAX  : 0x7fffffffL )
#define mt_UCRAND_MAX ((CHAR_BIT*sizeof(unsigned char)     <32) ? UCHAR_MAX : 0xffffffffU )
#define mt_UHRAND_MAX ((CHAR_BIT*sizeof(unsigned short int)<32) ? USHRT_MAX : 0xffffffffU )
#define mt_URAND_MAX  ((CHAR_BIT*sizeof(unsigned int)      <32) ? UINT_MAX  : 0xffffffffU )
#define mt_ULRAND_MAX ((CHAR_BIT*sizeof(unsigned long int) <32) ? ULONG_MAX : 0xffffffffUL)

BEGIN_C_DECLS

/********************************
 * Constructors and destructors *
 ********************************/

mt_rng_t *
new_mt_rng( unsigned int seed );

void
delete_mt_rng( mt_rng_t * rng );
 
/*********************************
 * Seeders and related functions *
 *********************************/

void
seed_mt_rng( mt_rng_t * rng,
             unsigned int seed );

size_t
get_mt_rng_size( mt_rng_t * rng );

void
get_mt_rng_state( mt_rng_t * rng,
                  void * state );

void
set_mt_rng_state( mt_rng_t * rng,
                  const void * state,
                  size_t sz );

/**********************
 * Integer generators *
 **********************/

#define int_gen_decl( name, type )                                      \
  signed   type mt_##name ( mt_rng_t * rng );                           \
  unsigned type mt_u##name( mt_rng_t * rng );                           \
  void mt_##name##_fill   ( mt_rng_t * rng, signed   type * x, size_t n ); \
  void mt_u##name##_fill  ( mt_rng_t * rng, unsigned type * x, size_t n )

int_gen_decl( crand, char      );
int_gen_decl( hrand, short int );
int_gen_decl( rand,  int       );
int_gen_decl( lrand, long int  );

#undef int_gen_decl

/*************************************
 * Uniform floating point generators *
 *************************************/

#define float_gen_decl( name, type )                                   \
  type mt_##name          ( mt_rng_t * rng );                          \
  type mt_##name##_c0     ( mt_rng_t * rng );                          \
  type mt_##name##_c1     ( mt_rng_t * rng );                          \
  type mt_##name##_c      ( mt_rng_t * rng );                          \
  void mt_##name##_fill   ( mt_rng_t * rng, type * x, size_t n );      \
  void mt_##name##_c0_fill( mt_rng_t * rng, type * x, size_t n );      \
  void mt_##name##_c1_fill( mt_rng_t * rng, type * x, size_t n );      \
  void mt_##name##_c_fill ( mt_rng_t * rng, type * x, size_t n )

float_gen_decl( frand,      float  );
float_gen_decl( fast_drand, double );
float_gen_decl( drand,      double );

#undef float_gen_decl

/*************************
 * Speciality generators *
 *************************/

double
mt_drandn( mt_rng_t * rng );

void
mt_drandn_fill( mt_rng_t * rng,
                double * x,
                size_t n );

float
mt_frandn( mt_rng_t * rng );

void
mt_frandn_fill( mt_rng_t * rng,
                float * x,
                size_t n );

double
mt_drande( mt_rng_t * rng );

void
mt_drande_fill( mt_rng_t * rng,
                double * x,
                size_t n );

void
mt_randperm( mt_rng_t * rng,
             int * x,
             int n );

#define MT_SHUFFLE(rng,x,n) mt_shuffle((rng),(x),(n),sizeof(*(x)))

void
mt_shuffle( mt_rng_t * rng,
            void * x,
            size_t n,
            size_t sz );

END_C_DECLS

#endif // _mtrand_h_

