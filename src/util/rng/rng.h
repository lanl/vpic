#ifndef _rng_h_
#define _rng_h_

#include "../util_base.h"

/* rng_t opaque handle */

struct rng;
typedef struct rng rng_t;

/* A rng_pool is a collection of random number generators. */

typedef struct rng_pool {
  rng_t ** rng; /* Random number generators (indexed 0:n_rng-1) */
  int n_rng;    /* Number of random number generators in pool */
} rng_pool_t;

#if !defined(__SPU__) /* The SPU version will define its equivalents */

BEGIN_C_DECLS

/* In rng_pool.c. */

rng_pool_t *              /* New pool (already seeded via seed_rng_pool) */
new_rng_pool( int n_rng,  /* Number of generators for pool */
              int seed,   /* Pool seed (different meaning from seed_rng) */
              int sync ); /* True for synchronized seeding */

void
delete_rng_pool( rng_pool_t * RESTRICT rp ); /* Pool to delete */

/* In seed_rng_pool, seeding is done such that:
     local_pool = seed_rng_pool( rp, seed, 0 );   
     sync_pool  = seed_rng_pool( rp, seed, 1 );
   gives each local_pool rng and each sync_pool rng has a unique seed
   on all calling processes and that the sync pool rngs are
   identically initialized on all calling processes. */

/* FIXME: WE NEED BIGGER SEEDS.  NOTE THAT THE EFFECT SEED SPACE FOR
   POOLS IS ROUGHLY FLOOR( UINT_MAX / (n_rng*(world_size+1)) )  */

rng_pool_t *                             /* Returns rp */
seed_rng_pool( rng_pool_t * RESTRICT rp, /* Pool to seed */
               int seed,                 /* Pool seed (different meaning from
                                            seed_rng) */
               int sync );               /* True for synchronized seeding */

/* In rng.c */

rng_t *              /* New generator (already seeded via seed_rng) */
new_rng( int seed ); /* Random number generator seed */

void
delete_rng( rng_t * RESTRICT r ); /* Generator to delete */

/* FIXME: WE NEED TO BIGGER SEEDS.  SEE ABOVE. */

rng_t *                            /* Returns r */
seed_rng( rng_t * RESTRICT r,      /* Generator to seed */
          int              seed ); /* Seed */

/* Integer random generators make uniform rands on [0,INTTYPE_MAX] for
   signed types and on [0,UINTTYPE_MAX] for unsigned types.  There are
   singleton generators for each primitive integral type (including
   both signed and unsigned).  Each singleton generator has a
   corresponding mass production generator.  The singleton generators
   are not error trapped to reduce their overhead as much as possible.
   Generators are named as follows:

     [c,h,i,l,i8,i16,i32,i64,uc,uh,ui,ul,u8,u16,u32,u64]rand[,_fill]

   where the generated data type is:
      c   => char,    uc  => unsigned char, 
      h   => short,   uh  => unsigned short,
      i   => int,     ui  => unsigned int,
      l   => long,    ul  => unsigned long
      i8  => int8_t,  u8  => uint8_t,
      i16 => int16_t, u16 => uint16_t,
      i32 => int32_t, u32 => uint32_t,
      i64 => int64_t, u64 => uint64_t
   and _fill indicates a mass production variant */

#define _( prefix, type )                                        \
type                                /* Returns sample deviate */ \
prefix##rand( rng_t * RESTRICT r ); /* Generator to use */       \
                                                                 \
type *                                 /* Returns x */           \
prefix##rand_fill( rng_t * RESTRICT r, /* Generator to use */    \
                   type  * RESTRICT x, /* Array to fill */       \
                   size_t str_ele,     /* Element stride */      \
                   size_t n_ele );     /* Number of elements */

_( c, char  ) _( uc, unsigned char  )
_( h, short ) _( uh, unsigned short )
_( i, int   ) _( ui, unsigned int   )
_( l, long  ) _( ul, unsigned long  )

_(  i8,  int8_t ) _(  u8,  uint8_t )
_( i16, int16_t ) _( u16, uint16_t )
_( i32, int32_t ) _( u32, uint32_t )
_( i64, int64_t ) _( u64, uint64_t )

#undef _

/* Floating point random generators make uniform rands on [0,1],
   [0,1), (0,1] or (0,1), depending on the variant selected.  Each
   variant has a rigorous interpretation.  For the open variants, if
   the interval [0,1) is divided into 2^N sub-intervals (where N=23
   and N=52 for single and double precision respectively), the value
   returned is as if a true uniform rand on [0,1) was produced and
   rounded to the midpoint of each sub-interval (with ties breaking to
   the even interval).  For the other variants, if the interval [0,1]
   is divided into 2^N intervals (where N=24 and N=53 for single and
   double precision respectively ... note these variants have an extra
   bit of precision), the value returned is as if a true uniform rand
   on [0,1) was produced and rounded to the nearest even / down / up
   location between each sub-intervals for the closed, half open at 1
   and half open at 0 variants respectively.

   Note then that:
   - In the open variant, 0 or 1 can never be returned.
   - In the closed variant, 0 or 1 can both be returned.
   - In the half open at 1 variant, 1 can never be returned.
   - In the half open at 0 variant, 0 can never be returned.
  
   There are single generators for each primitive floating point type
   and domain.  Each singleton generator has a corresponding mass
   production generator.  The singleton generators not error trapped
   to reduce their overhead as much as possible.  Generators are named
   as follows:

      [f,d]rand[,_c0,_c1,_c][,_fill]

   where the generated data type is:
      f => float, d => double
   the interval is:
      nothing => open
      _c0     => [0,1)
      _c1     => (0,1]
      _c      => [0,1]
   and _fill indicates a mass production variant */

#define _( type, prefix, variant )                                         \
type                                         /* Returns sample deviate */  \
prefix##rand##variant( rng_t * RESTRICT r ); /* Generator to use */        \
                                                                           \
type *                                            /* Returns x */          \
prefix##rand##variant##_fill( rng_t * RESTRICT r, /* Generator to use */   \
                              type  * RESTRICT x, /* Array to fill */      \
                              size_t str_ele,     /* Element stride */     \
                              size_t n_ele );     /* Number of elements */

_( float, f,     ) _( double, d,     )
_( float, f, _c0 ) _( double, d, _c0 )
_( float, f, _c1 ) _( double, d, _c1 )
_( float, f, _c  ) _( double, d, _c  )

#undef _

/* The normal generators generate a normally distributed random number
   (f(x) = exp( -x^2 / 2 ) / sqrt( 2*pi ) for x in (-inf,inf)).  Based
   on the Ziggurat method under the hood. */

float                         /* Returns sample deviate */
frandn( rng_t * RESTRICT r ); /* Generator to use */

float *                                /* Returns x */
frandn_fill( rng_t * RESTRICT r,       /* Generator to use */
             float * RESTRICT x,       /* Array to fill */
             size_t           str_ele, /* Element stride */
             size_t           n_ele ); /* Number of elements */

double                        /* Returns sample deviate */
drandn( rng_t * RESTRICT r ); /* Generator to use */

double *                                 /* Returns x */
drandn_fill( rng_t  * RESTRICT r,        /* Generator to use */
             double * RESTRICT x,        /* Array to fill */
             size_t            str_ele,  /* Element stride */
             size_t            n_ele );  /* Number of elements */

/* The exponential generators generate an exponentially distributed
   random number (f(x) = exp(-x) for x in [0,inf).  Based on the
   transformation method under the hood. */
   
float                         /* Returns sample deviate */
frande( rng_t * RESTRICT r ); /* Generator to use */

float *                                /* Returns x */
frande_fill( rng_t * RESTRICT r,       /* Generator to use */
             float * RESTRICT x,       /* Array to fill */
             size_t           str_ele, /* Element stride */
             size_t           n_ele ); /* Number of elements */

double                        /* Returns sample deviate */
drande( rng_t * RESTRICT r ); /* Generator to use */

double *                                /* Returns x */
drande_fill( rng_t  * RESTRICT r,       /* Generator to use */
             double * RESTRICT x,       /* Array to fill */
             size_t            str_ele, /* Element stride */
             size_t            n_ele ); /* Number of elements */

/* Specialty generators */

int *                           /* Returns x */
randperm( rng_t * RESTRICT r,   /* Generator to use */
          int   * RESTRICT x,   /* 0:n-1 indexed, holds a random
                                   permutation of 0:n-1 on output */
          int              n ); /* Permutation size */

/* This function is most efficient when str_ele is a integer multiple
   of sz_ele and sz_ele is either 0, 1, 2, 4 or 8. */
void *                             /* Returns x */
shuffle( rng_t * RESTRICT r,       /* Generator to use */
         void  * RESTRICT x,       /* Elements to shuffle */
         size_t           sz_ele,  /* Element _byte_ size */
         size_t           str_ele, /* Element _byte_ stride */
         size_t           n_ele ); /* Number of elements */

END_C_DECLS

#endif /* __SPU__ */
#endif /* _rng_h_ */
