/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Adapted from earlier V4PIC versions.
 *
 */

#include <stdlib.h> /* For malloc, NULL and size_t                     */
#include <math.h>   /* For sqrt, log and tan                           */
#include "mtrand_conv.h" /* For drand53_o, drand53_c0, drand53_c1, ... */
#include "mtrand.h" /* Assure consistency with the header              */

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029
#endif

#ifndef INT32_TYPE
#define INT32_TYPE int
#endif

#ifndef RESTRICT
#define RESTRICT
#endif

/*******************************************************
 * Mersenne-Twister 19937 Random Number Generator Core *
 *******************************************************/

/* The mt_uint32 datatype should be a 32-bit unsigned integer type. The 
   library should still work if it is even larger but it is not tested. */
typedef unsigned INT32_TYPE mt_uint32;

#define MT19937_N          624
#define MT19937_M          397
#define MT19937_TWIST(u,v) (((((u)&0x80000000UL)|((v)&0x7fffffffUL))>>1)^((-((v)&1))&0x9908b0dfUL))
#define MT19937_TEMPER(u)  (u) ^= ( (u) >> 11 );                \
                           (u) ^= ( (u) << 7  ) & 0x9d2c5680UL; \
                           (u) ^= ( (u) << 15 ) & 0xefc60000UL; \
                           (u) ^= ( (u) >> 18 )

typedef struct _mt_gen_t {
  mt_uint32 next;
  mt_uint32 state[MT19937_N];
} mt_gen_t;

static void mt19937_next_state( mt_gen_t *g ) {
  int j;
  mt_uint32 *p;

  g->next = 0;
  p = g->state;
  for( j=MT19937_N-MT19937_M+1; --j; p++ ) p[0] = p[MT19937_M]           ^ MT19937_TWIST( p[0], p[1] );
  for( j=MT19937_M            ; --j; p++ ) p[0] = p[MT19937_M-MT19937_N] ^ MT19937_TWIST( p[0], p[1] );
  p[0] = p[MT19937_M-MT19937_N] ^ MT19937_TWIST( p[0], g->state[0] );
}

#define urand32(g,y)					\
  if( (g)->next==MT19937_N ) mt19937_next_state(g);	\
  (y) = (g)->state[ (g)->next++ ];			\
  MT19937_TEMPER(y)

/*****************************************************************************
 * Constructors and destructors                                              *
 *****************************************************************************/

mt_handle mt_new_generator( unsigned int s ) {
  mt_handle h = (mt_handle)malloc(sizeof(mt_gen_t));
  if( h!=NULL ) mt_srand( h, s );
  return h;
}

void mt_delete_generator( mt_handle *h ) {
  if( h==NULL ) return;
  if( (*h)!=NULL ) free(*h);
  (*h) = NULL;
}

/*****************************************************************************
 * Seeders and related functions                                             *
 *****************************************************************************/

int mt_srand( mt_handle h, unsigned int s ) {
  mt_gen_t *g = (mt_gen_t *)h;
  int j;

  if( h==NULL ) return 1;
  g->next = MT19937_N;
  g->state[0] = s^0x900df00cUL;  /* state[0] on srand(1) is goodfood */
  g->state[0] &= 0xffffffffUL;   /* 64-bit compatibility */
  for( j=1; j<MT19937_N; j++ ) {
    g->state[j] = 1812433253UL*(g->state[j-1]^(g->state[j-1]>>30)) + j;
    g->state[j] &= 0xffffffffUL; /* 64-bit compatibility */
  }
  return 0;
}

#if CHAR_BIT==8
/* mt_getsize returns the number of chars required to hold the generator's
   internal state in the format used by mt_getstate function. This routine
   assumes 8-bit "char"s */
size_t mt_getsize( mt_handle h ) {
  return 4*(MT19937_N+1);
}

/* mt_getstate saves the state of the generator in a machine independent
   format on machines with 8-bits "char"s (true on virtually all hardware
   in the last 30 years) if the char array is large enough to hold it. */
int mt_getstate( mt_handle h, char *s, size_t n ) {
  mt_gen_t *g = (mt_gen_t *)h;
  size_t j;

  /* Error checking */
  if( h==NULL || s==NULL ) return -1;
  j = mt_getsize(h);
  if( n<j ) return j;

  /* Serialize g->next */
  *(s++) = (g->next & 0x000000ff) >> 0;
  *(s++) = (g->next & 0x0000ff00) >> 8;
  *(s++) = (g->next & 0x00ff0000) >> 16;
  *(s++) = (g->next & 0xff000000) >> 24;

  /* Serialize g->state */
  for( j=0; j<MT19937_N; j++ ) {
    *(s++) = (g->state[j] & 0x000000ff) >> 0;
    *(s++) = (g->state[j] & 0x0000ff00) >> 8;
    *(s++) = (g->state[j] & 0x00ff0000) >> 16;
    *(s++) = (g->state[j] & 0xff000000) >> 24;
  }

  return 0;
}

/* mt_setstate sets the state of the generator. It can take a state provided
   by mt_getstate or a user designed state of at least 5 bytes. A valid
   user designed state has at least one non-zero char among s[4] through
   s[min(n,mt_getsize(h))-1]. This routine assumes 8-bit chars. */
int mt_setstate( mt_handle h, const char *s, size_t n ) {
  mt_gen_t *g = (mt_gen_t *)h;
  size_t j, k;

  /* Error checking */
  if( h==NULL || s==NULL || n<5 ) return 1;
  j = 4*(MT19937_N+1); if( j>n ) j=n;
  for( k=4; k<j; k++ ) if( s[k]!=0 ) break;
  if( k==j ) return 2;

  /* Extract g->next */
  g->next  = ((mt_uint32)(s[0])) << 0;
  g->next |= ((mt_uint32)(s[1])) << 8;
  g->next |= ((mt_uint32)(s[2])) << 16;
  g->next |= ((mt_uint32)(s[3])) << 24;
  if( g->next<0 || g->next>MT19937_N ) g->next=MT19937_N;

  /* Extract g->state */
  k=4;
  for( j=0, k=4; j<MT19937_N; j++ ) {
    g->state[j]  = ((mt_uint32)(s[k++])) << 0;  if( k==n ) k=4;
    g->state[j] |= ((mt_uint32)(s[k++])) << 8;  if( k==n ) k=4;
    g->state[j] |= ((mt_uint32)(s[k++])) << 16; if( k==n ) k=4;
    g->state[j] |= ((mt_uint32)(s[k++])) << 24; if( k==n ) k=4;
  }

  return 0;
}
#endif

/*****************************************************************************
 * Generate integer random numbers                                           *
 *****************************************************************************/
/* The irf code-generator is very portable but not maximally efficient when
   dealing with short and char datatypes (for example, on a 32-bit system with
   8-bit chars, 4 chars could be made for every 32-bit rand). Note: any
   self-respecting compiler will optimize the bit shift at compile time */
#define irf( name, t, s )				                    \
t mt_##name( mt_handle h ) {						    \
  mt_uint32 y;								    \
  urand32((mt_gen_t *)h,y);						    \
  return (t)(y>>(s+((CHAR_BIT*sizeof(t)<32)?(32-CHAR_BIT*sizeof(t)):0)));   \
}									    \
int mt_##name##_fill( mt_handle h, t * RESTRICT x, size_t n ) {		    \
  mt_uint32 y;								    \
  if( h==NULL || x==NULL || n<1 ) return 1;				    \
  for(;n;n--) {								    \
    urand32((mt_gen_t *)h,y);						    \
    *(x++)=(t)(y>>(s+((CHAR_BIT*sizeof(t)<32)?(32-CHAR_BIT*sizeof(t)):0))); \
  }									    \
  return 0;                                                                 \
}                    
irf(crand,   signed char,        1) /* Force signed chars */
irf(shrand,  signed short int,   1) /* Can't use 's' due to srand */
irf(rand,    signed int,         1)
irf(lrand,   signed long int,    1)
irf(ucrand,  unsigned char,      0)
irf(ushrand, unsigned short int, 0)
irf(urand,   unsigned int,       0)
irf(ulrand,  unsigned long int,  0)
#undef irf

/*****************************************************************************
 * Generate floating point random numbers                                    *
 *****************************************************************************/
#define frf( name, type, which )				\
type mt_##name( mt_handle h ) {					\
  mt_uint32 a;							\
  urand32((mt_gen_t *)h,a);					\
  return which(a);						\
}								\
int mt_##name##_fill( mt_handle h, type * RESTRICT x, size_t n ) {	\
  mt_uint32 a;							\
  if( h==NULL || x==NULL || n<1 ) return 1;			\
  for(;n;n--) {							\
    urand32((mt_gen_t *)h,a);					\
    *(x++) = which(a);						\
  }								\
  return 0;                                                     \
}
#define drf( name, which )					\
double mt_##name( mt_handle h ) {				\
  mt_uint32 a, b;						\
  urand32((mt_gen_t *)h,a);					\
  urand32((mt_gen_t *)h,b);					\
  return which(a,b);						\
}								\
int mt_##name##_fill( mt_handle h, double * RESTRICT x, size_t n ) {	\
  mt_uint32 a, b;						\
  if( h==NULL || x==NULL || n<1 ) return 1;			\
  for(;n;n--) {							\
    urand32((mt_gen_t *)h,a);					\
    urand32((mt_gen_t *)h,b);					\
    *(x++) = which(a,b);					\
  }								\
  return 0;                                                     \
}
frf(frand,         float,  frand24_o )
frf(frand_c0,      float,  frand24_c0)
frf(frand_c1,      float,  frand24_c1)
frf(frand_c,       float,  frand24_c )
frf(fast_drand,    double, drand32_o )
frf(fast_drand_c0, double, drand32_c0)
frf(fast_drand_c1, double, drand32_c1)
frf(fast_drand_c,  double, drand32_c )
drf(drand,    drand53_o )
drf(drand_c0, drand53_c0)
drf(drand_c1, drand53_c1)
drf(drand_c,  drand53_c )
#undef frf
#undef drf

/********************************************
 * Normal random number based distributions *
 ********************************************/

#define normal_drand_core()			\
  do {						\
    urand32((mt_gen_t *)h,a);			\
    urand32((mt_gen_t *)h,b);			\
    d1 = drand53_o(a,b); d1 += d1 - 1;		\
    urand32((mt_gen_t *)h,a);			\
    urand32((mt_gen_t *)h,b);			\
    d2 = drand53_o(a,b); d2 += d2 - 1;		\
    d3 = d1*d1 + d2*d2;				\
  } while( d3>1 || d3==0 );			\
  d3 = -log(d3)/d3; d3 += d3; d3 = sqrt(d3)

/*****************************************************************************
 * Generate a normal random number. Range is (-inf,inf)                      *
 * f(x) = exp( -x^2 / 2 ) / sqrt( 2*pi )                                     *
 *****************************************************************************/

double mt_normal_drand( mt_handle h ) {
  double d1, d2, d3;
  mt_uint32 a, b;
  normal_drand_core();
  return d1*d3; /* d2*d3 is also a normal drand but it is wasted! */
}

int mt_normal_drand_fill( mt_handle h, double * RESTRICT x, size_t n ) {
  double d1, d2, d3;
  mt_uint32 a, b;
  if( h==NULL || x==NULL || n<1 ) return 1;
  if( n&1 ) *(x++) = mt_normal_drand(h);
  n>>=1;
  for(;n;n--) {
    normal_drand_core();
    x[0] = d1*d3;
    x[1] = d2*d3;
    x +=2;
  }
  return 0;
}

/*****************************************************************************
 * Generate a log normal random number. Range is (0,inf)                     *
 * f(x) = exp( -(ln x)^2 / (2 sigma^2) ) / (x*sigma*sqrt( 2*pi ))            *
 *****************************************************************************/

double mt_lognormal_drand( mt_handle h, double sigma ) {
  double d1, d2, d3;
  mt_uint32 a, b;
  normal_drand_core();
  return exp(sigma*d1*d3); /* d2*d3 is wasted! */
}

int mt_lognormal_drand_fill( mt_handle h, double sigma,
                             double * RESTRICT x, size_t n ) {
  double d1, d2, d3;
  mt_uint32 a, b;
  if( h==NULL || x==NULL || n<1 ) return 1;
  if( n&1 ) *(x++) = mt_lognormal_drand(h,sigma);
  n>>=1;
  for(;n;n--) {
    normal_drand_core();
    d3 *= sigma;
    x[0] = exp(d1*d3);
    x[1] = exp(d2*d3);
    x +=2;
  }
  return 0;
}

/*****************************************************************************
 * Generate a Birnbaum-Saunders random number. Range is (0,inf)              *
 * f(x) = ?                                                                  *
 *****************************************************************************/

double mt_bs_drand( mt_handle h, double gamma ) {
  double d1, d2, d3;
  mt_uint32 a, b;
  normal_drand_core();
  d1 *= 0.5*gamma*d3; d1 += sqrt( 1+d1*d1 ); d1 *= d1;
  return d1; /* d2 is wasted! */
}

int mt_bs_drand_fill( mt_handle h, double gamma,
                      double * RESTRICT x, size_t n ) {
  double d1, d2, d3;
  mt_uint32 a, b;
  if( h==NULL || x==NULL || n<1 ) return 1;
  if( n&1 ) *(x++) = mt_bs_drand(h,gamma);
  n>>=1;
  for(;n;n--) {
    normal_drand_core();
    d3 *= 0.5;
    d1 *= gamma*d3; d1 += sqrt( 1+d1*d1 ); d1 *= d1;
    d2 *= gamma*d3; d2 += sqrt( 1+d2*d2 ); d2 *= d2;
    x[0] = d1;
    x[1] = d2;
    x +=2;
  }
  return 0;
}

#undef normal_drand_core

/*************************************************
 * Exponential random number based distributions *
 *************************************************/

/*****************************************************************************
 * Generate an exponential random number. Range is (0,inf)                   *
 * f(x) = exp(-x) for x>0, 0 otherwise                                       *
 *****************************************************************************/

double mt_exp_drand( mt_handle h ) {
  mt_uint32 a, b;
  urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b);
  return -log(drand53_o(a,b));
}

int mt_exp_drand_fill( mt_handle h, double * RESTRICT x, size_t n ) {
  mt_uint32 a, b;
  if( h==NULL || x==NULL || n<1 ) return 1;
  for(;n;n--) {
    urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b);
    *(x++) = -log(drand53_o(a,b));
  }
  return 0;
}

/*****************************************************************************
 * Generate a double exponential distributed random number.                  *
 * Range is (-inf,inf) f(x) = 0.5*exp(-|x|)                                  *
 *****************************************************************************/

double mt_dblexp_drand( mt_handle h ) {
  double u;
  mt_uint32 a, b;
  urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b); u = drand53_o(a,b);
  u += u;
  return (u<=1) ? log(u) : -log(u-1);
}

int mt_dblexp_drand_fill( mt_handle h, double * RESTRICT x, size_t n ) {
  double u;
  mt_uint32 a, b;
  if( h==NULL || x==NULL || n<1 ) return 1;
  for(;n;n--) {
    urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b); u = drand53_o(a,b);
    u += u;
    *(x++) = (u<=1) ? log(u) : -log(u-1);
  }
  return 0;
}

/*****************************************************************************
 * Generate a Gumbel distributed random number (minimum case).               *
 * Range is (-inf,inf), f(x) = exp(x)exp(-exp(x))                            *
 *****************************************************************************/

double mt_gumbel_drand( mt_handle h ) {
  mt_uint32 a, b;
  urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b);
  return log(-log(drand53_o(a,b)));
}

int mt_gumbel_drand_fill( mt_handle h, double * RESTRICT x, size_t n ) {
  mt_uint32 a, b;
  if( h==NULL || x==NULL || n<1 ) return 1;
  for(;n;n--) {
    urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b);
    *(x++) = log(-log(drand53_o(a,b)));
  }
  return 0;
}

/*****************************************************************************
 * Generate a Weibull distributed random number. Range is (0,inf)            *
 * f(x) = gamma*(x^(gamma-1))*exp(-x^gamma)                                  *
 *****************************************************************************/

/* FIXME: Optimize integer and sqrt cases of pow?? */

double mt_weibull_drand( mt_handle h, double gamma ) {
  double rgamma = 1/gamma;
  mt_uint32 a, b;
  urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b);
  return pow(-log(drand53_o(a,b)),rgamma);
}

int mt_weibull_drand_fill( mt_handle h, double gamma,
                           double * RESTRICT x, size_t n ) {
  double rgamma = 1/gamma;
  mt_uint32 a, b;
  if( h==NULL || x==NULL || n<1 ) return 1;
  for(;n;n--) {
    urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b);
    *(x++) = pow(-log(drand53_o(a,b)),rgamma);
  }
  return 0;
}

/********************************
 * Other distribution functions *
 ********************************/

/*****************************************************************************
 * Generate a Cauchy distributed random number. Range is (-inf,inf)          *
 * f(x) = 1/ ( pi*(1+x^2) )                                                  *
 *****************************************************************************/

double mt_cauchy_drand( mt_handle h ) {
  mt_uint32 a, b;
  urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b);
  return tan( M_PI*( drand53_o(a,b)-0.5 ) );
}

int mt_cauchy_drand_fill( mt_handle h, double * RESTRICT x, size_t n ) {
  mt_uint32 a, b;
  if( h==NULL || x==NULL || n<1 ) return 1;
  for(;n;n--) {
    urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b);
    *(x++) = tan( M_PI*( drand53_o(a,b)-0.5 ) );
  }
  return 0;
}

/*****************************************************************************
 * Generate a Tukey-lambda distributed random number. Range is (-inf,inf)    *
 * f(x) has no simple closed form                                            *
 *****************************************************************************/

double mt_lambda_drand( mt_handle h, double lambda ) {
  double u;
  mt_uint32 a, b;
  urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b); u = drand53_o(a,b);
  if      ( lambda==-1   ) u = 1/(1-u) - 1/u;
  else if ( lambda==-0.5 ) u = 1/sqrt(1-u) - 1/sqrt(u), u += u;
  else if ( lambda==0    ) u = log(u/(1-u));
  else if ( lambda==0.5  ) u = sqrt(u) - sqrt(1-u), u += u;
  else if ( lambda==1    ) u += u - 1;
  else                     u = ( pow(u,lambda) - pow(1-u,lambda) )/lambda;
  return u;
}

int mt_lambda_drand_fill( mt_handle h, double lambda,
                          double * RESTRICT x, size_t n ) {
  double u;
  mt_uint32 a, b;
  if( h==NULL || x==NULL || n<1 ) return 1;
  if(        lambda==-1   ) {
    for(;n;n--) {
      urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b); u = drand53_o(a,b);
      *(x++) = 1/(1-u) - 1/u;
    }
  } else if( lambda==-0.5 ) {
    for(;n;n--) {
      urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b); u = drand53_o(a,b);
      u = 1/sqrt(1-u) - 1/sqrt(u);
      *(x++) = u + u;
    }
  } else if( lambda==0    ) {
    for(;n;n--) {
      urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b); u = drand53_o(a,b);
      *(x++) = log(u/(1-u));
    }
  } else if( lambda==0.5  ) {
    for(;n;n--) {
      urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b); u = drand53_o(a,b);
      u = sqrt(u) - sqrt(1-u);
      *(x++) = u + u;
    }
  } else if( lambda==1    ) {
    for(;n;n--) {
      urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b); u = drand53_o(a,b);
      *(x++) = u + u - 1;
    }
  } else {
    lambda = 1/lambda;
    for(;n;n--) {
      urand32((mt_gen_t *)h,a); urand32((mt_gen_t *)h,b); u = drand53_o(a,b);
      *(x++) = ( pow(u,lambda) - pow(1-u,lambda) ) * lambda;
    }
  }
  return 0;
}

/*****************************************************************************
 * Fill x with a random permutation of the integers {0,1,2,...,n-2,n-1}      *
 * Routine generally valid for 1 <= n <= min(2^32,int_max)                   *
 *****************************************************************************/

int mt_randperm( mt_handle h, int * RESTRICT x, int n ) {
  mt_uint32 a, d;
  int i, t, r;

  /* Error checking */
  if( h==NULL || x==NULL || n<1 ) return 1;

  /* Create the initial permutation */
  for( i=0; i<n; i++ ) x[i] = i;

  /* Apply a random swap to each element of the permutation such that any of
     the n! permutations could be generated with equal probability. Note:
     The method used to pick a random number on [0...t-1] ([0...n-i-1]) uses
     the highest quality (most significant) bits of a urand32 to determine an
     _exactly_ uniformly distributed rand. A simpler "mod n" method uses the
     least significant bits and is not exactly uniform unless 2^32 is an exact
     multiple of t. Note: the d calculation is done with double precision
     because 2^32 cannot be represented with 32-bit ints and not all 32-bit
     integers have an exact single precision representation. */
  for( i=0; i<n-1; i++ ) {
    t = n-i;
    d = (mt_uint32)((double)4294967296./(double)t);
    do {
      urand32((mt_gen_t *)h,a);
      r = a/d;
    } while( r>=t );
    r += i;
    t = x[i];
    x[i] = x[r];
    x[r] = t;
  }

  return 0;
}

/*****************************************************************************
 * Shuffle the array x. x has n members. Each member has size s              *
 * This routine works for 1 <= n <= min(2^32,size_t_max)                     *
 *****************************************************************************/

int mt_shuffle( mt_handle h, void * RESTRICT x, size_t n, size_t s ) {
  mt_uint32 a, d;
  size_t i, t, r;
  char *xi, *xr, c;

  /* Error checking */
  if( h==NULL || x==NULL || n<1 || s<1 ) return 1;

  /* See randperm comment */
  for( i=0; i<n-1; i++ ) {
    t = n-i;
    d = (size_t)((double)4294967296./(double)t);
    do {
      urand32((mt_gen_t *)h,a);
      r = a/d;
    } while( r>=t );
    r += i;
    xi = ((char *)x) + s*i;
    xr = ((char *)x) + s*r;
    for( t=s; t; t-- ) {
      c       = *xi;
      *(xi++) = *xr;
      *(xr++) = c;
    }
  }

  return 0;
}
