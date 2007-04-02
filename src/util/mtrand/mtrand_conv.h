/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Adapted from earlier V4PIC versions.
 *
 */

#ifndef _mtrand_conv_h_
#define _mtrand_conv_h_

/* Integer to floating point conversions */

/* The below is tricky. In drand53_c, multiplying by the inverse causes a
   0.5*eps error that prevents 1 from ever being attained. drand53_c uses a
   division instead, which is unfortunately slower. The multiply has an error
   because 1/(2^53-1) does not have an exact representation (e.r.) in floating
   point but 2^53-1 does. In the c0 and c1 cases, 1/2^53 does have an e.r. so
   multiplying is preferred for speed. drand53_o is even trickier. Ideally you
   should take a 53-bit rand, add 0.5 and multiply the result by 1/2^53. But
   2^53-0.5 doesn't have an e.r. 2^53 does though so we add alpha to the
   53-bit rand where 1<alpha<2. Also, we can't divide by 2^53+1 ... no e.r,
   but 2^53+2 does have an e.r. so we use that. alpha is picked to center the
   53-bit integer rand interval between [0, 2^53+2]. Some of these macros
   could be simplified if the FPU computes intermediates in higher precision
   (i.e. the x86 FPU) but this is not assumed. Similar issues apply for
   the frand24 and drand32 macros. */

/* 24-bit single precision rands */
#define frand24_c(a)   (((a)>>8)/16777215.0)                      /* [0,1] */
#define frand24_c0(a)  (((a)>>8)*(1.0/16777216.0))                /* [0,1) */
#define frand24_c1(a)  (1-frand24_c0(a))                          /* (0,1] */
#define frand24_o(a)   ((((a)>>8)+1.5)*(1.0/16777218.0))          /* (0,1) */

/* 32-bit double precision rands. Note: Because doubles have more precision,
   these macros are much less tricky. Still need a division in drand32_c. */
#define drand32_c(a)   ((a)/4294967295.0)
#define drand32_c0(a)  ((a)*(1.0/4294967296.0))
#define drand32_c1(a)  (1-drand32_c0(a))
#define drand32_o(a)   (((a)+0.5)*(1.0/4294967296.0))

/* 53-bit double precision rands */
#define drand53_i(a,b)  (((a)>>5)*67108864.0+((b)>>6))
#define drand53_c(a,b)  (drand53_i(a,b)/9007199254740991.0)
#define drand53_c0(a,b) (drand53_i(a,b)*(1.0/9007199254740992.0))
#define drand53_c1(a,b) (1-drand53_c0(a,b))
#define drand53_o(a,b)  ((drand53_i(a,b)+1.5)*(1.0/9007199254740994.0))

#endif /* _mtrand_conv_h_ */
