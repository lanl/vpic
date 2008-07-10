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

// Integer to random floating point conversions

// The c0, c1 and c variants have very rigorous interpretations.
// Imagine a real random number is generated on [0,1).  The result is
// then rounded to the 2^53+1 point uniform lattice that covers this
// interval (endpoints inclusive).  The c0 variant corresponds to
// rounding the real random down to the nearest lattice point (as
// such, the result 1 will never occur).  The c1 variant corresponds
// to rounding the real random up to the nearest lattice point (as
// such, the result 0 will never occur).  The c variant corresponds to
// rounding to even the real random (as such, both 0 and 1 could
// occur).
//
// The o variant does not have nearly as rigorous an interpretation.
// The simple approach would be to take a 53-bit rand, add 0.5 and
// multiply the result by 1/2^53. But 2^53-0.5 doesn't have an
// e.r. 2^53 does though so we add alpha to the 53-bit rand where
// 1<alpha<2. Also, we can't divide by 2^53+1 ... no e.r, but 2^53+2
// does have an e.r. so we use that. alpha is picked to center the
// 53-bit integer rand interval between [0, 2^53+2]. Some of these
// macros could be simplified if the FPU computes intermediates in
// higher precision (i.e. the x86 FPU) but this is not assumed
//
// Similar issues apply for the frand24 and drand32 macros.

// 24-bit single precision rands

#define frand24_c0(u32) (( (u32)>>8           )*(1.f/16777216.f))
#define frand24_c1(u32) ((((u32)>>8)+1        )*(1.f/16777216.f))
#define frand24_c(u32)  ((((u32)>>8)+((u32)&1))*(1.f/16777216.f))
#define frand24_o(u32)  ((((u32)>>8)+1.5f     )*(1.f/16777218.f) )

// 32-bit double precision rands. Note: Because doubles have more
// precision, these macros are much less tricky.  Note that because a
// is an unsigned 32-bit, we cannot safely add until we prompote a to
// a double.

#define drand32_c0(u32) ( ((double)(u32))           *(1./4294967296.))
#define drand32_c1(u32) ((((double)(u32))+1.0      )*(1./4294967296.))
#define drand32_c(u32)  ((((double)(u32))+((u32)&1))*(1./4294967296.))
#define drand32_o(u32)  ((((double)(u32))+0.5      )*(1./4294967297.))

// 53-bit double precision rands

#define drand53_c0(a,b) ((((a)>>5)*67108864.+ ((b)>>6)         )*(1./9007199254740992.))
#define drand53_c1(a,b) ((((a)>>5)*67108864.+(((b)>>6)+1      ))*(1./9007199254740992.))
#define drand53_c(a,b)  ((((a)>>5)*67108864.+(((b)>>6)+((b)&1)))*(1./9007199254740992.))
#define drand53_o(a,b)  ((((a)>>5)*67108864.+ ((b)>>6)+1.5     )*(1./9007199254740994.))

#endif // _mtrand_conv_h_
