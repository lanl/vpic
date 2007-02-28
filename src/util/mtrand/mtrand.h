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

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#  define __BEGIN_DECLS extern "C" {
#  define __END_DECLS }
#else
#  define __BEGIN_DECLS /* empty */
#  define __END_DECLS   /* empty */
#endif

#undef __ARGLIST
#if defined (__STDC__)                            || \
    defined (_AIX)                                || \
    (defined (__mips) && defined (_SYSTYPE_SVR4)) || \
    defined(WIN32)                                || \
    defined(__cplusplus)
#  define __ARGLIST(args) args
#else
#  define __ARGLIST(args) ()
#endif

#include <stddef.h> /* Need size_t */

/* Setup the mt_handle datatype */
typedef void *mt_handle;

/* Note: All these trinary operators should be evaluated
   at compile time for any self-respecting compiler. */
#include <limits.h>
#define mt_CRAND_MAX   ((CHAR_BIT*sizeof(signed char)<32) ?             \
                        SCHAR_MAX : 0x7fffffff)
#define mt_SHRAND_MAX  ((CHAR_BIT*sizeof(short int)<32) ?               \
                        SHRT_MAX : 0x7fffffff)
#define mt_RAND_MAX    ((CHAR_BIT*sizeof(int)<32) ?                     \
                        INT_MAX : 0x7fffffff)
#define mt_LRAND_MAX   ((CHAR_BIT*sizeof(long int)<32) ?                \
                        LONG_MAX : 0x7fffffffL)
#define mt_UCRAND_MAX  ((CHAR_BIT*sizeof(unsigned char)<32) ?           \
                        UCHAR_MAX : 0xffffffffU)
#define mt_USHRAND_MAX ((CHAR_BIT*sizeof(unsigned short int)<32) ?      \
                        USHRT_MAX : 0xffffffffU)
#define mt_URAND_MAX   ((CHAR_BIT*sizeof(unsigned int)<32) ?            \
                        UINT_MAX : 0xffffffffU)
#define mt_ULRAND_MAX  ((CHAR_BIT*sizeof(unsigned long int)<32) ?       \
                        ULONG_MAX : 0xffffffffUL)

__BEGIN_DECLS

/********************************
 * Constructors and destructors *
 ********************************/

extern mt_handle mt_new_generator    __ARGLIST((unsigned int));
extern void      mt_delete_generator __ARGLIST((mt_handle *));
 
/*********************************
 * Seeders and related functions *
 *********************************/
extern int    mt_srand    __ARGLIST((mt_handle, unsigned int));
extern size_t mt_getsize  __ARGLIST((mt_handle));
extern int    mt_getstate __ARGLIST((mt_handle, char *, size_t));
extern int    mt_setstate __ARGLIST((mt_handle, const char *, size_t));

/**********************
 * Integer generators *
 **********************/
#define int_gen_decl( name, type )                                            \
extern signed   type mt_##name  __ARGLIST((mt_handle));                       \
extern unsigned type mt_u##name __ARGLIST((mt_handle));                       \
extern int mt_##name##_fill  __ARGLIST((mt_handle, signed   type *, size_t)); \
extern int mt_u##name##_fill __ARGLIST((mt_handle, unsigned type *, size_t))
int_gen_decl(crand,  char     );
int_gen_decl(shrand, short int);
int_gen_decl(rand,   int      );
int_gen_decl(lrand,  long int );
#undef int_gen_decl

/*************************************
 * Uniform floating point generators *
 *************************************/
#define float_gen_decl( name, type )					\
extern type mt_##name           __ARGLIST((mt_handle));			\
extern type mt_##name##_c0      __ARGLIST((mt_handle));			\
extern type mt_##name##_c1      __ARGLIST((mt_handle));			\
extern type mt_##name##_c       __ARGLIST((mt_handle));			\
extern int  mt_##name##_fill    __ARGLIST((mt_handle, type *, size_t));	\
extern int  mt_##name##_c0_fill __ARGLIST((mt_handle, type *, size_t));	\
extern int  mt_##name##_c1_fill __ARGLIST((mt_handle, type *, size_t));	\
extern int  mt_##name##_c_fill  __ARGLIST((mt_handle, type *, size_t))
float_gen_decl(frand,      float );
float_gen_decl(fast_drand, double);
float_gen_decl(drand,      double);
#undef float_gen_decl

/*************************
 * Speciality generators *
 *************************/
extern double mt_normal_drand      __ARGLIST((mt_handle));
extern int    mt_normal_drand_fill __ARGLIST((mt_handle, double *, size_t));

extern double mt_lognormal_drand      __ARGLIST((mt_handle, double));
extern int    mt_lognormal_drand_fill __ARGLIST((mt_handle, double,
						 double *, size_t));

extern double mt_bs_drand      __ARGLIST((mt_handle, double));
extern int    mt_bs_drand_fill __ARGLIST((mt_handle, double,
					  double *, size_t));

extern double mt_exp_drand      __ARGLIST((mt_handle));
extern int    mt_exp_drand_fill __ARGLIST((mt_handle, double *, size_t));

extern double mt_dblexp_drand      __ARGLIST((mt_handle));
extern int    mt_dblexp_drand_fill __ARGLIST((mt_handle, double *, size_t));

extern double mt_gumbel_drand      __ARGLIST((mt_handle));
extern int    mt_gumbel_drand_fill __ARGLIST((mt_handle, double *, size_t));

extern double mt_weibull_drand      __ARGLIST((mt_handle, double));
extern int    mt_weibull_drand_fill __ARGLIST((mt_handle, double,
					       double *, size_t));

extern double mt_cauchy_drand __ARGLIST((mt_handle));
extern int    mt_cauchy_drand_fill __ARGLIST((mt_handle, double *, size_t));

extern double mt_lambda_drand      __ARGLIST((mt_handle, double));
extern int    mt_lambda_drand_fill __ARGLIST((mt_handle, double,
					      double *, size_t));

extern int mt_randperm __ARGLIST((mt_handle, int *, int));

#define MT_SHUFFLE(h,x,n) mt_shuffle((h),(x),(n),sizeof(*(x)))
extern int mt_shuffle __ARGLIST((mt_handle, void *, size_t, size_t));

__END_DECLS

#endif
