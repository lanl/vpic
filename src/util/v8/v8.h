#ifndef _v8_h_
#define _v8_h_
/* FIXME: STYLE */
#define IN_v8_h
/* FIXME: SHOULDN'T THIS INCLUDE UTIL_BASE.H? */
#ifdef __cplusplus
# if defined USE_V8_PORTABLE
#   include "v8_portable.h"
# elif defined USE_V8_AVX2
#   include "v8_avx2.h"
# elif defined USE_V8_AVX
#   include "v8_avx.h"
# endif
#endif
#undef IN_v8_h
#endif // _v8_h_
