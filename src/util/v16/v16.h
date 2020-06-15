#ifndef _v16_h_
#define _v16_h_
/* FIXME: STYLE */
#define IN_v16_h
/* FIXME: SHOULDN'T THIS INCLUDE UTIL_BASE.H? */
#ifdef __cplusplus
#if defined USE_V16_PORTABLE
#include "v16_portable.h"
#elif defined USE_V16_AVX512
#include "v16_avx512.h"
#endif
#endif
#undef IN_v16_h
#endif // _v16_h_
