#ifndef _v4_h_
#define _v4_h_

/* FIXME: STYLE */
#define IN_v4_h

/* FIXME: SHOULDN'T THIS INCLUDE UTIL_BASE.H? */

#ifdef __cplusplus

# if defined USE_V4_ALTIVEC
#   include "v4_altivec.h"

# elif defined USE_V4_PORTABLE
#   include "v4_portable.h"

# elif defined USE_V4_SSE
#   include "v4_sse.h"

# elif defined USE_V4_AVX
#   include "v4_avx.h"

# elif defined USE_V4_AVX2
#   include "v4_avx2.h"

# elif defined USE_V4_NEON
#   include "v4_neon.h"

# endif

#endif

#undef IN_v4_h

#endif // _v4_h_
