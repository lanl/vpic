#ifndef _v4_h_
#define _v4_h_
#define IN_v4_h
#ifdef __cplusplus
# if defined USE_V4_PORTABLE
#   include <v4_portable.hxx>
# elif defined USE_V4_SSE
#   include <v4_sse.hxx>
# elif defined USE_V4_ALTIVEC
#   include <v4_altivec.hxx>
# elif defined USE_V4_SPU
#   include <v4_spu.hxx>
# endif
#endif
#undef IN_v4_h
#endif // _v4_h_
