#ifndef v4_hxx
#define v4_hxx

/* pick the appropriate ISA specific implementation */
#if defined SCALAR
	#include <v4_portable.hxx>
#elif defined __SSE__
	#include <v4_sse.hxx>
#elif defined __ALTIVEC__
	#include <v4_altivec.hxx>
#elif defined __SPU__
	#include <v4_spu.hxx>
#endif // ISA

#endif // v4_hxx
