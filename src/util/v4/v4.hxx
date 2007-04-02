#ifndef v4_hxx
#define v4_hxx

#if defined __SSE__
	#include <v4_sse.hxx>
#elif defined __ALTIVEC__
#elif defined __SPU__
	#include <v4_spu.hxx>
#else
	#include <v4_portable.hxx>
#endif // ISA

#endif // v4_hxx
