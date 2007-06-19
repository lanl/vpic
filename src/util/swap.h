/*
	Byte-swapping utilities

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$
	vim: set ts=3 :
*/

#ifndef swap_h
#define swap_h

#if defined(__GNUC__)
	#include <byteswap.h>
#else
	#define bswap_16(x) ((((uint16_t)(x) & 0xff00u) >> 8) | \
		(((uint16_t)(x) & 0x00ffu) << 8))
	#define bswap_32(x) ((((uint32_t)(x) & 0xff000000u) >> 24) | \
		(((uint32_t)(x) & 0x00ff0000u) >> 8)  | \
		(((uint32_t)(x) & 0x0000ff00u) << 8)  | \
		(((uint32_t)(x) & 0x000000ffu) << 24))
	#define bswap_64(x) ((((uint64_t)(x) & 0xff00000000000000ull) >> 56) | \
		(((uint64_t)(x) & 0x00ff000000000000ull) >> 40) | \
		(((uint64_t)(x) & 0x0000ff0000000000ull) >> 24) | \
		(((uint64_t)(x) & 0x000000ff00000000ull) >> 8) | \
		(((uint64_t)(x) & 0x00000000ff000000ull) << 8) | \
		(((uint64_t)(x) & 0x0000000000ff0000ull) << 24) | \
		(((uint64_t)(x) & 0x000000000000ff00ull) << 40) | \
		(((uint64_t)(x) & 0x00000000000000ffull) << 56))
#endif // __GNUC__

#if defined __cplusplus
namespace utils {

template<typename T> void swap(T & element);

template<> void swap<double>(double & element) {

	union type64 {
		type64(double d_) : d(d_) {}

		double  d;
		uint64_t ui;
	} t64(element);

	t64.ui = bswap_64(t64.ui);
	element = t64.d;
} // swap

template<> void swap<uint64_t>(uint64_t & element) {
	element = bswap_64(element);
} // swap

template<> void swap<int64_t>(int64_t & element) {
	element = bswap_64(element);
} // swap

template<> void swap<float>(float & element) {

	union type32 {
		type32(float f_) : f(f_) {}

		float f;
		uint32_t ui;
	} t32(element);

	t32.ui = bswap_32(t32.ui);
	element = t32.f;
} // swap

template<> void swap<uint32_t>(uint32_t & element) {
	element = bswap_32(element);
} // swap

template<> void swap<int32_t>(int32_t & element) {
	element = bswap_32(element);
} // swap

} // namespace utils
#endif // __cplusplus

#endif // swap_h
