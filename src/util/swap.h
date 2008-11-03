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

#include <field_advance.h>
#include <sf_interface.h>
#include <species_advance.h>

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

//template<typename T> void inline swap(T & element);
void inline swap(char & element) {}

void inline swap(short int & element) {
	element = bswap_16(element);
} // swap

//template<> void inline swap<double>(double & element) {
void inline swap(double & element) {

	union type64 {
		type64(double d_) : d(d_) {}

		double  d;
		uint64_t ui;
	} t64(element);

	t64.ui = bswap_64(t64.ui);
	element = t64.d;
} // swap

union type32 {
	type32(float f_) : f(f_) {}

	float f;
	uint32_t ui;
};

float inline swap_float(float f) {
	union type32 {
		type32(float f_) : f(f_) {}

		float f;
		uint32_t ui;
	} t32(f);

	t32.ui = bswap_32(t32.ui);
	return t32.f;
} // swap_float

//template<> void inline swap<uint64_t>(uint64_t & element) {
void inline swap(uint64_t & element) {
	element = bswap_64(element);
} // swap

//template<> void inline swap<int64_t>(int64_t & element) {
void inline swap(int64_t & element) {
	element = bswap_64(element);
} // swap

//template<> void inline swap<float>(float & element) {
void inline swap(float & element) {

	union type32 {
		type32(float f_) : f(f_) {}

		float f;
		uint32_t ui;
	} t32(element);

	t32.ui = bswap_32(t32.ui);
	element = t32.f;
} // swap

void inline swap(uint16_t & element) {
	element = bswap_16(element);
} // swap

//template<> void inline swap<uint32_t>(uint32_t & element) {
void inline swap(uint32_t & element) {
	element = bswap_32(element);
} // swap

//template<> void inline swap<int32_t>(int32_t & element) {
void inline swap(int32_t & element) {
	element = bswap_32(element);
} // swap

void inline swap(field_t & element) {
	// electric field
	utils::swap(element.ex);
	utils::swap(element.ey);
	utils::swap(element.ez);
	utils::swap(element.div_e_err);

	// magnetic field
	utils::swap(element.cbx);
	utils::swap(element.cby);
	utils::swap(element.cbz);
	utils::swap(element.div_b_err);

	// tca field
	utils::swap(element.tcax);
	utils::swap(element.tcay);
	utils::swap(element.tcaz);
	utils::swap(element.rhob);

	// tca field
	utils::swap(element.jfx);
	utils::swap(element.jfy);
	utils::swap(element.jfz);
	utils::swap(element.rhof);

	// material
	utils::swap(element.ematx);
	utils::swap(element.ematy);
	utils::swap(element.ematz);
	utils::swap(element.nmat);

	// material
	utils::swap(element.fmatx);
	utils::swap(element.fmaty);
	utils::swap(element.fmatz);
	utils::swap(element.cmat);
} // swap

void inline swap(hydro_t & element) {
	// current and charge
	utils::swap(element.jx);
	utils::swap(element.jy);
	utils::swap(element.jz);
	utils::swap(element.rho);

	// current and charge
	utils::swap(element.px);
	utils::swap(element.py);
	utils::swap(element.pz);
	utils::swap(element.ke);

	// stress diag
	utils::swap(element.txx);
	utils::swap(element.tyy);
	utils::swap(element.tzz);

	// stress off-diag
	utils::swap(element.tyz);
	utils::swap(element.tzx);
	utils::swap(element.txy);
} // swap

void inline swap(particle_t & element) {
	// position
	utils::swap(element.dx);
	utils::swap(element.dy);
	utils::swap(element.dz);

	// id
	utils::swap(element.i);

	// momentum
	utils::swap(element.ux);
	utils::swap(element.uy);
	utils::swap(element.uz);

	// charge
	utils::swap(element.q);
} // swap

} // namespace utils
#endif // __cplusplus

#endif // swap_h
