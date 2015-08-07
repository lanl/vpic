/*------------------------------------------------------------------------------
 * Definition of SystemRAM class
 *
 * PLBM_DISTRIBUTION_HEADER
 *
 * $Revision:$
 * $Date:$
 * $Author:$
 *
 * vim: set ts=3 :
------------------------------------------------------------------------------*/

#ifndef SystemRAM_h
#define SystemRAM_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

#include <util_base.h>

// String to type conversion
template <typename T>
bool from_string(T & t, const std::string & s,
	std::ios_base & (*f)(std::ios_base&)) {
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
} // from_string


/*!
	\struct SystemRAM SystemRAM.h
	\brief SystemRAM provides...
*/
struct SystemRAM
	{
		static inline void print_available() {
			MESSAGE(("Available RAM (kilobytes): %ld", available()));
		} // print_available

		//! Report the available RAM on the system in kilobytes.
		static inline uint64_t available() {
			
			#if !__linux__
				ERROR(("SystemRAM: Unsupported Operating System!!!"));
			#endif

			char buffer[81];
			std::ifstream meminfo("/proc/meminfo", std::ifstream::in);
			
			// Make sure that we were able to open the file
			if(meminfo.fail()) {
				ERROR(("Failed opening /proc/meminfo file!!!"));
			} // if

			// Get the MemFree line
			meminfo.getline(buffer, 81);
			meminfo.getline(buffer, 81);

			meminfo.close();

			// Parse out the free mem in kilobytes
			std::string memfree = buffer;
			size_t begin = memfree.find_first_not_of("MemFr: ");
			size_t end = memfree.find(" ", begin);

			// Convert to size_t
			uint64_t kilobytes;
			if(!from_string<uint64_t>(kilobytes,
				memfree.substr(begin, end-begin),
				std::dec)) {
				ERROR(("String conversion to size_t failed!!!"));
			} // if

			return kilobytes;
		} // available
	}; // class SystemRAM

#endif // SystemRAM_h
