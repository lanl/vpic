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

#ifndef SystemRAM_hxx
#define SystemRAM_hxx

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

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
			std::cerr << "Available RAM (kilobytes): " <<
				available() << std::endl;
		} // print_available

		//! Report the available RAM on the system in kilobytes.
		static inline size_t available() {
			char buffer[81];
			std::ifstream meminfo("/proc/meminfo", std::ifstream::in);
			
			// Make sure that we were able to open the file
			if(meminfo.fail()) {
				std::cerr << "Failed opening /proc/meminfo file!!!" << std::endl;
				exit(1);
			} // if

			// Get the MemFree line
			meminfo.getline(buffer, 81);
			meminfo.getline(buffer, 81);

			// Parse out the free mem in kilobytes
			std::string memfree = buffer;
			size_t begin = memfree.find_first_not_of("MemFr: ");
			size_t end = memfree.find(" ", begin);

			// Convert to size_t
			size_t kilobytes;
			if(!from_string<size_t>(kilobytes, memfree.substr(begin, end-begin),
				std::dec)) {
				std::cerr << "Failed opening /proc/meminfo file!!!" << std::endl;
				exit(1);
			} // if

			return kilobytes;
		} // available
	}; // class SystemRAM

#endif // SystemRAM_hxx
