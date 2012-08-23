#ifndef DaCSUtils_h
#define DaCSUtils_h

#include <cstdlib>

#if defined(cplusplus__)
	#include <iostream>
	#if defined HOST_BUILD
		#define DACS_PRINT_ERR(e, f, l) std::cerr << "ERROR HOST " << \
			(e) << " " << (f) << " " << (l) << std::endl
		#define DACS_PRINT(s) std::cout << "HOST " << (s) << std::endl
	#else
		#define DACS_PRINT_ERR(e, f, l) std::cerr << "ERROR ACCEL " << \
			(e) << " " << (f) << " " << (l) << std::endl
		#define DACS_PRINT(s) std::cout << "ACCEL " << (s) << std::endl
	#endif // BUILD
#else
	#include <stdio.h>
	#if defined HOST_BUILD
		#define DACS_PRINT_ERR(e, f, l) \
			printf("ERROR HOST %s %s %d\n", (e), (f), (l))
		#define DACS_PRINT(s) printf("HOST %s\n", (s))
	#else
		#define DACS_PRINT_ERR(e, f, l) \
			printf("ERROR ACCEL %s %s %d\n", (e), (f), (l))
		#define DACS_PRINT(s) printf("ACCEL %s\n", (s))
	#endif // BUIL
#endif // cplusplus__

inline void output_swap_type(DACS_BYTE_SWAP_T type) {
	switch(type) {
		case DACS_BYTE_SWAP_DISABLE:
			DACS_PRINT("DACS_BYTE_SWAP_DISABLE");
			break;
		case DACS_BYTE_SWAP_HALF_WORD:
			DACS_PRINT("DACS_BYTE_SWAP_WORD");
			break;
		case DACS_BYTE_SWAP_WORD:
			DACS_PRINT("DACS_BYTE_SWAP_WORD");
			break;
		case DACS_BYTE_SWAP_DOUBLE_WORD:
			DACS_PRINT("DACS_BYTE_SWAP_DOUBLE_WORD");
			break;
		default:
			DACS_PRINT("Unknown Swap Type");
			break;
	} // switch
} // output_swap_type

inline void process_dacs_errcode(DACS_ERR_T errcode,
	const char * file, int line) {
    if(errcode <= -35000) {
        DACS_PRINT_ERR(dacs_strerror(errcode), file, line);
        exit(errcode);
    } // if
} // process_dacs_errcode

#endif // DaCSUtils_h
