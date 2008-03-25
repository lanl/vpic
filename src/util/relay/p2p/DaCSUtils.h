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
	char * file, int line) {
	
	switch(errcode) {
		case DACS_SUCCESS:
		case DACS_STS_PROC_FINISHED:
			return;
		case DACS_WID_BUSY:
			DACS_PRINT_ERR("DACS_WID_BUSY", file, line);
			break;
		case DACS_STS_PROC_RUNNING:
			DACS_PRINT_ERR("DACS_STS_PROC_RUNNING", file, line);
			break;
		case DACS_STS_PROC_FAILED:
			DACS_PRINT_ERR("DACS_STS_PROC_FAILED", file, line);
			break;
		case DACS_STS_PROC_ABORTED:
			DACS_PRINT_ERR("DACS_STS_PROC_ABORTED", file, line);
			break;
		case DACS_LAST_STATUS:
			DACS_PRINT_ERR("DACS_LAST_STATUS", file, line);
			break;
		case DACS_ERR_FIRST_ERROR:
			DACS_PRINT_ERR("DACS_ERR_FIRST_ERROR", file, line);
			break;
		case DACS_ERR_INTERNAL:
			DACS_PRINT_ERR("DACS_ERR_INTERNAL", file, line);
			break;
		case DACS_ERR_SYSTEM:
			DACS_PRINT_ERR("DACS_ERR_SYSTEM", file, line);
			break;
		case DACS_ERR_INVALID_ARGV:
			DACS_PRINT_ERR("DACS_ERR_INVALID_ARGV", file, line);
			break;
		case DACS_ERR_INVALID_ENV:
			DACS_PRINT_ERR("DACS_ERR_INVALID_ENV", file, line);
			break;
		case DACS_ERR_INVALID_HANDLE:
			DACS_PRINT_ERR("DACS_ERR_INVALID_HANDLE", file, line);
			break;
		case DACS_ERR_INVALID_ADDR:
			DACS_PRINT_ERR("DACS_ERR_INVALID_ADDR", file, line);
			break;
		case DACS_ERR_INVALID_ATTR:
			DACS_PRINT_ERR("DACS_ERR_INVALID_ATTR", file, line);
			break;
		case DACS_ERR_INVALID_DE:
			DACS_PRINT_ERR("DACS_ERR_INVALID_DE", file, line);
			break;
		case DACS_ERR_INVALID_PID:
			DACS_PRINT_ERR("DACS_ERR_INVALID_PID", file, line);
			break;
		case DACS_ERR_INVALID_TARGET:
			DACS_PRINT_ERR("DACS_ERR_INVALID_TARGET", file, line);
			break;
		case DACS_ERR_BUF_OVERFLOW:
			DACS_PRINT_ERR("DACS_ERR_BUF_OVERFLOW", file, line);
			break;
		case DACS_ERR_NOT_ALIGNED:
			DACS_PRINT_ERR("DACS_ERR_NOT_ALIGNED", file, line);
			break;
		case DACS_ERR_INVALID_SIZE:
			DACS_PRINT_ERR("DACS_ERR_INVALID_SIZE", file, line);
			break;
		case DACS_ERR_BYTESWAP_MISMATCH:
			DACS_PRINT_ERR("DACS_ERR_BYTESWAP_MISMATCH", file, line);
			break;
		case DACS_ERR_NO_RESOURCE:
			DACS_PRINT_ERR("DACS_ERR_NO_RESOURCE", file, line);
			break;
		case DACS_ERR_PROC_LIMIT:
			DACS_PRINT_ERR("DACS_ERR_PROC_LIMIT", file, line);
			break;
		case DACS_ERR_NO_PERM:
			DACS_PRINT_ERR("DACS_ERR_NO_PERM", file, line);
			break;
		case DACS_ERR_OWNER:
			DACS_PRINT_ERR("DACS_ERR_OWNER", file, line);
			break;
		case DACS_ERR_NOT_OWNER:
			DACS_PRINT_ERR("DACS_ERR_NOT_OWNER", file, line);
			break;
		case DACS_ERR_RESOURCE_BUSY:
			DACS_PRINT_ERR("DACS_ERR_RESOURCE_BUSY", file, line);
			break;
		case DACS_ERR_GROUP_CLOSED:
			DACS_PRINT_ERR("DACS_ERR_GROUP_CLOSED", file, line);
			break;
		case DACS_ERR_GROUP_OPEN:
			DACS_PRINT_ERR("DACS_ERR_GROUP_OPEN", file, line);
			break;
		case DACS_ERR_GROUP_DUPLICATE:
			DACS_PRINT_ERR("DACS_ERR_GROUP_DUPLICATE", file, line);
			break;
		case DACS_ERR_INVALID_WID:
			DACS_PRINT_ERR("DACS_ERR_INVALID_WID", file, line);
			break;
		case DACS_ERR_INVALID_STREAM:
			DACS_PRINT_ERR("DACS_ERR_INVALID_STREAM", file, line);
			break;
		case DACS_ERR_NO_WIDS:
			DACS_PRINT_ERR("DACS_ERR_NO_WIDS", file, line);
			break;
		case DACS_ERR_WID_ACTIVE:
			DACS_PRINT_ERR("DACS_ERR_WID_ACTIVE", file, line);
			break;
		case DACS_ERR_WID_NOT_ACTIVE:
			DACS_PRINT_ERR("DACS_ERR_WID_NOT_ACTIVE", file, line);
			break;
		case DACS_ERR_INITIALIZED:
			DACS_PRINT_ERR("DACS_ERR_INITIALIZED", file, line);
			break;
		case DACS_ERR_NOT_INITIALIZED:
			DACS_PRINT_ERR("DACS_ERR_NOT_INITIALIZED", file, line);
			break;
		case DACS_ERR_MUTEX_BUSY:
			DACS_PRINT_ERR("DACS_ERR_MUTEX_BUSY", file, line);
			break;
		case DACS_ERR_NOT_SUPPORTED_YET:
			DACS_PRINT_ERR("DACS_ERR_NOT_SUPPORTED_YET", file, line);
			break;
		case DACS_ERR_VERSION_MISMATCH:
			DACS_PRINT_ERR("DACS_ERR_VERSION_MISMATCH", file, line);
			break;
		case DACS_ERR_DACSD_FAILURE:
			DACS_PRINT_ERR("DACS_ERR_DACSD_FAILURE", file, line);
			break;
		case DACS_ERR_INVALID_PROG:
			DACS_PRINT_ERR("DACS_ERR_INVALID_PROG", file, line);
			break;
		case DACS_ERR_LAST_ERROR:
			DACS_PRINT_ERR("DACS_ERR_LAST_ERROR", file, line);
			break;
		default:
			DACS_PRINT_ERR("Unknown Error", file, line);
			break;
	} // switch

	exit(1);
} // process_dacs_errcode

#endif // DaCSUtils_h
