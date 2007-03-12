#ifndef spe_thread_h
#define spe_thread_h

#if defined __PPU__

#include <libspe2.h>

/* arguments for running an SPE thread */
typedef struct {
	struct spe_context * spe;
	void * argp;
	void * envp;
} thread_args;

#endif // __PPU__

#endif // spe_thread_h
