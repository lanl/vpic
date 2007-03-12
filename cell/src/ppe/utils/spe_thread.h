#ifndef spe_thread_h
#define spe_thread_h

#include <libspe2.h>

/* arguments for running an SPE thread */
typedef struct {
	struct spe_context * spe;
	void * argp;
	void * envp;
} thread_args;

/* function for running an SPE thread using pthreads */
void * spe_thread(void * vargs);

#endif // spe_thread_h
