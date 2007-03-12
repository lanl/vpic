#if defined __PPU__

#include <stdlib.h>
#include <pthread.h>
#include <spe_thread.h>

/*
	Run the context pointed to by args->spe using the arguments in
	args->argp and args->envp.  The entry pointed SPE_DEFAULT_ENTRY just means
	to start at the beginning of the text segment.
*/
void * spe_thread(void * vargs) {
	thread_args * args = (thread_args *)vargs;
	unsigned int runflags = 0;
	unsigned int entry = SPE_DEFAULT_ENTRY;

	spe_context_run(args->spe, &entry, runflags, args->argp, args->envp, NULL);
	pthread_exit(NULL);
} /* spe_thread */

#endif // __PPU__
