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

/* sent integer message to the SPE pointed to by handle */
static inline void signal_spe(spe_context_ptr_t handle, unsigned int signal) {
	spe_in_mbox_write(handle, &signal, 1, SPE_MBOX_ALL_BLOCKING);
} // signal_spe

#endif // __PPU__

#endif // spe_thread_h
