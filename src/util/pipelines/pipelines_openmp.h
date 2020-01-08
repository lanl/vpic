#ifndef _pipelines_openmp_h_
#define _pipelines_openmp_h_

#ifndef THREAD_REROUTE
#error "Do not include pipelines_openmp.h directly; use pipelines.h"
#endif

// #include "../util_base.h"

#include <omp.h>

//----------------------------------------------------------------------------//
// Generic macros that are used for all cases of vector acceleration as well
// as the standard case that does not use vector acceleration.
//----------------------------------------------------------------------------//

#define N_PIPELINE omp_helper.n_pipeline

//----------------------------------------------------------------------------//
// A container object to mimic 'serial' and 'thread' objects.  Currently just
// useful for avoiding global var 'n_pipeline'.
//----------------------------------------------------------------------------//

typedef struct omp_container
{
    int n_pipeline;
    int dispatch_to_host;

    // const char * f_dump;
    // const char * e_dump;

    // boot gets the number of pipelines from the cmd line
    // and passes it on to the EXEC_PIPELINS macro eventually
    void ( *boot )( int* pargc, char*** pargv );
} omp_container_t;

BEGIN_C_DECLS

extern omp_container_t omp_helper;

END_C_DECLS

#endif // _pipelines_openmp_h_
