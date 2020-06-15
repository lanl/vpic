#ifndef _pipelines_h_
#define _pipelines_h_

//----------------------------------------------------------------------------//
// Include some stuff that is common to both Pthreads and OpenMP.
//----------------------------------------------------------------------------//

#include "../util_base.h"

enum
{
    MAX_PIPELINE = 272
};

// Is this even related to pipelines.  Maybe this should be in util_base.h.
#define PAD_STRUCT( sz )

//----------------------------------------------------------------------------//
// Make sure that pipelines_pthreads.h and pipelines_openmp.h can only be
// included via this header file.
//----------------------------------------------------------------------------//

#define THREAD_REROUTE

//----------------------------------------------------------------------------//
// If using Pthreads, include pipelines_pthreads.h.
//----------------------------------------------------------------------------//

#if defined( VPIC_USE_PTHREADS )

#include "pipelines_pthreads.h"

//----------------------------------------------------------------------------//
// If using OpenMP, include pipelines_openmp.h.
//----------------------------------------------------------------------------//

#elif defined( VPIC_USE_OPENMP )

#include "pipelines_openmp.h"

//----------------------------------------------------------------------------//
// I wonder if VPIC will actually run without a threading model.
//----------------------------------------------------------------------------//

#else

// Need to figure out how to handle this case.

#error "VPIC_USE_OPENMP or VPIC_USE_PTHREADS must be specified"

#endif

//----------------------------------------------------------------------------//
// Make sure that pipelines_pthreads.h and pipelines_openmp.h can only be
// included via this header file.
//----------------------------------------------------------------------------//

#undef THREAD_REROUTE

#endif // _pipelines_h_
