#ifndef _pipelines_h_
#define _pipelines_h_

//Redirect pipelines macros depending on pthreads vs. openmp

//----------------------------------------------------------------------------//
// Make sure that pipelines_pthreads.h and pipelines_openmp.h can only be
// included via this header file.
//----------------------------------------------------------------------------//

#define THREAD_REROUTE

//----------------------------------------------------------------------------//
// If using Pthreads, include pipelines_pthreads.h.
//----------------------------------------------------------------------------//

#if defined(VPIC_USE_PTHREADS)

#include "pipelines_pthreads.h"

//----------------------------------------------------------------------------//
// If using OpenMP, include pipelines_openmp.h.
//----------------------------------------------------------------------------//

#elif defined(VPIC_USE_OPENMP)

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
