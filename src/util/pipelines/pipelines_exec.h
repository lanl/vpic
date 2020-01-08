#ifndef _pipelines_exec_h_
#define _pipelines_exec_h_

//----------------------------------------------------------------------------//
// Include some stuff that is common to both Pthreads and OpenMP.
//----------------------------------------------------------------------------//

#include "pipelines.h"

#include "../v16/v16.h"
#include "../v4/v4.h"
#include "../v8/v8.h"

//----------------------------------------------------------------------------//
// Make sure that pipelines_exec_pth.h and pipelines_exec_omp.h can only be
// included via this header file.
//----------------------------------------------------------------------------//

#define THREAD_REROUTE

//----------------------------------------------------------------------------//
// If using Pthreads, include pipelines_exec_pth.h.
//----------------------------------------------------------------------------//

#if defined( VPIC_USE_PTHREADS )

#include "pipelines_exec_pth.h"

//----------------------------------------------------------------------------//
// If using OpenMP, include pipelines_exec_omp.h.
//----------------------------------------------------------------------------//

#elif defined( VPIC_USE_OPENMP )

#include "pipelines_exec_omp.h"

//----------------------------------------------------------------------------//
// I wonder if VPIC will actually run without a threading model.
//----------------------------------------------------------------------------//

#else

// Need to figure out how to handle this case.

#error "VPIC_USE_OPENMP or VPIC_USE_PTHREADS must be specified"

#endif

//----------------------------------------------------------------------------//
// Undefine local macros.
//----------------------------------------------------------------------------//

#undef THREAD_REROUTE

#endif // _pipelines_exec_h_
