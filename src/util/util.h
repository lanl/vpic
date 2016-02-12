#ifndef _util_h_
#define _util_h_

// Expose all public functionality in util.  The below includes bring
// in util_base.h and other low level includes automatically.

#include "v4/v4.h" // Must be first (FIXME: REALLY?)
#include "v8/v8.h"
#include "checkpt/checkpt.h"
#include "mp/mp.h"
#include "rng/rng.h"
#include "pipelines/pipelines.h"
#include "profile/profile.h"

BEGIN_C_DECLS

// Boot all util functionality (should be the first thing in the program)

void
boot_services( int * pargc,
               char *** pargv );

// Halt all util functionality (should be the last thing in the program)

void
halt_services( void );

// Give an estimate of the time when boot_timestamp was called
// (in seconds since the epoch).  All processes agree on this.

#define boot_timestamp (double)_boot_timestamp
extern double _boot_timestamp;

// Give an estimate of how many seconds since boot_services was
// called (in seconds).  This call must be loosly synchronous over
// all processes; all processes agree on the result.

double
uptime( void );

END_C_DECLS

#endif // _util_h_
