#ifndef _vfa_v4_private_h_
#define _vfa_v4_private_h_

// Vacuum field advance implementation

#ifndef IN_vfa_v4
#error "Do not include vfa_v4_private.h; include field_advance.h"
#endif

// Note: sfa_v4_pipeline include must come after above vfa_private.h
// include to get appropriate pipeline dispatcher

#define IN_vfa
#include "../vfa_private.h"

#define IN_sfa_v4
#include "../../v4/sfa_v4_private.h"

// Note: See field methods for a detailed description of these functions

// Note: TCA IS _NOT_ UPDATED UNDER VFA UNDER NORMAL OPERATION ONLY
// THE COMPUTE_CURL_B MECHANISM SETS IT (AND THAT HAPPENS DURING
// INITIALIZATION).

// FIXME: THE HOST PROCESSED FIELD KERNELS SHOULD BE UPDATED TO USE
// SCALAR FMA INSTRUCTIONS WITH COMMENSURATE ROUND-OFF PROPERTIES TO
// THE FMA INSTRUCTIONS USED ON THE PIPELINE PROCESSED FIELDS!

BEGIN_C_DECLS

// In vfa_v4_advance_e.cxx

void
vfa_v4_advance_e( field_t                      * ALIGNED(128) f,
                  const material_coefficient_t * ALIGNED(128) m, // Ignored
                  const grid_t                 *              g );

END_C_DECLS

#endif // _vfa_v4_private_h_
