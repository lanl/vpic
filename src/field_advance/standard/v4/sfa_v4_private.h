#ifndef _sfa_v4_private_h_
#define _sfa_v4_private_h_

// Standard v4 accelerated field advance implementation

#ifndef IN_sfa_v4
#error "Do not include sfa_v4_private.h; include field_advance.h"
#endif

#define IN_sfa
#include "../sfa_private.h"

// FIXME: THE HOST PROCESSED FIELD KERNELS SHOULD BE UPDATED TO USE
// SCALAR FMA INSTRUCTIONS WITH COMMENSURATE ROUND-OFF PROPERTIES TO
// THE FMA INSTRUCTIONS USED ON THE PIPELINE PROCESSED FIELDS!

BEGIN_C_DECLS

// In v4_advance_b.cxx

void
v4_advance_b( field_t      * ALIGNED(128) f,
              const grid_t *              g,
              float                       frac );

// In v4_advance_e.c

void
v4_advance_e( field_t                      * ALIGNED(128) f,
              const material_coefficient_t * ALIGNED(128) m,
              const grid_t                 *              g );

// In v4_compute_curl_b.c

void
v4_compute_curl_b( field_t                      * ALIGNED(128) f,
                   const material_coefficient_t * ALIGNED(128) m,
                   const grid_t                 *              g );

// In v4_compute_div_b_err.c

void
v4_compute_div_b_err( field_t      * ALIGNED(128) f,
                      const grid_t *              g );

// In v4_clean_div_b.c

void
v4_clean_div_b( field_t      * ALIGNED(128) f,
                const grid_t *              g );

END_C_DECLS

#endif // _sfa_v4_private_h_
