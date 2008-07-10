#ifndef _vfa_private_h_
#define _vfa_private_h_

// Vacuum field advance implementation

// FIMXE: IT WOULD BE BETTER IF VACUUM WAS THE BASE FIELD SOLVER
// AND STANDARD INHERITED FROM IT

#ifndef IN_vfa
#error "Do not include vfa_private.h; include field_advance.h"
#endif

#define IN_sfa
#include "../sfa_private.h"

// Note: See field methods for a detailed description of these functions

// Note: struct material_coefficient is _not_ used for
// vacuum_field_advance

// Note: TCA IS _NOT_ UPDATED UNDER VFA UNDER NORMAL OPERATION ONLY
// THE COMPUTE_CURL_B MECHANISM SETS IT (AND THAT HAPPENS DURING
// INITIALIZATION).

BEGIN_C_DECLS

// In vfa.c

material_coefficient_t * ALIGNED(128)
vfa_new_material_coefficients( grid_t * g,
                               material_t * m_list );

void
vfa_delete_material_coefficients( material_coefficient_t * ALIGNED(128) mc );

// In vfa_advance_e.cxx

void
vfa_advance_e( field_t                      * ALIGNED(128) f,
               const material_coefficient_t * ALIGNED(128) m, // Ignored
               const grid_t                 *              g );

// In vfa_energy_f.cxx

void
vfa_energy_f( double                       *              energy, // 6 elem
              const field_t                * ALIGNED(128) f,
              const material_coefficient_t * ALIGNED(128) m, // Ignored
              const grid_t                 *              g );

// In vfa_compute_curl_b.cxx

void
vfa_compute_curl_b( field_t                      * ALIGNED(128) f,
                    const material_coefficient_t * ALIGNED(128) m,
                    const grid_t                 *              g );

// In vfa_compute_rhob.cxx

void
vfa_compute_rhob( field_t                      * ALIGNED(128) f,
                  const material_coefficient_t * ALIGNED(128) m, // Ignored
                  const grid_t                 *              g ); 

// In vfa_compute_div_e_err.cxx

void
vfa_compute_div_e_err( field_t                      * ALIGNED(128) f,
                       const material_coefficient_t * ALIGNED(128) m, // Ignored
                       const grid_t                 *              g );

// In vfa_compute_rms_div_e_err.cxx

double
vfa_compute_rms_div_e_err( field_t      * ALIGNED(128) f,
                           const grid_t *              g );

// In clean_div_e.cxx

void
vfa_clean_div_e( field_t                      * ALIGNED(128) f,
                 const material_coefficient_t * ALIGNED(128) m, // Ignored
                 const grid_t                 *              g );

END_C_DECLS

#endif // _vfa_private_h_
