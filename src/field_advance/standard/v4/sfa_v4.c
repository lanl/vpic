#define IN_sfa_v4
#include "sfa_v4_private.h"

field_advance_methods_t
_standard_v4_field_advance[1] = { {

  // Field array construction / destruction

  new_field,
  delete_field,

  // Materal coefficient construction / destruction

  new_material_coefficients,
  delete_material_coefficients,

  // Time stepping interfaces

  v4_advance_b,
  v4_advance_e,

  // Diagnostic interfaces

  energy_f,

  // Accumulator interfaces

  clear_jf,   synchronize_jf,
  clear_rhof, synchronize_rho,

  // Initialize interface

  compute_rhob,
  v4_compute_curl_b,

  // Shared face cleaning interface

  synchronize_tang_e_norm_b,
  
  // Electric field divergence cleaning interface

  compute_div_e_err,
  compute_rms_div_e_err,
  clean_div_e,

  // Magnetic field divergence cleaning interface

  v4_compute_div_b_err,
  compute_rms_div_b_err,
  v4_clean_div_b

} };

