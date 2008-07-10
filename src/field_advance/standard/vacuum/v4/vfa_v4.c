#define IN_vfa_v4
#include "vfa_v4_private.h"

field_advance_methods_t
_vacuum_v4_field_advance[1] = { {

  // Field array construction / destruction

  new_field,
  delete_field,

  // Materal coefficient construction / destruction

  vfa_new_material_coefficients,
  vfa_delete_material_coefficients,

  // Time stepping interfaces

  v4_advance_b,
  vfa_v4_advance_e,

  // Diagnostic interfaces

  vfa_energy_f,

  // Accumulator interfaces

  clear_jf,   synchronize_jf,
  clear_rhof, synchronize_rho,

  // Initialize interface

  vfa_compute_rhob,
  vfa_compute_curl_b, // FIXME: WRITE THIS EVENTUALLY

  // Shared face cleaning interface

  synchronize_tang_e_norm_b,
  
  // Electric field divergence cleaning interface

  vfa_compute_div_e_err,
  compute_rms_div_e_err,
  vfa_clean_div_e,

  // Magnetic field divergence cleaning interface

  v4_compute_div_b_err,
  compute_rms_div_b_err,
  v4_clean_div_b

} };
