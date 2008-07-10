#define IN_vfa
#include "vfa_private.h"

field_advance_methods_t
_vacuum_field_advance[1] = { {

  // Field array construction / destruction

  new_field,
  delete_field,

  // Materal coefficient construction / destruction

  vfa_new_material_coefficients,
  vfa_delete_material_coefficients,

  // Time stepping interfaces

  advance_b,
  vfa_advance_e,

  // Diagnostic interfaces

  vfa_energy_f,

  // Accumulator interfaces

  clear_jf,   synchronize_jf,
  clear_rhof, synchronize_rho,

  // Initialize interface

  vfa_compute_rhob,
  vfa_compute_curl_b,

  // Shared face cleaning interface

  synchronize_tang_e_norm_b,
  
  // Electric field divergence cleaning interface

  vfa_compute_div_e_err,
  compute_rms_div_e_err,
  vfa_clean_div_e,

  // Magnetic field divergence cleaning interface

  compute_div_b_err,
  compute_rms_div_b_err,
  clean_div_b

} };

/*****************************************************************************/

material_coefficient_t * ALIGNED(128)
vfa_new_material_coefficients( grid_t * g,
                               material_t * m_list ) {
  const material_t *m;

  // Check the input parameters
  if( g==NULL ) ERROR(("Invalid grid."));

  // FIXME: TCA damping should really be a property of the standard
  // field solver and not that grid
  if( g->damp!=0 )
    ERROR(( "Vacuum field advance does not support TCA radiation damping" ));

  // Run sanity checks on the material list
  // Namely, that there aren't any non-trivial materials.
  LIST_FOR_EACH(m,m_list)
    if( m->epsx!=1   || m->epsy!=1   || m->epsz!=1   ||
        m->mux!=1    || m->muy!=1    || m->muz!=1    ||
        m->sigmax!=0 || m->sigmay!=0 || m->sigmaz!=0 ||
        m->zetax!=0  || m->zetay!=0  || m->zetaz!=0  )
      ERROR(( "Material %s is not supported by vacuum (hint) field advance",
              m->name ));

  return NULL;
}

void
vfa_delete_material_coefficients( material_coefficient_t * ALIGNED(128) mc ) {
}

