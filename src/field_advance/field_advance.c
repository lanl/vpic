#define IN_field_advance
#include <field_advance_private.h>

void
delete_field_array( field_array_t * fa ) {
  if( !fa ) return;
  fa->kernel->delete_fa( fa );
}

void
checkpt_field_advance_kernels( const field_advance_kernels_t * kernel ) {
  CHECKPT_SYM( kernel->delete_fa                 );
  CHECKPT_SYM( kernel->advance_b                 );
  CHECKPT_SYM( kernel->advance_e                 );
  CHECKPT_SYM( kernel->energy_f                  );
  CHECKPT_SYM( kernel->clear_jf                  );
  CHECKPT_SYM( kernel->synchronize_jf            );
  CHECKPT_SYM( kernel->clear_rhof                );
  CHECKPT_SYM( kernel->synchronize_rho           );
  CHECKPT_SYM( kernel->compute_rhob              );
  CHECKPT_SYM( kernel->compute_curl_b            );
  CHECKPT_SYM( kernel->synchronize_tang_e_norm_b );
  CHECKPT_SYM( kernel->compute_div_e_err         );
  CHECKPT_SYM( kernel->compute_rms_div_e_err     );
  CHECKPT_SYM( kernel->clean_div_e               );
  CHECKPT_SYM( kernel->compute_div_b_err         );
  CHECKPT_SYM( kernel->compute_rms_div_b_err     );
  CHECKPT_SYM( kernel->clean_div_b               );
}

void
restore_field_advance_kernels( field_advance_kernels_t * kernel ) {
  RESTORE_SYM( kernel->delete_fa                 );
  RESTORE_SYM( kernel->advance_b                 );
  RESTORE_SYM( kernel->advance_e                 );
  RESTORE_SYM( kernel->energy_f                  );
  RESTORE_SYM( kernel->clear_jf                  );
  RESTORE_SYM( kernel->synchronize_jf            );
  RESTORE_SYM( kernel->clear_rhof                );
  RESTORE_SYM( kernel->synchronize_rho           );
  RESTORE_SYM( kernel->compute_rhob              );
  RESTORE_SYM( kernel->compute_curl_b            );
  RESTORE_SYM( kernel->synchronize_tang_e_norm_b );
  RESTORE_SYM( kernel->compute_div_e_err         );
  RESTORE_SYM( kernel->compute_rms_div_e_err     );
  RESTORE_SYM( kernel->clean_div_e               );
  RESTORE_SYM( kernel->compute_div_b_err         );
  RESTORE_SYM( kernel->compute_rms_div_b_err     );
  RESTORE_SYM( kernel->clean_div_b               );
}

