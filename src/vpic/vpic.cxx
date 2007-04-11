/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Heavily revised and extended from earlier V4PIC versions
 *
 */

#include <vpic.hxx>

vpic_simulation::vpic_simulation() {
  step = 0;
  num_step = 0;
  status_interval = 0;
  clean_div_e_interval = 0;
  clean_div_b_interval = 0;
  sync_shared_interval = 0;

  quota=0;
  restart_interval=0;
  hydro_interval=0;
  field_interval=0;
  particle_interval=0;

  rng = NULL;
  material_list = NULL;
  material_coefficient = NULL;
  grid = NULL;
  field = NULL;
  interpolator = NULL;
  accumulator = NULL;
  hydro = NULL;
  species_list = NULL;
  species_lookup = NULL;
  emitter_list = NULL;
  p_time = 0;
  g_time = 0;
  f_time = 0;
  u_time = 0;

  /* Allow run-time modification of variables */ 
  quota=11;

  memset( user_global, 0, USER_GLOBAL_SIZE );
}

vpic_simulation::~vpic_simulation() {
  step = 0;
  num_step = 0;
  status_interval = 0;
  clean_div_e_interval = 0;
  clean_div_b_interval = 0;
  sync_shared_interval = 0;

  quota=0;
  restart_interval=0;
  hydro_interval=0;
  field_interval=0;
  particle_interval=0;

  mt_delete_generator( &rng );
  delete_material_list( &material_list );
  delete_material_coefficients( &material_coefficient );
  delete_grid( &grid );
  delete_field( &field );
  delete_interpolator( &interpolator );
  delete_accumulators( &accumulator );
  delete_hydro( &hydro );
  delete_species_list( &species_list );
  delete_species_lookup( &species_lookup );
  delete_emitter_list( &emitter_list );
  p_time = 0;
  g_time = 0;
  f_time = 0;
  u_time = 0;
  memset( user_global, 0, USER_GLOBAL_SIZE );
}

