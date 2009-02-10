/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "vpic.hxx"

void vpic_simulation::initialize( int argc, char **argv ) {
  double tmp;
  species_t *sp;
  int rank;

  // Create an empty grid (creates the communicator too)
  grid = new_grid();

  // Create a random number generator seeded with the rank
  rank = mp_rank(grid->mp);
  rng = new_mt_rng( rank );

  // Call the user initialize the simulation
  user_initialization(argc,argv);

  // Do some consistency checks on user initialized fields
  if( rank==0 ) MESSAGE(("Checking interdomain synchronization"));
  tmp = field_advance->method->synchronize_tang_e_norm_b( field_advance->f, field_advance->g );
  if( rank==0 ) MESSAGE(("Error = %e (arb units)",tmp));

  if( rank==0 ) MESSAGE(("Checking magnetic field divergence"));
  field_advance->method->compute_div_b_err( field_advance->f, field_advance->g );
  tmp = field_advance->method->compute_rms_div_b_err( field_advance->f, field_advance->g );
  if( rank==0 ) MESSAGE(("RMS error = %e (charge/volume)",tmp));
  if( rank==0 ) MESSAGE(("trying clean_div_b"));
  field_advance->method->clean_div_b( field_advance->f, field_advance->g );
  if( rank==0 ) MESSAGE(("succeeded clean_div_b"));

  // Load fields not initialized by the user
  if( rank==0 ) MESSAGE(("Initializing radiation damping fields"));
  field_advance->method->compute_curl_b( field_advance->f, field_advance->m, field_advance->g );

  if( rank==0 ) MESSAGE(("Initializing bound charge density"));
  field_advance->method->clear_rhof( field_advance->f, field_advance->g );
  LIST_FOR_EACH(sp,species_list)
    accumulate_rho_p( field_advance->f, sp->p, sp->np, field_advance->g );
  field_advance->method->synchronize_rho( field_advance->f, field_advance->g );
  field_advance->method->compute_rhob( field_advance->f, field_advance->m, field_advance->g );

# if 1 // Internal sanity checks
  if( rank==0 ) MESSAGE(("Checking electric field divergence"));
  field_advance->method->compute_div_e_err( field_advance->f, field_advance->m, field_advance->g );
  tmp = field_advance->method->compute_rms_div_e_err( field_advance->f, field_advance->g );
  if( rank==0 ) MESSAGE(("RMS error = %e (charge/volume)",tmp));
  if( tmp>0 ) field_advance->method->clean_div_e( field_advance->f, field_advance->m, field_advance->g );

  if( rank==0 ) MESSAGE(("Double checking interdomain synchronization"));
  tmp = field_advance->method->synchronize_tang_e_norm_b( field_advance->f, field_advance->g );
  if( rank==0 ) MESSAGE(("Error = %e (arb units)",tmp));
# endif
    
  if( species_list!=NULL ) {
    if( rank==0 ) MESSAGE(("Uncentering particles"));
    load_interpolator( interpolator, field_advance->f, field_advance->g );
  }
  LIST_FOR_EACH(sp, species_list)
    uncenter_p( sp->p, sp->np, sp->q_m, interpolator, grid );

  if( rank==0 ) MESSAGE(("Initialization complete"));

  // Let the user to perform diagnostics on the initial condition
  // field(i,j,k).jfx, jfy, jfz will not be valid at this point.
  tmp = mp_time00(grid->mp);
  user_diagnostics();
  u_time += mp_time00(grid->mp) - tmp;
}
