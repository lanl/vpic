/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <vpic.hxx>

void vpic_simulation::initialize( int argc, char **argv ) {
  double tmp;
  species_t *sp;
  int rank;

  // Create an empty grid (creates the communicator too)
  grid = new_grid();
  if( grid==NULL ) ERROR(("Could not create the grid"));

  // Create a random number generator seeded with the rank
  rank = mp_rank(grid->mp);
  rng = mt_new_generator(rank);
  if( rng==NULL ) ERROR(("Could not create the random number generator"));

  // Call the user initialize the simulation
  user_initialization(argc,argv);

  // Do some consistency checks on user initialized fields
  if( rank==0 ) MESSAGE(("Checking interdomain synchronization"));
  tmp = synchronize_tang_e_norm_b( field, grid );
  if( rank==0 ) MESSAGE(("Error = %e (arb units)",tmp));

  if( rank==0 ) MESSAGE(("Checking magnetic field divergence"));
  compute_div_b_err( field, grid );
  tmp = compute_rms_div_b_err( field, grid );
  if( rank==0 ) MESSAGE(("RMS error = %e (charge/volume)",tmp));
  if( rank==0 ) MESSAGE(("trying clean_div_b"));
  clean_div_b( field, grid );
  if( rank==0 ) MESSAGE(("succeeded clean_div_b"));

  // Load fields not initialized by the user
  if( rank==0 ) MESSAGE(("Initializing radiation damping fields"));
  compute_curl_b( field, material_coefficient, grid );

  if( rank==0 ) MESSAGE(("Initializing bound charge density"));
  clear_rhof( field, grid );
  LIST_FOR_EACH(sp,species_list)
    accumulate_rho_p( field, sp->p, sp->np, grid );
  synchronize_rhof( field, grid );
  compute_rhob( field, material_coefficient, grid );

#if 1 // Internal sanity checks
  if( rank==0 ) MESSAGE(("Checking electric field divergence"));
  compute_div_e_err( field, material_coefficient, grid );
  tmp = compute_rms_div_e_err( field, grid );
  if( rank==0 ) MESSAGE(("RMS error = %e (charge/volume)",tmp));
  if( tmp>0 ) clean_div_e( field, material_coefficient, grid );

  if( rank==0 ) MESSAGE(("Double checking interdomain synchronization"));
  tmp = synchronize_tang_e_norm_b( field, grid );
  if( rank==0 ) MESSAGE(("Error = %e (arb units)",tmp));
#endif
    
  if( species_list!=NULL ) {
    if( rank==0 ) MESSAGE(("Uncentering particles"));
    load_interpolator( interpolator, field, grid );
  }
  LIST_FOR_EACH(sp, species_list)
    uncenter_p( sp->p, sp->np, sp->q_m, interpolator, grid );

  // FIXME!!!
  // start pipelines
  // This is a hack to get us through the RR assessment.  At some
  // point this will have to be re-worked to use overlays to allow
  // for multple accelerated implementations
  advance_p_initialize();

  if( rank==0 ) MESSAGE(("Initialization complete"));

  // Let the user to perform diagnostics on the initial condition
  // field(i,j,k).jfx, jfy, jfz will not be valid at this point.
  tmp = mp_time00(grid->mp);
  user_diagnostics();
  u_time += mp_time00(grid->mp) - tmp;
}
