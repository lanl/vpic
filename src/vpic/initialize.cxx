#include "vpic.hxx"

#define FA field_advance

void vpic_simulation::initialize( int argc, char **argv ) {
  double err;
  species_t * sp;
  int rank;

  // Create an empty grid (creates the communicator too)

  grid = new_grid();

  // Create a random number generator seeded with the rank

  rank = mp_rank( grid->mp );
  rng = new_mt_rng( rank );

  // Call the user initialize the simulation

  TIC { user_initialization( argc, argv ); } TOC( user_initialization, 1 );

  // Do some consistency checks on user initialized fields

  if( rank==0 ) MESSAGE(( "Checking interdomain synchronization" ));
  TIC { err = FA->method->synchronize_tang_e_norm_b( FA->f, FA->g ); } TOC( synchronize_tang_e_norm_b, 1 );
  if( rank==0 ) MESSAGE(( "Error = %e (arb units)", err ));

  if( rank==0 ) MESSAGE(( "Checking magnetic field divergence" ));
  TIC { FA->method->compute_div_b_err( FA->f, FA->g ); } TOC( compute_div_b_err, 1 );
  TIC { err = FA->method->compute_rms_div_b_err( FA->f, FA->g ); } TOC( compute_rms_div_b_err, 1 );
  if( rank==0 ) MESSAGE(( "RMS error = %e (charge/volume)", err ));
  TIC { FA->method->clean_div_b( FA->f, FA->g ); } TOC( clean_div_b, 1 );

  // Load fields not initialized by the user

  if( rank==0 ) MESSAGE(( "Initializing radiation damping fields" ));
  TIC { FA->method->compute_curl_b( FA->f, FA->m, FA->g ); } TOC( compute_curl_b, 1 );

  if( rank==0 ) MESSAGE(( "Initializing bound charge density" ));
  TIC { FA->method->clear_rhof( FA->f, FA->g ); } TOC( clear_rhof, 1 );
  LIST_FOR_EACH( sp, species_list ) TIC {
    accumulate_rho_p( FA->f, sp->p, sp->np, FA->g );
  } TOC( accumulate_rho_p, 1 );
  TIC { FA->method->synchronize_rho( FA->f, FA->g ); } TOC( synchronize_rho, 1 );
  TIC { FA->method->compute_rhob( FA->f, FA->m, FA->g ); } TOC( compute_rhob, 1 );

  // Internal sanity checks

  if( rank==0 ) MESSAGE(( "Checking electric field divergence" ));

  TIC { FA->method->compute_div_e_err( FA->f, FA->m, FA->g ); } TOC( compute_div_e_err, 1 );
  TIC { err = FA->method->compute_rms_div_e_err( FA->f, FA->g ); } TOC( compute_rms_div_e_err, 1 );
  if( rank==0 ) MESSAGE(( "RMS error = %e (charge/volume)", err ));
  TIC { FA->method->clean_div_e( FA->f, FA->m, FA->g ); } TOC( clean_div_e, 1 );

  if( rank==0 ) MESSAGE(( "Rechecking interdomain synchronization" ));
  TIC { err = FA->method->synchronize_tang_e_norm_b( FA->f, FA->g ); } TOC( synchronize_tang_e_norm_b, 1 );
  if( rank==0 ) MESSAGE(( "Error = %e (arb units)", err ));
    
  if( species_list!=NULL ) {
    if( rank==0 ) MESSAGE(( "Uncentering particles" ));
    TIC { load_interpolator( interpolator, FA->f, FA->g ); } TOC( load_interpolator, 1 );
  }
  LIST_FOR_EACH( sp, species_list ) TIC {
    uncenter_p( sp->p, sp->np, sp->q_m, interpolator, grid );
  } TOC( uncenter_p, 1 );

  if( rank==0 ) MESSAGE(( "Performing initial diagnostics" ));

  // Let the user to perform diagnostics on the initial condition
  // field(i,j,k).jfx, jfy, jfz will not be valid at this point.
  TIC { user_diagnostics(); } TOC( user_diagnostics, 1 );

  if( rank==0 ) MESSAGE(( "Initialization complete" ));
  update_profile( rank==0 ); // Let the user know how initialization went
}

