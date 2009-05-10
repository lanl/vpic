/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Heavily revised and extended from earlier V4PIC versions
 *
 */

#include "vpic.hxx"

#define FA field_advance

int vpic_simulation::advance(void) {
  int rank;
  species_t *sp;
  emitter_t *emitter;
  double err;

  // Determine if we are done ... see note below why this is done here

  if( num_step>0 && step>=num_step ) return 0;

  rank = mp_rank( grid->mp );

  // At this point, fields are at E_0 and B_0 and the particle positions
  // are at r_0 and u_{-1/2}. Further the mover lists for the particles should
  // empty and all particles should be inside the local computational domain.
  // Advance the particle lists.

  if( species_list!=NULL ) TIC { clear_accumulators( accumulator, grid ); } TOC( clear_accumulators, 1 );

  LIST_FOR_EACH( sp, species_list )
    if( (sp->sort_interval>0) && ((step % sp->sort_interval)==0) ) {
      if( rank==0 ) MESSAGE(( "Performance sorting \"%s\".", sp->name ));
      TIC { sort_p( sp, grid ); } TOC( sort_p, 1 );
    } 

  TIC { user_particle_collisions(); } TOC( user_particle_collisions, 1 );

  LIST_FOR_EACH( sp, species_list ) TIC {
    sp->nm = advance_p( sp->p, sp->np, sp->q_m, sp->pm, sp->max_nm, accumulator, interpolator, grid );
  } TOC( advance_p, 1 );

  if( species_list!=NULL ) TIC { reduce_accumulators( accumulator, grid ); } TOC( reduce_accumulators, 1 );

  // Because the partial position push when injecting aged particles might
  // place those particles onto the guard list (boundary interaction) and
  // because advance_p requires an empty guard list, particle injection must
  // be done after advance_p and before guard list processing. Note:
  // user_particle_injection should be a stub if species_list is empty.

  LIST_FOR_EACH( emitter, emitter_list ) TIC {
    (emitter->emission_model)( emitter, interpolator, FA->f, accumulator, FA->g, rng );
  } TOC( emission_model, 1 );

  TIC { user_particle_injection(); } TOC( user_particle_injection, 1 );

  // At this point, most particle positions are at r_1 and u_{1/2}. Particles
  // that had boundary interactions are now on the guard list. Process the
  // guard lists. Particles that absorbed are added to rhob (using a corrected
  // local accumulation).

  TIC {
    for( int round=0; round<num_comm_round; round++ )
      boundary_p( species_list, FA->f, accumulator, FA->g, rng );
  } TOC( boundary_p, num_comm_round );

  LIST_FOR_EACH( sp, species_list ) {
    if( sp->nm!=0 && !verbose ) WARNING(( "Ignoring %i unprocessed %s movers (increase num_comm_round)",
                                          sp->nm, sp->name ));
    sp->nm = 0;
  }

  // At this point, all particle positions are at r_1 and u_{1/2}, the
  // guard lists are empty and the accumulators on each processor are current.
  // Convert the accumulators into currents.

  TIC { FA->method->clear_jf( FA->f, FA->g ); } TOC( clear_jf, 1 );
  if( species_list!=NULL ) TIC { unload_accumulator( FA->f, accumulator, FA->g ); } TOC( unload_accumulator, 1 );
  TIC { FA->method->synchronize_jf( FA->f, FA->g ); } TOC( synchronize_jf, 1 );

  // At this point, the particle currents are known at jf_{1/2}.
  // Let the user add their own current contributions. It is the users
  // responsibility to insure injected currents are consistent across domains.
  // It is also the users responsibility to update rhob according to
  // rhob_1 = rhob_0 + div juser_{1/2} (corrected local accumulation) if
  // the user wants electric field divergence cleaning to work.

  TIC { user_current_injection(); } TOC( user_current_injection, 1 );

  // Half advance the magnetic field from B_0 to B_{1/2}

  TIC { FA->method->advance_b( FA->f, FA->g, 0.5 ); } TOC( advance_b, 1 );

  // Advance the electric field from E_0 to E_1

  TIC { FA->method->advance_e( FA->f, FA->m, FA->g ); } TOC( advance_e, 1 );

  // Let the user add their own contributions to the electric field. It is the
  // users responsibility to insure injected electric fields are consistent
  // across domains.

  TIC { user_field_injection(); } TOC( user_field_injection, 1 );

  // Half advance the magnetic field from B_{1/2} to B_1

  TIC { FA->method->advance_b( FA->f, FA->g, 0.5 ); } TOC( advance_b, 1 );

  // Divergence clean e

  if( (clean_div_e_interval>0) && ((step % clean_div_e_interval)==0) ) {
    if( rank==0 ) MESSAGE(( "Divergence cleaning electric field" ));

    TIC { FA->method->clear_rhof( FA->f, FA->g ); } TOC( clear_rhof,1 );
    LIST_FOR_EACH( sp, species_list ) TIC { accumulate_rho_p( FA->f, sp->p, sp->np, FA->g ); } TOC( accumulate_rho_p, 1 );
    TIC { FA->method->synchronize_rho( FA->f, FA->g ); } TOC( synchronize_rho, 1 );

    for( int round=0; round<num_div_e_round; round++ ) {
      TIC { FA->method->compute_div_e_err( FA->f, FA->m, FA->g ); } TOC( compute_div_e_err, 1 );
      if( round==0 || round==num_div_e_round-1 ) {
        TIC { err = FA->method->compute_rms_div_e_err( FA->f, FA->g ); } TOC( compute_rms_div_e_err, 1 );
        if( rank==0 ) MESSAGE(( "%s rms error = %e (charge/volume)", round==0 ? "Initial" : "Cleaned", err ));
      }
      TIC { FA->method->clean_div_e( FA->f, FA->m, FA->g ); } TOC( clean_div_e, 1 );
    }
  }

  // Divergence clean b

  if( (clean_div_b_interval>0) && ((step % clean_div_b_interval)==0) ) {
    if( rank==0 ) MESSAGE(( "Divergence cleaning magnetic field" ));

    for( int round=0; round<num_div_b_round; round++ ) {
      TIC { FA->method->compute_div_b_err( FA->f, FA->g ); } TOC( compute_div_b_err, 1 );
      if( round==0 || round==num_div_b_round-1 ) {
        TIC { err = FA->method->compute_rms_div_b_err( FA->f, FA->g ); } TOC( compute_rms_div_b_err, 1 );
        if( rank==0 ) MESSAGE(( "%s rms error = %e (charge/volume)", round==0 ? "Initial" : "Cleaned", err ));
      }
      TIC { FA->method->clean_div_b( FA->f, FA->g ); } TOC( clean_div_b, 1 );
    }
  }

  // Synchronize the shared faces

  if( (sync_shared_interval>0) && ((step % sync_shared_interval)==0) ) {
    if( rank==0 ) MESSAGE(( "Synchronizing shared tang e, norm b, rho_b" ));
    TIC { err = FA->method->synchronize_tang_e_norm_b( FA->f, FA->g ); } TOC( synchronize_tang_e_norm_b, 1 );
    if( rank==0 ) MESSAGE(( "Domain desynchronization error = %e (arb units)", err ));
  }

  // Fields are updated ... load the interpolator for next time step and
  // particle diagnostics in user_diagnostics if there are any particle
  // species to worry about

  if( species_list!=NULL ) TIC { load_interpolator( interpolator, FA->f, FA->g ); } TOC( load_interpolator, 1 );

  step++;

  // Print out status

  if( (status_interval>0) && ((step % status_interval)==0) ) {
    if( rank==0 ) MESSAGE(( "Completed step %i of %i", step, num_step ));
    update_profile( rank==0 );
  }

  // Let the user compute diagnostics

  TIC { user_diagnostics(); } TOC( user_diagnostics, 1 );

  // "return step!=num_step" is more intuitive. But if a restart dump,
  // saved in the call to user_diagnostics() above, is done on the final step
  // (silly but it might happen), the test will be skipped on the restart. We
  // return true here so that the first call to advance after a restart dump
  // will act properly for this edge case.

  return 1;
}
