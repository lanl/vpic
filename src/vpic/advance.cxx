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

int vpic_simulation::advance(void) {
  int rank;
  species_t *sp;
  emitter_t *emitter;
  double overhead, err;

  // Determine if we are done ... see note below why this is done here

  if( num_step>0 && step>=num_step ) return 0;

  overhead = mp_time00(grid->mp);

  rank = mp_rank(grid->mp);

  // At this point, fields are at E_0 and B_0 and the particle positions
  // are at r_0 and u_{-1/2}. Further the mover lists for the particles should
  // empty and all particles should be inside the local computational domain.
  // Advance the particle lists.

  if( species_list!=NULL ) clear_accumulators( accumulator, grid );
  LIST_FOR_EACH(sp,species_list) {
    if( sp->sort_interval>0 && step%sp->sort_interval==0 ) {
      if( rank==0 ) MESSAGE(("Performance sorting \"%s\".",sp->name));
      p_time += mp_time00(grid->mp) - overhead; overhead = mp_time00(grid->mp);
      sort_p( sp, grid );
      s_time += mp_time00(grid->mp) - overhead; overhead = mp_time00(grid->mp);
    } 
    sp->nm = advance_p( sp->p, sp->np, sp->q_m, sp->pm, sp->max_nm,
                        accumulator, interpolator, grid );
  }
  if( species_list!=NULL ) reduce_accumulators( accumulator, grid );

  p_time += mp_time00(grid->mp) - overhead; overhead = mp_time00(grid->mp);

  // BJA - for particle collisions (commented out to not collide with Kevin's
  //       and Ben's work on overlays). 
  // 
  // Collisions need to be done between sort and particle advance
  // Tally collision time as user time since collision models live in the 
  // input deck.   Collisions presently are implemented in user input 
  // decks; count their time against "user" time. 
  // 
  // FIXME:  A real interface for collisions? 
  // FIXME:  Give collisions their own time.  
  user_particle_collisions();

  u_time += mp_time00(grid->mp) - overhead; overhead = mp_time00(grid->mp);

  // Because the partial position push when injecting aged particles might
  // place those particles onto the guard list (boundary interaction) and
  // because advance_p requires an empty guard list, particle injection must
  // be done after advance_p and before guard list processing. Note:
  // user_particle_injection should be a stub if species_list is empty.

  LIST_FOR_EACH(emitter,emitter_list)
    (emitter->emission_model)( emitter, interpolator, field_advance->f, accumulator, field_advance->g, rng );
  user_particle_injection();

  u_time += mp_time00(grid->mp) - overhead; overhead = mp_time00(grid->mp);

  // At this point, most particle positions are at r_1 and u_{1/2}. Particles
  // that had boundary interactions are now on the guard list. Process the
  // guard lists. Particles that absorbed are added to rhob (using a corrected
  // local accumulation).

  for( int round=0; round<num_comm_round; round++ )
    boundary_p( species_list,
                field_advance->f, accumulator, field_advance->g, rng );

  LIST_FOR_EACH( sp, species_list ) {
    if( sp->nm!=0 && verbose )
      WARNING(( "Ignoring %i unprocessed %s movers (increase num_comm_round)",
                sp->nm, sp->name ));
    sp->nm = 0;
  }

  // At this point, all particle positions are at r_1 and u_{1/2}, the
  // guard lists are empty and the accumulators on each processor are current.
  // Convert the accumulators into currents.

  field_advance->method->clear_jf( field_advance->f, field_advance->g );
  if( species_list!=NULL )
    unload_accumulator( field_advance->f, accumulator, field_advance->g );
  field_advance->method->synchronize_jf( field_advance->f, field_advance->g );

  g_time += mp_time00(grid->mp) - overhead; overhead = mp_time00(grid->mp);

  // At this point, the particle currents are known at jf_{1/2}.
  // Let the user add their own current contributions. It is the users
  // responsibility to insure injected currents are consistent across domains.
  // It is also the users responsibility to update rhob according to
  // rhob_1 = rhob_0 + div juser_{1/2} (corrected local accumulation) if
  // the user wants electric field divergence cleaning to work.

  user_current_injection();

  u_time += mp_time00(grid->mp) - overhead; overhead = mp_time00(grid->mp);

  // Half advance the magnetic field from B_0 to B_{1/2}

  field_advance->method->advance_b( field_advance->f, field_advance->g, 0.5 );

  // Advance the electric field from E_0 to E_1

  field_advance->method->advance_e( field_advance->f, field_advance->m, field_advance->g );

  f_time += mp_time00(grid->mp) - overhead; overhead = mp_time00(grid->mp);

  // Let the user add their own contributions to the electric field. It is the
  // users responsibility to insure injected electric fields are consistent
  // across domains.

  user_field_injection();

  u_time += mp_time00(grid->mp) - overhead; overhead = mp_time00(grid->mp);

  // Half advance the magnetic field from B_{1/2} to B_1

  field_advance->method->advance_b( field_advance->f, field_advance->g, 0.5 );

  // Divergence clean e

  if( clean_div_e_interval>0 && step%clean_div_e_interval==0 ) {
    if( rank==0 ) MESSAGE(("Divergence cleaning electric field"));
    field_advance->method->clear_rhof( field_advance->f, field_advance->g );
    LIST_FOR_EACH(sp,species_list)
      accumulate_rho_p( field_advance->f, sp->p, sp->np, field_advance->g );
    field_advance->method->synchronize_rho( field_advance->f, field_advance->g );
    field_advance->method->compute_div_e_err( field_advance->f, field_advance->m, field_advance->g );
    err = field_advance->method->compute_rms_div_e_err( field_advance->f, field_advance->g );
    if( rank==0 ) MESSAGE(("Initial rms error = %e (charge/volume)",err));
    if( err>0 ) {
      field_advance->method->clean_div_e( field_advance->f, field_advance->m, field_advance->g );
      field_advance->method->compute_div_e_err( field_advance->f, field_advance->m, field_advance->g );
      err = field_advance->method->compute_rms_div_e_err( field_advance->f, field_advance->g );
      if( rank==0 ) MESSAGE(("Cleaned rms error = %e (charge/volume)",err));
      if( err>0 ) field_advance->method->clean_div_e( field_advance->f, field_advance->m, field_advance->g );
    }
  }

  // Divergence clean b

  if( clean_div_b_interval>0 && step%clean_div_b_interval==0 ) {
    if( rank==0 ) MESSAGE(("Divergence cleaning magnetic field"));
    field_advance->method->compute_div_b_err( field_advance->f, field_advance->g );
    err = field_advance->method->compute_rms_div_b_err( field_advance->f, field_advance->g );
    if( rank==0 ) MESSAGE(("Initial rms error = %e (charge/volume)",err));
    if( err>0 ) {
      field_advance->method->clean_div_b( field_advance->f, field_advance->g );
      field_advance->method->compute_div_b_err( field_advance->f, field_advance->g );
      err = field_advance->method->compute_rms_div_b_err( field_advance->f, field_advance->g );
      if( rank==0 ) MESSAGE(("Cleaned rms error = %e (charge/volume)",err));
      if( err>0 ) field_advance->method->clean_div_b( field_advance->f, field_advance->g );
    }
  }

  // Synchronize the shared faces

  if( sync_shared_interval>0 && step%sync_shared_interval==0 ) {
    if( rank==0 ) MESSAGE(("Synchronizing shared tang e, norm b, rho_b"));
    err = field_advance->method->synchronize_tang_e_norm_b( field_advance->f, field_advance->g );
    if( rank==0 )
      MESSAGE(("Domain desynchronization error = %e (arb units)",err));
  }

  // Fields are updated ... load the interpolator for next time step and
  // particle diagnostics in user_diagnostics if there are any particle
  // species to worry about

  if( species_list!=NULL ) load_interpolator( interpolator, field_advance->f, field_advance->g );

  f_time += mp_time00(grid->mp) - overhead;

  step++;

  // Print out status

  if( status_interval>0 && step%status_interval==0 ) {
    if(rank==0)
      MESSAGE(("Completed step %i of %i (p=%.2e,s=%.2e,g=%.2e,f=%.2e,u=%.2e)",
               step, num_step, p_time-s_time, s_time, g_time, f_time, u_time));
    p_time = s_time = g_time = f_time = u_time = 0;
  }

  overhead = mp_time00(grid->mp);

  // Let the user compute diagnostics

  user_diagnostics();

  u_time += mp_time00(grid->mp) - overhead;

  // "return step!=num_step" is more intuitive. But if a restart dump,
  // saved in the call to user_diagnostics() above, is done on the final step
  // (silly but it might happen), the test will be skipped on the restart. We
  // return true here so that the first call to advance after a restart dump
  // will act properly for this edge case.

  return 1;
}
