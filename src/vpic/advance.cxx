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
#include <string.h>

int vpic_simulation::advance(void) {
  int rank;
  species_t *sp;
  emitter_t *emitter;
  double overhead, err;

  overhead = mp_time00(grid->mp);
  // Synchronize all the processors at this timestep and determine if we are
  // done ... see note at end as to why this is done here
  mp_barrier(grid->mp); // FIXME: THIS SHOULD GO!!!!
  if( num_step>0 && step>=num_step ) return 0;
  rank = mp_rank(grid->mp);
  g_time += mp_time00(grid->mp) - overhead;

  overhead = mp_time00(grid->mp);
  // At this point, fields are at E_0 and B_0 and the particle positions
  // are at r_0 and u_{-1/2}. Further the mover lists for the particles should
  // empty and all particles should be inside the local computational domain.
  // Advance the particle lists.
  if( species_list!=NULL ) clear_accumulators( accumulator, grid );
  LIST_FOR_EACH(sp,species_list) {
    if( sp->sort_interval>0 && step%sp->sort_interval==0 ) {
      if( rank==0 ) MESSAGE(("Performance sorting \"%s\".",sp->name));
      sort_p( sp, grid );
    } 
    sp->nm = advance_p( sp->p, sp->np, sp->q_m, sp->pm, sp->max_nm,
                        accumulator, interpolator, grid );
  }
  if( species_list!=NULL ) reduce_accumulators( accumulator, grid );
  p_time += mp_time00(grid->mp) - overhead;

  overhead = mp_time00(grid->mp);
  // Because the partial position push when injecting aged particles might
  // place those particles onto the guard list (boundary interaction) and
  // because advance_p requires an empty guard list, particle injection must
  // be done after advance_p and before guard list processing. Note:
  // user_particle_injection should be a stub if species_list is empty.
  LIST_FOR_EACH(emitter,emitter_list)
    (emitter->emission_model)( emitter, interpolator, field, accumulator, grid, rng );
  user_particle_injection();
  u_time += mp_time00(grid->mp) - overhead;

  overhead = mp_time00(grid->mp);
  // At this point, most particle positions are at r_1 and u_{1/2}. Particles
  // that had boundary interactions are now on the guard list. Process the
  // guard lists. Particles that absorbed are added to rhob (using a corrected
  // local accumulation).
  LIST_FOR_EACH(sp,species_list) 
    sp->np = boundary_p( sp->pm, sp->nm, sp->max_nm,
                         sp->p,  sp->np, sp->max_np,
                         field, accumulator, grid, sp, rng );
  // At this point, all particle positions are at r_1 and u_{1/2}, the
  // guard lists are empty and the accumulators on each processor are current.
  // Convert the accumulators into currents.
  clear_jf( field, grid );
  if( species_list!=NULL )
    unload_accumulator( field, accumulator, grid );
  synchronize_jf( field, grid );
  g_time += mp_time00(grid->mp) - overhead;

  overhead = mp_time00(grid->mp);
  // At this point, the particle currents are known at jf_{1/2}.
  // Let the user add their own current contributions. It is the users
  // responsibility to insure injected currents are consistent across domains.
  // It is also the users responsibility to update rhob according to
  // rhob_1 = rhob_0 + div juser_{1/2} (corrected local accumulation) if
  // the user wants electric field divergence cleaning to work.
  user_current_injection();
  u_time += mp_time00(grid->mp) - overhead;

  overhead = mp_time00(grid->mp);
  // Half advance the magnetic field from B_0 to B_{1/2}
  advance_b( field, grid, 0.5 );
  // Advance the electric field from E_0 to E_1
  advance_e( field, material_coefficient, grid );
  f_time += mp_time00(grid->mp) - overhead;

  overhead = mp_time00(grid->mp);
  // Let the user add their own contributions to the electric field. It is the
  // users responsibility to insure injected electric fields are consistent
  // across domains.
  user_field_injection();
  u_time += mp_time00(grid->mp) - overhead;

  overhead = mp_time00(grid->mp);
  // Half advance the magnetic field from B_{1/2} to B_1
  advance_b( field, grid, 0.5 );
  // Divergence clean e
  if( clean_div_e_interval>0 && step%clean_div_e_interval==0 ) {
    if( rank==0 ) MESSAGE(("Divergence cleaning electric field"));
    clear_rhof( field, grid );
    LIST_FOR_EACH(sp,species_list)
      accumulate_rho_p( field, sp->p, sp->np, grid );
    synchronize_rhof( field, grid );
    synchronize_rhob( field, grid );
    compute_div_e_err( field, material_coefficient, grid );
    err = compute_rms_div_e_err( field, grid );
    if( rank==0 ) MESSAGE(("Initial rms error = %e (charge/volume)",err));
    if( err>0 ) {
      clean_div_e( field, material_coefficient, grid );
      compute_div_e_err( field, material_coefficient, grid );
      err = compute_rms_div_e_err( field, grid );
      if( rank==0 ) MESSAGE(("Cleaned rms error = %e (charge/volume)",err));
      if( err>0 ) clean_div_e( field, material_coefficient, grid );
    }
  }
  // Divergence clean b
  if( clean_div_b_interval>0 && step%clean_div_b_interval==0 ) {
    if( rank==0 ) MESSAGE(("Divergence cleaning magnetic field"));
    compute_div_b_err( field, grid );
    err = compute_rms_div_b_err( field, grid );
    if( rank==0 ) MESSAGE(("Initial rms error = %e (charge/volume)",err));
    if( err>0 ) {
      clean_div_b( field, grid );
      compute_div_b_err( field, grid );
      err = compute_rms_div_b_err( field, grid );
      if( rank==0 ) MESSAGE(("Cleaned rms error = %e (charge/volume)",err));
      if( err>0 ) clean_div_b( field, grid );
    }
  }
  // Synchronize the shared faces
  if( sync_shared_interval>0 && step%sync_shared_interval==0 ) {
    if( rank==0 ) MESSAGE(("Synchronizing shared tang e, norm b, rho_b"));
    synchronize_rhob( field, grid );
    err = synchronize_tang_e_norm_b( field, grid );
    if( rank==0 )
      MESSAGE(("Domain desynchronization error = %e (arb units)",err));
  }
  // Fields are updated ... load the interpolator for next time step and
  // particle diagnostics in user_diagnostics if there are any particle
  // species to worry about
  if( species_list!=NULL ) load_interpolator( interpolator, field, grid );
  step++;
  f_time += mp_time00(grid->mp) - overhead;

  // Print out status
  if( status_interval>0 && step%status_interval==0 ) {
    if(rank==0)
      MESSAGE(("Completed step %i of %i (p=%.2e,g=%.2e,f=%.2e,u=%.2e)",
               step, num_step, p_time, g_time, f_time, u_time));
    p_time = g_time = f_time = u_time = 0;
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
