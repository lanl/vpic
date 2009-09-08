/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version (data structures were revised from
 *                    earlier V4PIC versions)
 *
 */

// FIXME: THIS FILE SHOULD EVENTUALLY BE DEFLATED AND REMOVED

#ifndef _spa_h_
#define _spa_h_

#include "species_advance.h"

BEGIN_C_DECLS

// In move_p.cxx

// Note: changes to move_p likely need to be reflected in the SPU
// move_p implementation as well!

int
move_p( particle_t       * ALIGNED(128) p0,    // Particle array
        particle_mover_t * ALIGNED(16)  m,     // Particle mover to apply
        accumulator_t    * ALIGNED(128) a0,    // Accumulator to use
        const grid_t     *              g,     // Grid parameters
        const float                     qsp ); // Species particle charge

// In boundary_p.cxx

void
accumulate_rhob( field_t          * RESTRICT ALIGNED(128) f,     // Field data
                 const particle_t * RESTRICT ALIGNED(32)  p,     // Particle to remove
                 const grid_t     * RESTRICT              g,     // Grid
                 const float                              qsp ); // Species charge

// FIXME: THIS FUNCTION BELONGS IN BOUNDARY
struct particle_bc;
typedef struct particle_bc particle_bc_t;
void
boundary_p( species_t           * RESTRICT sp_list,  // Species to handle
            particle_bc_t       * RESTRICT pbc_list, // Particle boundary conds
            field_array_t       * RESTRICT fa,       // For rhob accum and/or
                                                     // custom pbcs i.e. field
                                                     // emission
            accumulator_array_t * RESTRICT aa );     // For j accum

// In sort_p.c

void
sort_p( species_t * RESTRICT sp );

// In advance_p.cxx

void
advance_p( /**/  species_t            * RESTRICT sp,
           /**/  accumulator_array_t  * RESTRICT aa,
           const interpolator_array_t * RESTRICT ia );

// In center_p.cxx

// This does a half advance field advance and a half Boris rotate on
// the particles.  As such particles with r at the time step and u
// half a step stale is moved second order accurate to have r and u on
// the time step.

void
center_p( /**/  species_t            * RESTRICT sp,
          const interpolator_array_t * RESTRICT ia );

// In uncenter_p.cxx

// This is the inverse of center_p.  Thus, particles with r and u at
// the time step are adjusted to have r at the time step and u half a
// step stale.

void
uncenter_p( /**/  species_t            * RESTRICT sp,
            const interpolator_array_t * RESTRICT ia );

// In energy.cxx

// This computes the kinetic energy stored in the particles.  The
// calculation is done numerically robustly.  All nodes get the same
// result.

double
energy_p( const species_t            * RESTRICT sp,
          const interpolator_array_t * RESTRICT ia );

// In rho_p.cxx

void
accumulate_rho_p( /**/  field_array_t * RESTRICT fa,
                  const species_t     * RESTRICT sp );

// In hydro_p.c

void
accumulate_hydro_p( /**/  hydro_array_t        * RESTRICT ha,
                    const species_t            * RESTRICT sp,
                    const interpolator_array_t * RESTRICT ia );

END_C_DECLS

#endif // _spa_h_
