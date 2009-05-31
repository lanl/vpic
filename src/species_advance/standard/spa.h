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

// In sort_p.c

void
sort_p( species_t    * sp, 
        const grid_t * g );

// In boundary_p.cxx

void
accumulate_rhob( field_t          * ALIGNED(128) f0,  // Field data
                 const particle_t * ALIGNED(32)  p,   // Particle to remove
                 const grid_t     *              g ); // Grid params

void
boundary_p( species_t     *              sp_list, // Species params for use
            field_t       * ALIGNED(128) f0,      // For rhob accum and/or
                                                  // custom pbcs i.e. field
                                                  // emission
            accumulator_t * ALIGNED(128) a0,      // For j accum
            const grid_t  *              g,       // Local grid params
                                                  // in custom pbcs
            mt_rng_t      *              rng );   // Number of passes

// In move_p.cxx

// Note: changes to move_p likely need to be reflected in the SPU
// move_p implementation as well!

int
move_p( particle_t       * ALIGNED(128) p0,  // Particle array
        particle_mover_t * ALIGNED(16)  m,   // Particle mover to apply
        accumulator_t    * ALIGNED(128) a0,  // Accumulator to use
        const grid_t     *              g ); // Grid parameters

// In advance_p.cxx

int // Number of particles had a boundary interaction
advance_p( particle_t           * ALIGNED(128) p0,
           int                                 np,
           const float                         q_m,
           particle_mover_t     * ALIGNED(128) pm,
           int                                 max_nm,
           accumulator_t        * ALIGNED(128) a0,
           const interpolator_t * ALIGNED(128) f0,
           const grid_t         *              g );

// In center_p.cxx

// This does a half advance field advance and a half Boris rotate on
// the particles.  As such particles with r at the time step and u
// half a step stale is moved second order accurate to have r and u on
// the time step.

void
center_p( particle_t           * ALIGNED(128) p0,
          int                                 np,
          const float                         q_m,
          const interpolator_t * ALIGNED(128) f0,
          const grid_t         *              g );

// In uncenter_p.cxx

// This is the inverse of center_p.  Thus, particles with r and u at
// the time step are adjusted to have r at the time step and u half a
// step stale.

void
uncenter_p( particle_t           * ALIGNED(128) p0,
            int                                 np,
            const float                         q_m,
            const interpolator_t * ALIGNED(128) f0,
            const grid_t         *              g );

// In energy.cxx

// This computes the kinetic energy stored in the particles.  The
// calculation is done numerically robustly.  All nodes get the same
// result.

double
energy_p( const particle_t     * ALIGNED(128) p0,
          int                                 np,
          float                               q_m,
          const interpolator_t * ALIGNED(128) f0,
          const grid_t         *              g );

// In rho_p.cxx

void
accumulate_rho_p( field_t          * ALIGNED(128)  f,
                  const particle_t * ALIGNED(128) p0,
                  int                             np,
                  const grid_t     *              g );

// In hydro_p.c

void
accumulate_hydro_p( hydro_t              * ALIGNED(16)  h0,
                    const particle_t     * ALIGNED(128) p0,
                    int                                 np,
                    float                               q_m,
                    const interpolator_t * ALIGNED(128) f0,
                    const grid_t         *              g );

END_C_DECLS

#endif // _spa_h_
