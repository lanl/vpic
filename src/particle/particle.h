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

#ifndef _particle_h_
#define _particle_h_

#include <field.h> 

typedef struct particle {
  float dx, dy, dz; // Particle position in cell coordinates (on [-1,1])
  int32_t i;        // Index of cell containing the particle
  float ux, uy, uz; // Particle normalized momentum
  float q;          // Particle charge
} particle_t;

// WARNING: FUNCTIONS THAT USE A PARTICLE_MOVER ASSUME THAT EVERYBODY
// WHO USES THAT PARTICLE MOVER WILL BE GIVEN A POINTER TO PARTICLE 0

typedef struct particle_mover {
  float dispx, dispy, dispz; // Displacement of particle
  int32_t i;                 // Index of the particle to move
} particle_mover_t;

typedef struct particle_injector {
  float dx, dy, dz; int32_t i;
  float ux, uy, uz, q;
  float dispx, dispy, dispz;
  // FIXME: consider 16-byte aligning this!
} particle_injector_t;

BEGIN_C_DECLS

struct species; 

// In boundary_p.c

int
boundary_p( particle_mover_t * ALIGNED(16)  pm,
            int                             nm,
            int                             max_nm,
            particle_t       * ALIGNED(128) p0,
            int                             np,
            int                             max_np,
            field_t          * ALIGNED(16)  f,
            accumulator_t    * ALIGNED(128) a0,
            const grid_t     *              g,
            struct species   *              sp,
            mt_handle                       rng );

// In hydro_p.c

void
accumulate_hydro_p( hydro_t              * ALIGNED(16)  h0,
                    const particle_t     * ALIGNED(128) p0,
                    int                                 np,
                    float                               q_m,
                    const interpolator_t * ALIGNED(128) f0,
                    const grid_t         *              g );

// In move_p.c

int
inject_p( particle_t                * ALIGNED(128) p0, // Array to inject into
          int                                      np, // Where to inject
                                                       // Caller promises
                                                       // <max_np.
          particle_mover_t          * ALIGNED(16)  pm, // Free mover
          field_t                   * ALIGNED(16)  f,  // rhob accum
          accumulator_t             * ALIGNED(128) a0, // j accum
          const particle_injector_t *              pi,
          const grid_t              *              g );

// In sort_p.c

void
sort_p( struct species * sp, 
        const grid_t * g );

// In rho_p.c

void
accumulate_rho_p( field_t          * ALIGNED(16)  f,
                  const particle_t * ALIGNED(128) p0,
                  int                             np,
                  const grid_t     *              g );

// In particle_structors.c

particle_t * ALIGNED(128)
new_particle_array( int np );

particle_mover_t * ALIGNED(16)
new_particle_mover( int nm );

void
delete_particle_array( particle_t ** ALIGNED(128) p );

void
delete_particle_mover( particle_mover_t ** ALIGNED(16) pm );

// In advance_p.c

// Returns the number particle movers in use

int
advance_p( particle_t           * ALIGNED(128) p0,
           int                                 np,
           const float                         q_m,
           particle_mover_t     * ALIGNED(16)  pm,
           int                                 nm,
           accumulator_t        * ALIGNED(128) a0,
           const interpolator_t * ALIGNED(128) f0,
           const grid_t         *              g );

// In energy.c

// This computes the kinetic energy stored in the particles.  The
// calculation is done numerically robustly.  All nodes get the same
// result.

// FIXME: SHOULD THIS FUNCTION DO THE REDUCE OVER ALL NODES??

double
energy_p( const particle_t     * ALIGNED(128) p0,
          int                                 np,
          float                               q_m,
          const interpolator_t * ALIGNED(128) f0,
          const grid_t         *              g );

// In center_p.c

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

// In uncenter_p.c

// This is the inverse of center_p.  Thus, particles with r and u at
// the time step are adjusted to have r at the time step and u half a
// step stale.

void
uncenter_p( particle_t           * ALIGNED(128) p0,
            int                                 np,
            const float                         q_m,
            const interpolator_t * ALIGNED(128) f0,
            const grid_t         *              g );

// INTERNAL USE ONLY FUNCTIONS

// In move_p.c

// Note: changes to move_p likely need to be reflected in the SPU
// move_p implementation as well!

int
move_p( particle_t       * ALIGNED(128) p0, 
        particle_mover_t * ALIGNED(16)  m,
        accumulator_t    * ALIGNED(128) a0,
        const grid_t     *              g );

int
remove_p( particle_t   * ALIGNED(16)  r,
          particle_t   * ALIGNED(128) p0,
          int                         np,
          field_t      * ALIGNED(16)  f,
          const grid_t *              g );

END_C_DECLS

#endif
