/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version (data structures based on earlier
 *                    V4PIC versions)
 *
 */

#ifndef _species_advance_h_
#define _species_advance_h_

#include <algorithm>
#include <functional>
#include <string>

#include "../sf_interface/sf_interface.h"

//----------------------------------------------------------------------------//
// Choose between using AoSoA or AoS data layout for the particles.
//----------------------------------------------------------------------------//

#include "species_advance_aos.h"

//----------------------------------------------------------------------------//
// Declare methods.
//----------------------------------------------------------------------------//

BEGIN_C_DECLS

// In species_advance.cc

int
num_species( const species_t * sp_list );

void
delete_species_list( species_t * sp_list );

species_t *
find_species_id( species_id id,
                 species_t * sp_list );

species_t *
find_species_name( const char * name,
                   species_t * sp_list );


std::string
make_tracer_name_unique( const std::string prefix,
                         species_t* sp_list );

species_t *
append_species( species_t * sp,
                species_t ** sp_list );

species_t *
species( const char * name,
         float q,
         float m,
         size_t max_local_np,
         size_t max_local_nm,
         int sort_interval,
         int sort_out_of_place,
         grid_t * g
       );

enum class Tracertype { Copy, Move };
species_t *
tracerspecies_by_skip( species_t* parentspecies,
                       const float skip,
                       const Tracertype copyormove,
                       std::string tracername,
                       species_t* sp_list,
                       grid_t* grid
                     );
species_t *
tracerspecies_by_predicate( species_t* parentspecies,
                       std::function <bool (particle_t)> f,
                       const Tracertype copyormove,
                       std::string tracername,
                       species_t* sp_list,
                       grid_t* grid
                     );

// FIXME: TEMPORARY HACK UNTIL THIS SPECIES_ADVANCE KERNELS
// CAN BE CONSTRUCTED ANALOGOUS TO THE FIELD_ADVANCE KERNELS
// (THESE FUNCTIONS ARE NECESSARY FOR HIGHER LEVEL CODE)

// In sort_p.cc

void
sort_p( species_t * RESTRICT sp );

void
sort_p_pipeline( species_t * sp );

// In advance_p.cc

void
advance_p( species_t * RESTRICT sp,
           accumulator_array_t * RESTRICT aa,
           const interpolator_array_t * RESTRICT ia );

void
advance_p_pipeline( species_t * RESTRICT sp,
                    accumulator_array_t * RESTRICT aa,
                    const interpolator_array_t * RESTRICT ia );

// In center_p.cc

// This does a half advance field advance and a half Boris rotate on
// the particles.  As such particles with r at the time step and u
// half a step stale is moved second order accurate to have r and u on
// the time step.

void
center_p( species_t * RESTRICT sp,
          const interpolator_array_t * RESTRICT ia );

void
center_p_pipeline( species_t * RESTRICT sp,
                   const interpolator_array_t * RESTRICT ia );

// In uncenter_p.cc

// This is the inverse of center_p.  Thus, particles with r and u at
// the time step are adjusted to have r at the time step and u half a
// step stale.

void
uncenter_p( species_t * RESTRICT sp,
            const interpolator_array_t * RESTRICT ia );

void
uncenter_p_pipeline( species_t * RESTRICT sp,
                     const interpolator_array_t * RESTRICT ia );

// In energy.cc

// This computes the kinetic energy stored in the particles.  The
// calculation is done numerically robustly.  All nodes get the same
// result.

double
energy_p( const species_t * RESTRICT sp,
          const interpolator_array_t * RESTRICT ia );

double
energy_p_pipeline( const species_t * RESTRICT sp,
                   const interpolator_array_t * RESTRICT ia );

// In rho_p.cc

void
accumulate_rho_p( field_array_t * RESTRICT fa,
                  const species_t * RESTRICT sp );

void
accumulate_rhob( field_t * RESTRICT ALIGNED(128) f,
                 const particle_t * RESTRICT ALIGNED(32)  p,
                 const grid_t * RESTRICT g,
                 const float qsp );

// In hydro_p.cc

void
accumulate_hydro_p( hydro_array_t * RESTRICT ha,
                    const species_t * RESTRICT sp,
                    const interpolator_array_t * RESTRICT ia );

void
accumulate_hydro_p_pipeline( hydro_array_t * RESTRICT ha,
                             const species_t * RESTRICT sp,
                             const interpolator_array_t * RESTRICT ia );

// In move_p.cc

int
move_p( particle_t       * ALIGNED(128) p0,    // Particle array
        particle_mover_t * ALIGNED(16)  m,     // Particle mover to apply
        accumulator_t    * ALIGNED(128) a0,    // Accumulator to use
        const grid_t     *              g,     // Grid parameters
        const float                     qsp ); // Species particle charge

END_C_DECLS

#endif // _species_advance_h_
