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

#include "../sf_interface/sf_interface.h"

typedef int32_t species_id; // Must be 32-bit wide for particle_injector_t

// FIXME: Eventually particle_t (definitely) and their other formats
// (maybe) should be opaque and specific to a particular
// species_advance implementation

typedef struct particle {
  float dx, dy, dz; // Particle position in cell coordinates (on [-1,1])
  int32_t i;        // Voxel containing the particle.  Note that
  /**/              // particled awaiting processing by boundary_p
  /**/              // have actually set this to 8*voxel + face where
  /**/              // face is the index of the face they interacted
  /**/              // with (on 0:5).  This limits the local number of
  /**/              // voxels to 2^28 but emitter handling already
  /**/              // has a stricter limit on this (2^26).
  float ux, uy, uz; // Particle normalized momentum
  float w;          // Particle weight (number of physical particles)
} particle_t;

// WARNING: FUNCTIONS THAT USE A PARTICLE_MOVER ASSUME THAT EVERYBODY
// WHO USES THAT PARTICLE MOVER WILL HAVE ACCESS TO PARTICLE ARRAY

typedef struct particle_mover {
  float dispx, dispy, dispz; // Displacement of particle
  int32_t i;                 // Index of the particle to move
} particle_mover_t;

// NOTE: THE LAYOUT OF A PARTICLE_INJECTOR _MUST_ BE COMPATIBLE WITH
// THE CONCATENATION OF A PARTICLE_T AND A PARTICLE_MOVER!

typedef struct particle_injector {
  float dx, dy, dz;          // Particle position in cell coords (on [-1,1])
  int32_t i;                 // Index of cell containing the particle
  float ux, uy, uz;          // Particle normalized momentum
  float w;                   // Particle weight (number of physical particles)
  float dispx, dispy, dispz; // Displacement of particle
  species_id sp_id;          // Species of particle
} particle_injector_t;

typedef struct species {
  char * name;                        // Species name
  float q;                            // Species particle charge
  float m;                            // Species particle rest mass

  int np, max_np;                     // Number and max local particles
  particle_t * ALIGNED(128) p;        // Array of particles for the species

  int nm, max_nm;                     // Number and max local movers in use
  particle_mover_t * ALIGNED(128) pm; // Particle movers

  int64_t last_sorted;                // Step when the particles were last
                                      // sorted.
  int sort_interval;                  // How often to sort the species
  int sort_out_of_place;              // Sort method
  int * ALIGNED(128) partition;       // Static array indexed 0:
  /**/                                // (nx+2)*(ny+2)*(nz+2).  Each value
  /**/                                // corresponds to the associated particle
  /**/                                // array index of the first particle in
  /**/                                // the cell.  Array is allocated and
  /**/                                // values computed in sort_p.  Purpose is
  /**/                                // for implementing collision models
  /**/                                // This is given in terms of the
  /**/                                // underlying's grids space filling
  /**/                                // curve indexing.  Thus, immediately
  /**/                                // after a sort:
  /**/                                //   sp->p[sp->partition[g->sfc[i]  ]:
  /**/                                //         sp->partition[g->sfc[i]+1]-1]
  /**/                                // are all the particles in voxel
  /**/                                // with local index i, while:
  /**/                                //   sp->p[ sp->partition[ j   ]:
  /**/                                //          sp->partition[ j+1 ] ]
  /**/                                // are all the particles in voxel
  /**/                                // with space filling curve index j.
  /**/                                // Note: SFC NOT IN USE RIGHT NOW THUS
  /**/                                // g->sfc[i]=i ABOVE.

  grid_t * g;                         // Underlying grid
  species_id id;                      // Unique identifier for a species
  struct species *next;               // Next species in the list
} species_t;

BEGIN_C_DECLS

// In species_advance.c

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

species_t *
append_species( species_t * sp,
                species_t ** sp_list );

species_t *
species( const char * name,
         float q,
         float m,
         int max_local_np,
         int max_local_nm,
         int sort_interval,
         int sort_out_of_place,
         grid_t * g );

// FIXME: TEMPORARY HACK UNTIL THIS SPECIES_ADVANCE KERNELS
// CAN BE CONSTRUCTED ANALOGOUS TO THE FIELD_ADVANCE KERNELS
// (THESE FUNCTIONS ARE NECESSARY FOR HIGHER LEVEL CODE)

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

void
accumulate_rhob( field_t          * RESTRICT ALIGNED(128) f,
                 const particle_t * RESTRICT ALIGNED(32)  p,
                 const grid_t     * RESTRICT              g,
                 const float                              qsp );

// In hydro_p.c

void
accumulate_hydro_p( /**/  hydro_array_t        * RESTRICT ha,
                    const species_t            * RESTRICT sp,
                    const interpolator_array_t * RESTRICT ia );

// In move_p.cxx

int
move_p( particle_t       * ALIGNED(128) p0,    // Particle array
        particle_mover_t * ALIGNED(16)  m,     // Particle mover to apply
        accumulator_t    * ALIGNED(128) a0,    // Accumulator to use
        const grid_t     *              g,     // Grid parameters
        const float                     qsp ); // Species particle charge

END_C_DECLS

#endif // _species_advance_h_
