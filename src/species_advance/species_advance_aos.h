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

#ifndef _species_advance_aos_h_
#define _species_advance_aos_h_

typedef int32_t species_id; // Must be 32-bit wide for particle_injector_t

// FIXME: Eventually particle_t (definitely) and their other formats
// (maybe) should be opaque and specific to a particular
// species_advance implementation

typedef struct particle
{
    float dx, dy, dz; // Particle position in cell coordinates (on [-1,1])
    int32_t i;        // Voxel containing the particle.  Note that
    /**/                // particles awaiting processing by boundary_p
    /**/                // have actually set this to 8*voxel + face where
    /**/                // face is the index of the face they interacted
    /**/                // with (on 0:5).  This limits the local number of
    /**/                // voxels to 2^28 but emitter handling already
    /**/                // has a stricter limit on this (2^26).
    float ux, uy, uz; // Particle normalized momentum
    float w;          // Particle weight (number of physical particles)
} particle_t;

// WARNING: FUNCTIONS THAT USE A PARTICLE_MOVER ASSUME THAT EVERYBODY
// WHO USES THAT PARTICLE MOVER WILL HAVE ACCESS TO PARTICLE ARRAY

typedef struct particle_mover
{
    float dispx, dispy, dispz; // Displacement of particle
    int32_t i;                 // Index of the particle to move
} particle_mover_t;

// NOTE: THE LAYOUT OF A PARTICLE_INJECTOR _MUST_ BE COMPATIBLE WITH
// THE CONCATENATION OF A PARTICLE_T AND A PARTICLE_MOVER!

typedef struct particle_injector
{
    float dx, dy, dz;          // Particle position in cell coords (on [-1,1])
    int32_t i;                 // Index of cell containing the particle
    float ux, uy, uz;          // Particle normalized momentum
    float w;                   // Particle weight (number of physical particles)
    float dispx, dispy, dispz; // Displacement of particle
    species_id sp_id;          // Species of particle
} particle_injector_t;

typedef struct species
{
    char *name; // Species name
    float q;    // Species particle charge
    float m;    // Species particle rest mass

    int np, max_np;               // Number and max local particles
    particle_t *ALIGNED( 128 ) p; // Array of particles for the species

    int nm, max_nm;                      // Number and max local movers in use
    particle_mover_t *ALIGNED( 128 ) pm; // Particle movers

    int64_t last_sorted;           // Step when the particles were last
                                   // sorted.
    int sort_interval;             // How often to sort the species
    int sort_out_of_place;         // Sort method
    int *ALIGNED( 128 ) partition; // Static array indexed 0:
    /**/                             // (nx+2)*(ny+2)*(nz+2).  Each value
    /**/                             // corresponds to the associated particle
    /**/                             // array index of the first particle in
    /**/                             // the cell.  Array is allocated and
    /**/                             // values computed in sort_p.  Purpose is
    /**/                             // for implementing collision models
    /**/                             // This is given in terms of the
    /**/                             // underlying's grids space filling
    /**/                             // curve indexing.  Thus, immediately
    /**/                             // after a sort:
    /**/                             //   sp->p[sp->partition[g->sfc[i]  ]:
    /**/                             //         sp->partition[g->sfc[i]+1]-1]
    /**/                             // are all the particles in voxel
    /**/                             // with local index i, while:
    /**/                             //   sp->p[ sp->partition[ j   ]:
    /**/                             //          sp->partition[ j+1 ] ]
    /**/                             // are all the particles in voxel
    /**/                             // with space filling curve index j.
    /**/                             // Note: SFC NOT IN USE RIGHT NOW THUS
    /**/                             // g->sfc[i]=i ABOVE.

    grid_t *g;            // Underlying grid
    species_id id;        // Unique identifier for a species
    struct species *next; // Next species in the list
} species_t;

#endif // _species_advance_aos_h_
