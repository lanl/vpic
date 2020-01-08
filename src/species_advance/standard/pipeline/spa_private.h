#ifndef _spa_private_h_
#define _spa_private_h_

#ifndef IN_spa
#error "Do not include spa_private.h; include species_advance.h"
#endif

#include "../../species_advance.h"

///////////////////////////////////////////////////////////////////////////////
// advance_p_pipeline interface

typedef struct particle_mover_seg
{
    MEM_PTR( particle_mover_t, 16 ) pm; // First mover in segment
    int max_nm;                         // Maximum number of movers
    int nm;                             // Number of movers used
    int n_ignored;                      // Number of movers ignored

    PAD_STRUCT( SIZEOF_MEM_PTR + 3 * sizeof( int ) )

} particle_mover_seg_t;

typedef struct advance_p_pipeline_args
{
    MEM_PTR( particle_t, 128 ) p0;            // Particle array
    MEM_PTR( particle_mover_t, 128 ) pm;      // Particle mover array
    MEM_PTR( accumulator_t, 128 ) a0;         // Accumulator arrays
    MEM_PTR( const interpolator_t, 128 ) f0;  // Interpolator array
    MEM_PTR( particle_mover_seg_t, 128 ) seg; // Dest for return values
    MEM_PTR( const grid_t, 1 ) g;             // Local domain grid params

    float qdt_2mc; // Particle/field coupling
    float cdt_dx;  // x-space/time coupling
    float cdt_dy;  // y-space/time coupling
    float cdt_dz;  // z-space/time coupling
    float qsp;     // Species particle charge

    int np;     // Number of particles
    int max_nm; // Number of movers
    int nx;     // x-mesh resolution
    int ny;     // y-mesh resolution
    int nz;     // z-mesh resolution

    PAD_STRUCT( 6 * SIZEOF_MEM_PTR + 5 * sizeof( float ) + 5 * sizeof( int ) )

} advance_p_pipeline_args_t;

// PROTOTYPE_PIPELINE( advance_p, advance_p_pipeline_args_t );

void advance_p_pipeline_scalar( advance_p_pipeline_args_t *args,
                                int pipeline_rank, int n_pipeline );

void advance_p_pipeline_v4( advance_p_pipeline_args_t *args, int pipeline_rank,
                            int n_pipeline );

void advance_p_pipeline_v8( advance_p_pipeline_args_t *args, int pipeline_rank,
                            int n_pipeline );

void advance_p_pipeline_v16( advance_p_pipeline_args_t *args, int pipeline_rank,
                             int n_pipeline );

///////////////////////////////////////////////////////////////////////////////
// center_p_pipeline and uncenter_p_pipeline interface

typedef struct center_p_pipeline_args
{
    MEM_PTR( particle_t, 128 ) p0;           // Particle array
    MEM_PTR( const interpolator_t, 128 ) f0; // Interpolator array
    float qdt_2mc;                           // Particle/field coupling
    int np;                                  // Number of particles

    PAD_STRUCT( 2 * SIZEOF_MEM_PTR + sizeof( float ) + sizeof( int ) )

} center_p_pipeline_args_t;

// PROTOTYPE_PIPELINE( center_p,   center_p_pipeline_args_t );

void center_p_pipeline_scalar( center_p_pipeline_args_t *args,
                               int pipeline_rank, int n_pipeline );

void center_p_pipeline_v4( center_p_pipeline_args_t *args, int pipeline_rank,
                           int n_pipeline );

void center_p_pipeline_v8( center_p_pipeline_args_t *args, int pipeline_rank,
                           int n_pipeline );
void center_p_pipeline_v16( center_p_pipeline_args_t *args, int pipeline_rank,
                            int n_pipeline );

// PROTOTYPE_PIPELINE( uncenter_p, center_p_pipeline_args_t );

void uncenter_p_pipeline_scalar( center_p_pipeline_args_t *args,
                                 int pipeline_rank, int n_pipeline );

void uncenter_p_pipeline_v4( center_p_pipeline_args_t *args, int pipeline_rank,
                             int n_pipeline );

void uncenter_p_pipeline_v8( center_p_pipeline_args_t *args, int pipeline_rank,
                             int n_pipeline );

void uncenter_p_pipeline_v16( center_p_pipeline_args_t *args, int pipeline_rank,
                              int n_pipeline );

///////////////////////////////////////////////////////////////////////////////
// energy_p_pipeline interface

typedef struct energy_p_pipeline_args
{
    MEM_PTR( const particle_t, 128 ) p;     // Particle array
    MEM_PTR( const interpolator_t, 128 ) f; // Interpolator array
    MEM_PTR( double, 128 ) en;              // Return values
    float qdt_2mc;                          // Particle/field coupling
    float msp;                              // Species particle rest mass
    int np;                                 // Number of particles

    PAD_STRUCT( 3 * SIZEOF_MEM_PTR + 2 * sizeof( float ) + sizeof( int ) )

} energy_p_pipeline_args_t;

// PROTOTYPE_PIPELINE( energy_p, energy_p_pipeline_args_t );

void energy_p_pipeline_scalar( energy_p_pipeline_args_t *RESTRICT args,
                               int pipeline_rank, int n_pipeline );

void energy_p_pipeline_v4( energy_p_pipeline_args_t *args, int pipeline_rank,
                           int n_pipeline );

void energy_p_pipeline_v8( energy_p_pipeline_args_t *args, int pipeline_rank,
                           int n_pipeline );

void energy_p_pipeline_v16( energy_p_pipeline_args_t *args, int pipeline_rank,
                            int n_pipeline );

///////////////////////////////////////////////////////////////////////////////
// sort_p_pipeline interface

// Given the voxel index, compute which subsort is responsible for
// sorting particles within that voxel.  This takes into account
// that v*P might overflow 32-bits and that only voxels [vl,vh]
// may contain particles.  This macro is mostly robust.

#define V2P( v, P, vl, vh )                                                    \
    ( ( ( ( int64_t )( ( v ) - ( vl ) ) ) * ( ( int64_t )( P ) ) ) /           \
      ( ( int64_t )( ( vh ) - ( vl ) + 1 ) ) )

// Given the pipeline rank, compute the first voxel a subsort is
// responsible for handling.  This is based on:
//   p = floor(vP/V) =>
//   p <= vP/V < p+1 =>
//   pV/P <= v < (p+1)V/P
// The range of voxels which satisfy this inequality is then:
//   [ ceil(pV/P), ceil((p+1)V/P) )
// In integer math, the lower bound is thus:
//   v = (p*V + P-1)/P
// where v above is v-vl and V = vh-vl+1.  This takes into account
// that p*V might overflow 32-bits.  This macro is mostly robust.

#define P2V( p, P, vl, vh )                                                    \
    ( ( vl ) +                                                                 \
      ( ( ( ( int64_t )( p ) ) * ( ( int64_t )( ( vh ) - ( vl ) + 1 ) ) +      \
          ( ( int64_t )( (P)-1 ) ) ) /                                         \
        ( ( int64_t )( P ) ) ) )

// FIXME: safe to remove? enum { max_subsort_voxel = 26624 };

typedef struct sort_p_pipeline_args
{
    MEM_PTR( particle_t, 128 ) p;         // Particles (0:n-1)
    MEM_PTR( particle_t, 128 ) aux_p;     // Aux particle atorage (0:n-1)
    MEM_PTR( int, 128 ) coarse_partition; // Coarse partition storage
    /**/                                    // (0:max_subsort-1,0:MAX_PIPELINE-1)
    MEM_PTR( int, 128 ) partition;        // Partitioning (0:n_voxel)
    MEM_PTR( int, 128 ) next;             // Aux partitioning (0:n_voxel)
    int n;                                // Number of particles
    int n_subsort; // Number of pipelines to be used for subsorts
    int vl, vh;    // Particles may be contained in voxels [vl,vh].
    int n_voxel;   // Number of voxels total (including ghosts)

    PAD_STRUCT( 5 * SIZEOF_MEM_PTR + 5 * sizeof( int ) )

} sort_p_pipeline_args_t;

// PROTOTYPE_PIPELINE( coarse_count, sort_p_pipeline_args_t );
// PROTOTYPE_PIPELINE( coarse_sort,  sort_p_pipeline_args_t );
// PROTOTYPE_PIPELINE( subsort,      sort_p_pipeline_args_t );

void coarse_count_pipeline_scalar( sort_p_pipeline_args_t *args,
                                   int pipeline_rank, int n_pipeline );

void coarse_sort_pipeline_scalar( sort_p_pipeline_args_t *args,
                                  int pipeline_rank, int n_pipeline );

void subsort_pipeline_scalar( sort_p_pipeline_args_t *args, int pipeline_rank,
                              int n_pipeline );

#endif // _spa_private_h_
