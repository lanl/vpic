#ifndef _spa_private_h_
#define _spa_private_h_

#ifndef IN_spa
#error "Do not include spa_private.h; include species_advance.h"
#endif

#include "spa.h"

#define FOR_SPU ( defined(CELL_SPU_BUILD)        || \
                  ( defined(CELL_PPU_BUILD)    &&   \
                     defined(USE_CELL_SPUS)    &&   \
                     defined(HAS_SPU_PIPELINE) ) )

#if FOR_SPU

# if defined(CELL_PPU_BUILD) 

    // Use SPU dispatcher on the SPU pipeline
    // PPU will do straggler cleanup with scalar pipeline

#   define EXEC_PIPELINES(name,args,sz_args)                      \
    spu.dispatch( (pipeline_func_t)( (size_t)                     \
        ( root_segment_##name##_pipeline_spu) ), args, sz_args ); \
    name##_pipeline( args, spu.n_pipeline, spu.n_pipeline )

#   define WAIT_PIPELINES() spu.wait()

#   define N_PIPELINE spu.n_pipeline

#   define PROTOTYPE_PIPELINE( name, args_t )                      \
    extern uint32_t root_segment_##name##_pipeline_spu;            \
                                                                   \
    void                                                           \
    name##_pipeline( args_t * args,                                \
                     int pipeline_rank,                            \
                     int n_pipeline )

#   define PAD_STRUCT( sz ) char _pad[ PAD( (sz), 16 ) ];

# else

    // SPUs cannot dispatch pipelines

#   define PROTOTYPE_PIPELINE( name, args_t )                   \
    void                                                        \
    _SPUEAR_##name##_pipeline_spu( MEM_PTR( args_t, 128 ) argp, \
                                   int pipeline_rank,           \
                                   int n_pipeline )

#   define PAD_STRUCT( sz ) char _pad[ PAD( (sz), 16 ) ];

# endif

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  // Use thread dispatcher on the v4 pipeline
  // Caller will do straggler cleanup with scalar pipeline

# define EXEC_PIPELINES(name,args,sz_args)                               \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v4, args, sz_args ); \
  name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

# define WAIT_PIPELINES() thread.wait()

# define N_PIPELINE thread.n_pipeline

# define PROTOTYPE_PIPELINE( name, args_t ) \
  void                                      \
  name##_pipeline_v4( args_t * args,        \
                      int pipeline_rank,    \
                      int n_pipeline );     \
                                            \
  void                                      \
  name##_pipeline( args_t * args,           \
                   int pipeline_rank,       \
                   int n_pipeline )

# define PAD_STRUCT( sz )

#else

  // Use thread dispatcher on the scalar pipeline
  // Caller will do straggler cleanup with scalar pipeline

# define EXEC_PIPELINES(name,args,sz_args)                              \
  thread.dispatch( (pipeline_func_t)name##_pipeline, args, sz_args );   \
  name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

# define WAIT_PIPELINES() thread.wait()

# define N_PIPELINE thread.n_pipeline

# define PROTOTYPE_PIPELINE( name, args_t ) \
  void                                      \
  name##_pipeline( args_t * args,           \
                   int pipeline_rank,       \
                   int n_pipeline )

# define PAD_STRUCT( sz )

#endif

///////////////////////////////////////////////////////////////////////////////
// advance_p_pipeline interface

typedef struct particle_mover_seg {

  MEM_PTR( particle_mover_t, 16 ) pm; // First mover in segment
  int max_nm;                         // Maximum number of movers
  int nm;                             // Number of movers used
  int n_ignored;                      // Number of movers ignored

  PAD_STRUCT( SIZEOF_MEM_PTR+3*sizeof(int) )

} particle_mover_seg_t;

typedef struct advance_p_pipeline_args {

  MEM_PTR( particle_t,           128 ) p0;       // Particle array
  MEM_PTR( particle_mover_t,     128 ) pm;       // Particle mover array
  MEM_PTR( accumulator_t,        128 ) a0;       // Accumulator arrays
  MEM_PTR( const interpolator_t, 128 ) f0;       // Interpolator array
  MEM_PTR( particle_mover_seg_t, 128 ) seg;      // Dest for return values
  MEM_PTR( const grid_t,         1   ) g;        // Local domain grid params

# if FOR_SPU

  // For move_p_spu; it is easier to have the PPU unpack these grid_t
  // quantities for the SPUs than to have the SPUs pointer chase
  // through the above grid_t to extract these quantities.

  MEM_PTR( const int64_t,        128 ) neighbor; // Global voxel indices of
  /**/                                           // voxels adjacent to local
  /**/                                           // voxels
  int64_t                              rangel;   // First global voxel here
  int64_t                              rangeh;   // Last global voxel here

# endif

  float                                qdt_2mc;  // Particle/field coupling
  float                                cdt_dx;   // x-space/time coupling
  float                                cdt_dy;   // y-space/time coupling
  float                                cdt_dz;   // z-space/time coupling

  int                                  np;       // Number of particles
  int                                  max_nm;   // Number of movers
  int                                  nx;       // x-mesh resolution
  int                                  ny;       // y-mesh resolution
  int                                  nz;       // z-mesh resolution
 
# if FOR_SPU
  PAD_STRUCT( 7*SIZEOF_MEM_PTR + 4*sizeof(float) + 5*sizeof(int) + 2*sizeof(int64_t))
# else
  PAD_STRUCT( 6*SIZEOF_MEM_PTR + 4*sizeof(float) + 5*sizeof(int) )
# endif

} advance_p_pipeline_args_t;

PROTOTYPE_PIPELINE( advance_p, advance_p_pipeline_args_t );

///////////////////////////////////////////////////////////////////////////////
// center_p_pipeline and uncenter_p_pipeline interface

typedef struct center_p_pipeline_args {

  MEM_PTR( particle_t,           128 ) p0;      // Particle array
  MEM_PTR( const interpolator_t, 128 ) f0;      // Interpolator array
  float                                qdt_2mc; // Particle/field coupling
  int                                  np;      // Number of particles

  PAD_STRUCT( 2*SIZEOF_MEM_PTR + sizeof(float) + sizeof(int) )

} center_p_pipeline_args_t;

PROTOTYPE_PIPELINE( center_p,   center_p_pipeline_args_t );
PROTOTYPE_PIPELINE( uncenter_p, center_p_pipeline_args_t );

///////////////////////////////////////////////////////////////////////////////
// energy_p_pipeline interface

typedef struct energy_p_pipeline_args {

  MEM_PTR( const particle_t,     128 ) p0;      // Particle array
  MEM_PTR( const interpolator_t, 128 ) f0;      // Interpolator array
  MEM_PTR( double,               128 ) en;      // Return values
  float                                qdt_2mc; // Particle/field coupling
  int                                  np;      // Number of particles

  PAD_STRUCT( 3*SIZEOF_MEM_PTR + sizeof(float) + sizeof(int) )

} energy_p_pipeline_args_t;

PROTOTYPE_PIPELINE( energy_p, energy_p_pipeline_args_t );

///////////////////////////////////////////////////////////////////////////////
// sort_p_pipeline interface

enum {
  coarse_block_size = 8, // Larger than 8 saturates MFC queues in SPU
                         // course_sort
  max_coarse_bucket = 1+4*(MAX_PIPELINE/4) // Must be a multiple of 4
};

typedef struct sort_p_pipeline_args {

  MEM_PTR( particle_t, 128 ) p;                // Particles (0:n-1)
  MEM_PTR( particle_t, 128 ) aux_p;            // Aux particle atorage (0:n-1)
  MEM_PTR( int,        128 ) coarse_partition; // Coarse partition storage
  /**/ // (0:max_coarse_bucket-1,0:MAX_PIPELINE)
  MEM_PTR( int,        128 ) partition;        // Partitioning (0:n_voxel)
  MEM_PTR( int,        128 ) next;             // Aux partitioning (0:n_voxel)
  int n;               // Number of particles
  int n_voxel;         // Number of local voxels (including ghost voxels)
  int n_coarse_bucket; // Number of coarse sort buckets

  PAD_STRUCT( 5*SIZEOF_MEM_PTR + 3*sizeof(int) )

} sort_p_pipeline_args_t;

// FIXME: Temporary hack while SPU accelerating the sort
void coarse_sort_p( sort_p_pipeline_args_t * ALIGNED(128) args );

PROTOTYPE_PIPELINE( coarse_count, sort_p_pipeline_args_t );
PROTOTYPE_PIPELINE( coarse_sort,  sort_p_pipeline_args_t );
PROTOTYPE_PIPELINE( subsort,      sort_p_pipeline_args_t );

#undef FOR_SPU
#undef PAD_STRUCT
#undef PROTOTYPE_PIPELINE

#endif // _spa_private_h_
