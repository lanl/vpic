#ifndef _particle_pipelines_h_
#define _particle_pipelines_h_

#ifndef IN_particle_pipeline
#error "Do not include particle_pipelines.h; include particle.h"
#endif

#include <particle.h>

///////////////////////////////////////////////////////////////////////////////
// advance_p_pipeline interface

typedef struct particle_mover_seg {

  MEM_PTR( particle_mover_t, 16 ) pm; // First mover in segment
  int max_nm;                         // Maximum number of movers
  int nm;                             // Number of movers used
  int n_ignored;                      // Number of movers ignored

# ifdef USE_CELL_SPUS // Pad to 16-byte boundary
  char _pad[ PAD( SIZEOF_MEM_PTR+3*sizeof(int), 16 ) ];
# endif

} particle_mover_seg_t;

typedef struct advance_p_pipeline_args {

  MEM_PTR( particle_t,           128 ) p0;       // Particle array
  MEM_PTR( particle_mover_t,     128 ) pm;       // Particle mover array
  MEM_PTR( accumulator_t,        128 ) a0;       // Accumulator arrays
  MEM_PTR( const interpolator_t, 128 ) f0;       // Interpolator array
  MEM_PTR( particle_mover_seg_t, 128 ) seg;      // Dest for return values
# ifdef USE_CELL_SPUS // For move_p_spu
  MEM_PTR( const int64_t,        128 ) neighbor; // Global voxel indices of
  /**/                                           // voxels adjacent to local
  /**/                                           // voxels
  int                                  rangel;   // First global voxel here
  int                                  rangeh;   // Last global voxel here
# else                // For move_p
  MEM_PTR( const grid_t,         1   ) g;        // Local domain grid params
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
 
# ifdef USE_CELL_SPUS // Align to 16-bytes
  char _pad[ PAD( 6*SIZEOF_MEM_PTR + 4*sizeof(float) + 7*sizeof(int), 16 ) ];
# endif

} advance_p_pipeline_args_t;

///////////////////////////////////////////////////////////////////////////////
// center_p_pipeline and uncenter_p_pipeline interface

typedef struct center_p_pipeline_args {

  MEM_PTR( particle_t,           128 ) p0;      // Particle array
  MEM_PTR( const interpolator_t, 128 ) f0;      // Interpolator array
  float                                qdt_2mc; // Particle/field coupling
  int                                  np;      // Number of particles

# ifdef USE_CELL_SPUS // Align to 16-bytes
  char _pad[ PAD( 2*SIZEOF_MEM_PTR + sizeof(float) + sizeof(int), 16 ) ];
# endif

} center_p_pipeline_args_t;

///////////////////////////////////////////////////////////////////////////////
// energy_p_pipeline interface

typedef struct energy_p_pipeline_args {

  MEM_PTR( const particle_t,     128 ) p0;      // Particle array
  MEM_PTR( const interpolator_t, 128 ) f0;      // Interpolator array
  MEM_PTR( double,               128 ) en;      // Return values
  float                                qdt_2mc; // Particle/field coupling
  int                                  np;      // Number of particles

# ifdef USE_CELL_SPUS // Align to 16-bytes
  char _pad[ PAD( 3*SIZEOF_MEM_PTR + sizeof(float) + sizeof(int), 16 ) ];
# endif

} energy_p_pipeline_args_t;

#endif // _particle_pipelines_h_
