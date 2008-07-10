#ifndef _spa_private_h_
#define _spa_private_h_

#ifndef IN_spa
#error "Do not include spa_private.h; include species_advance.h"
#endif

#include "spa.h"
//#include "pipeline_control.h"

#define FOR_SPU ( defined(CELL_SPU_BUILD)        || \
                  ( defined(CELL_PPU_BUILD)    &&   \
                     defined(USE_CELL_SPUS)    &&   \
                     defined(HAS_SPU_PIPELINE) ) )

#if FOR_SPU

//    #include <spe_events.h>
    #define READ_ARGS_AND_ADVANCE 0
    #define ADVANCE_COMPLETE 1
    #define END_EVENT_LOOP 3

# if defined(CELL_PPU_BUILD) 

    // Use SPU dispatcher on the SPU pipeline

#   define INIT_PIPELINES(name,args,sz_args) \
    spu.dispatch( SPU_PIPELINE(name##_pipeline_spu), args, sz_args );

#   define EXEC_PIPELINES(name,args,sz_args) \
	spu.signal(READ_ARGS_AND_ADVANCE); \
    name##_pipeline( args, spu.n_pipeline, spu.n_pipeline )

#   define WAIT_PIPELINES() spu.sync(ADVANCE_COMPLETE)

#   define FINALIZE_PIPELINES() \
    spu.signal(END_EVENT_LOOP); \
    spu.wait()

#   define N_PIPELINE       spu.n_pipeline

# else

    // SPUs cannot dispatch pipelines

# endif

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  // Use thread dispatcher on the v4 pipeline

# define INIT_PIPELINES(name,args,sz_args)

# define EXEC_PIPELINES(name,args,sz_args)                               \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v4, args, sz_args ); \
  name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

# define WAIT_PIPELINES() thread.wait()

# define FINALIZE_PIPELINES()

# define N_PIPELINE       thread.n_pipeline

#else

  // Use thread dispatcher on the scalar pipeline

# define INIT_PIPELINES(name,args,sz_args)

# define EXEC_PIPELINES(name,args,sz_args)                              \
  thread.dispatch( (pipeline_func_t)name##_pipeline, args, sz_args );   \
  name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

# define WAIT_PIPELINES() thread.wait()

# define FINALIZE_PIPELINES()

# define N_PIPELINE       thread.n_pipeline

#endif

///////////////////////////////////////////////////////////////////////////////
// advance_p_pipeline interface

typedef struct particle_mover_seg {

  MEM_PTR( particle_mover_t, 16 ) pm; // First mover in segment
  int max_nm;                         // Maximum number of movers
  int nm;                             // Number of movers used
  int n_ignored;                      // Number of movers ignored

# if FOR_SPU // Align to 16-bytes
  char _pad[ PAD( SIZEOF_MEM_PTR+3*sizeof(int), 16 ) ];
# endif

} particle_mover_seg_t;

typedef struct advance_p_pipeline_args {

  MEM_PTR( particle_t,           128 ) p0;       // Particle array
  MEM_PTR( particle_mover_t,     128 ) pm;       // Particle mover array
  MEM_PTR( accumulator_t,        128 ) a0;       // Accumulator arrays
  MEM_PTR( const interpolator_t, 128 ) f0;       // Interpolator array
  MEM_PTR( particle_mover_seg_t, 128 ) seg;      // Dest for return values
  MEM_PTR( const grid_t,         1   ) g;        // Local domain grid params

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

  // For move_p_spu; it is easier to have the PPU unpack these grid_t
  // quantities for the SPUs than to have the SPUs pointer chase
  // through the above grid_t to extract these quantities.

  MEM_PTR( const int64_t,        128 ) neighbor; // Global voxel indices of
  /**/                                           // voxels adjacent to local
  /**/                                           // voxels
  int                                  rangel;   // First global voxel here
  int                                  rangeh;   // Last global voxel here

  // Align to 16-bytes

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

# if FOR_SPU // Align to 16-bytes
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

# if FOR_SPU // Align to 16-bytes
  char _pad[ PAD( 3*SIZEOF_MEM_PTR + sizeof(float) + sizeof(int), 16 ) ];
# endif

} energy_p_pipeline_args_t;

#undef FOR_SPU

#endif // _spa_private_h_
