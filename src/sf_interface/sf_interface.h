#ifndef _sf_interface_h_
#define _sf_interface_h_

// FIXME: THE HOST PROCESSED FIELD KERNELS SHOULD BE UPDATED TO USE
// SCALAR FMA INSTRUCTIONS WITH COMMENSURATE ROUND-OFF PROPERTIES TO
// THE FMA INSTRUCTIONS USED ON THE PIPELINE PROCESSED FIELDS!

// FIXME: (nx>1) ? (1/dx) : 0 TYPE LOGIC SHOULD BE FIXED SO THAT NX
// REFERS TO THE GLOBAL NUMBER OF CELLS IN THE X-DIRECTION (NOT THE
// _LOCAL_ NUMBER OF CELLS).  THIS LATENT BUG IS NOT EXPECTED TO
// AFFECT ANY PRACTICAL SIMULATIONS.

// FIXME: EXTERNAL DIAGNOSTICS THAT READ THESE WILL NEED TO BE UPDATED
// TO REFLECT SPU USAGE ALIGNMENT CHANGES.
  
#include "../field_advance/field_advance.h"
// FIXME: SHOULD INCLUDE SPECIES_ADVANCE TOO ONCE READY

// Hydro arrays shall be a (nx+2) x (ny+2) x (nz+2) allocation indexed
// FORTRAN style from (0:nx+1,0:ny+1,0:nz+1).  Hydros for voxels on
// the surface of the local domain (for example h(0,:,:) or
// h(nx+1,:,:)) are not used.

typedef struct hydro {
  float jx, jy, jz, rho; // Current and charge density => <q v_i f>, <q f>
  float px, py, pz, ke;  // Momentum and K.E. density  => <p_i f>, <m c^2 (gamma-1) f>
  float txx, tyy, tzz;   // Stress diagonal            => <p_i v_j f>, i==j
  float tyz, tzx, txy;   // Stress off-diagonal        => <p_i v_j f>, i!=j
# if ( defined(CELL_PPU_BUILD) || defined(CELL_SPU_BUILD) ) && defined(USE_CELL_SPUS)
  float _pad[2];         // 64-byte align (next power of two)
# else
  float _pad[2];         // 16-byte align
# endif
} hydro_t;

// Interpolator arrays shall be a (nx+2) x (ny+2) x (nz+2) allocation
// indexed FORTRAN style from (0:nx+1,0:ny+1,0:nz+1). Interpolators
// for voxels on the surface of the local domain (for example
// fi(0,:,:) or fi(nx+1,:,:)) are not used.

typedef struct interpolator {
  float ex, dexdy, dexdz, d2exdydz;
  float ey, deydz, deydx, d2eydzdx;
  float ez, dezdx, dezdy, d2ezdxdy;
  float cbx, dcbxdx;
  float cby, dcbydy;
  float cbz, dcbzdz;
# if ( defined(CELL_PPU_BUILD) || defined(CELL_SPU_BUILD) ) && defined(USE_CELL_SPUS)
  float _pad[2];
  int64_t neighbor[6]; // 128-byte align.  Use padding for neighbor data
# else
  float _pad[2];  // 16-byte align
# endif
} interpolator_t;

// Accumulator arrays shall be a
//   POW2_CEIL((nx+2)x(ny+2)x(nz+2),2)x(1+n_pipeline)
// allocation indexed FORTRAN style.  That is, the accumulator array
// is a 4d array.  a(:,:,:,0) is the accumulator used by the host
// processor.  a(:,:,:,1:n_pipeline) are the accumulators used by
// pipelines during operations.  Like the interpolator, accumualtors
// on the surface of the local domain are not used.

typedef struct accumulator {
  float jx[4];   // jx0@(0,-1,-1),jx1@(0,1,-1),jx2@(0,-1,1),jx3@(0,1,1)
  float jy[4];   // jy0@(-1,0,-1),jy1@(-1,0,1),jy2@(1,0,-1),jy3@(1,0,1)
  float jz[4];   // jz0@(-1,-1,0),jz1@(1,-1,0),jz2@(-1,1,0),jz3@(1,1,0)
# if ( defined(CELL_PPU_BUILD) || defined(CELL_SPU_BUILD) ) && defined(USE_CELL_SPUS)
  float _pad[4]; // 64-byte align (next power of two )
# else
  /**/           // 16-byte align
# endif
} accumulator_t;

BEGIN_C_DECLS

// In sf_structors.c

hydro_t * ALIGNED(128)
new_hydro( grid_t * g );

void
delete_hydro( hydro_t * ALIGNED(128) h );

void
clear_hydro( hydro_t       * ALIGNED(128) h,
             const grid_t  *              g );

interpolator_t * ALIGNED(128)
new_interpolator( grid_t * g );

void
delete_interpolator( interpolator_t * ALIGNED(128) fi );

accumulator_t * ALIGNED(128)
new_accumulators( grid_t * g );

void
delete_accumulators( accumulator_t * ALIGNED(128) a );

// In load_interpolator.c

// Going into load_interpolator, the field array f contains the
// current information such that the fields can be interpolated to
// particles within the local domain.  Load interpolate computes the
// field array into a set of interpolation coefficients for each voxel
// inside the local domain suitable for use by the particle update
// functions.

void
load_interpolator( interpolator_t * ALIGNED(128) fi,
                   const field_t  * ALIGNED(128) f,
                   const grid_t   *              g );

// In clear_accumulators.c

// This zeros out all the accumulator arrays in a pipelined fashion.

void
clear_accumulators( accumulator_t * ALIGNED(128) a,
                    const grid_t  *              g );

// In reduce_accumulators.c

// Going into reduce_accumulators, the host processor and the pipeline
// processors have each accumulated values to their personal
// accumulators.  This reduces the pipeline accumulators into the host
// accumulator with a pipelined horizontal reduction (a deterministic
// reduction).

void
reduce_accumulators( accumulator_t * ALIGNED(128) a,
                     const grid_t  *              g );

// In unload_accumulator.c

// Going into unload_accumulator, the accumulator contains 4 times the
// net amount of charge that crossed the quarter face associated with
// each accumulator component (has units of physical charge, i.e. C)
// computed by the advance_p functions.  unload_accumulator computes
// the physical current density (A/m^2 in MKS units) associated with
// all local quarter faces and accumulates the local quarter faces to
// local field array jf.  unload_accumulator assumes all the pipeline
// accumulators have been reduced into the host accumulator.

void
unload_accumulator( field_t             * ALIGNED(128) f, 
                    const accumulator_t * ALIGNED(128) a,
                    const grid_t        *              g );

// In hydro.c

void
synchronize_hydro( hydro_t      * ALIGNED(128) hydro,
                   const grid_t *              g );
             
void
local_adjust_hydro( hydro_t      * ALIGNED(128) h,
                    const grid_t *              g );

END_C_DECLS

#endif // _sf_interface_h_
