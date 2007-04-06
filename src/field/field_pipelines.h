#ifndef _field_pipelines_h_
#define _field_pipelines_h_

#include <field.h>
#include <pipeline.h>

/* THE FUNCTIONS IN THIS FILE ARE INTENDED TO RUN ON THE PROCESSORS
   THAT EXECUTE PIPELINES (E.G. SPE's ON THE CELL PROCESSOR). */

BEGIN_C_DECLS

#define DECLARE_PIPELINE(name)                                              \
void name##_pipeline(    name##_pipeline_args_t *args, int pipeline_rank ); \
void name##_pipeline_v4( name##_pipeline_args_t *args, int pipeline_rank ) 

/* In distribute_voxels.c */

/* Given a block of voxels to be processed by the pipelines, determine
   the number of voxels and the first voxel a particular job assigned
   to a pipeline should process.  The return voxel is the number of
   voxels to process.

   It is assumed that the pipelines will process voxels in FORTRAN
   ordering (e.g. inner loop increments x-index). */
   
int
distribute_voxels( int x0,  int x1,    /* range of x-indices (inclusive) */
                   int y0,  int y1,    /* range of y-indices (inclusive) */
                   int z0,  int z1,    /* range of z-indices (inclusive) */
                   int job, int n_job, /* job ... on [0,n_job-1] */
                   int * _x, int * _y, int * _z ); /* first voxel to process */

int
distribute_voxels_v4( int x0,  int x1,    /* range of x-indices (inclusive) */
                      int y0,  int y1,    /* range of y-indices (inclusive) */
                      int z0,  int z1,    /* range of z-indices (inclusive) */
                      int job, int n_job, /* job ... on [0,n_job-1] */
                      int * _x, int * _y, int * _z ); /* first voxel to process */

/* Time stepping pipelines ***************************************************/

/* In interpolator_pipeline.c */

/* This runs the interpolator kernel on a slice of the voxels in
   (1:nx,1:ny,1:nz).  This function is matched to the
   load_interpolator function in the parent directory. */

typedef struct load_interpolator_pipeline_args {
  interpolator_t * ALIGNED fi;
  const field_t  * ALIGNED f;
  const grid_t   *         g;
} load_interpolator_pipeline_args_t;

DECLARE_PIPELINE(load_interpolator);

/* In clear_accumulators_pipeline.c */

/* This has each pipeline clear the accumulator assigned to it.  This function
   is matched to the clear_accumulators function in the parent directory. */
   
typedef struct clear_accumulators_pipeline_args {
  accumulator_t * ALIGNED a; /* Base of all the accumulators */
  int n_voxel;
} clear_accumulators_pipeline_args_t;

DECLARE_PIPELINE(clear_accumulators);

/* In reduce_accumulators_pipeline.c */

/* This runs the interpolator kernel on a slice of the voxels in
   (1:nx,1:ny,1:nz).  This function is matched to the
   load_interpolator function in the parent directory. */

typedef struct reduce_accumulators_pipeline_args {
  accumulator_t * ALIGNED a;
  const grid_t  *         g;
} reduce_accumulators_pipeline_args_t;

DECLARE_PIPELINE(reduce_accumulators);

/* In unload_accumulator_pipeline.c */

/* This runs the accumulator kernel on a slice of the voxels in
   (1:nx,1:ny,1:nz).  This function is matched to the
   unload_accumulator function in the parent directory. */

typedef struct unload_accumulator_pipeline_args {
  field_t             * ALIGNED f;
  const accumulator_t * ALIGNED a;
  const grid_t        *         g;
} unload_accumulator_pipeline_args_t;

DECLARE_PIPELINE(unload_accumulator);

/* In advance_b_pipeline.c */

/* This runs the advance_b kernel on a slice of the voxels in
   (1:nx,1:ny,1:nz).  This function is matched to the advance_b
   function in the parent directory. */

typedef struct advance_b_pipeline_args {
  field_t      * ALIGNED f;
  const grid_t *         g;
  float frac;
} advance_b_pipeline_args_t;

DECLARE_PIPELINE(advance_b);

/* In advance_e_pipeline.c */

/* This runs the advance_e kernel on a slice of the voxels in
   (2:nx,2:ny,2:nz).  This function is matched to the advance_e
   function in the parent directory.  compute_curl_b is structurally
   identical to advance_e.  Changes in one likely need to be
   reflected in the other. */

typedef struct advance_e_pipeline_args {
  field_t                      * ALIGNED f;
  const material_coefficient_t * ALIGNED m;
  const grid_t                 *         g;
} advance_e_pipeline_args_t;

DECLARE_PIPELINE(advance_e);

/* Diagnostic pipelines ******************************************************/

/* In energy_f_pipeline.c */

/* This runs the energy_f kernel on a slice of the voxels in
   (1:nx,1:ny,1:nz).  This function is matched to the energy_f
   function in the parent directory.  Each pipelines gets its
   own args to allow them to return energies. */

typedef struct energy_f_pipeline_args {
  const field_t                * ALIGNED f;
  const material_coefficient_t * ALIGNED m;
  const grid_t                 *         g;
  double en[MAX_PIPELINE][6];
} energy_f_pipeline_args_t;

DECLARE_PIPELINE(energy_f);

/* In compute_curl_b_pipeline.c */

/* This runs the compute_curl_b kernel on a slice of the voxels in
   (2:nx,2:ny,2:nz).  This function is matched to the compute_curl_b
   function in the parent directory.  compute_curl_b is structurally
   identical to advance_e.  Changes in one likely need to be
   reflected in the other. */

typedef struct compute_curl_b_pipeline_args {
  field_t                      * ALIGNED f;
  const material_coefficient_t * ALIGNED m;
  const grid_t                 *         g;
} compute_curl_b_pipeline_args_t;

DECLARE_PIPELINE(compute_curl_b);

/* Divergence cleaning pipelines *********************************************/

/* In compute_rhob_pipeline.c */

/* This runs the compute_rhob kernel on a slice of the voxels in
   (2:nx,2:ny,2:nz).  This function is matched to the compute_rhob
   function in the parent directory.  Note that this pipeline is
   structurally the same as compute_div_e_err pipeline; updates to
   one likely should be reflected in the other. */

typedef struct compute_rhob_pipeline_args {
  field_t                      * ALIGNED f;
  const material_coefficient_t * ALIGNED m;
  const grid_t                 *         g;
} compute_rhob_pipeline_args_t;

DECLARE_PIPELINE(compute_rhob);

/* In compute_div_e_err_pipeline.c */

/* This runs the compute_div_e_err kernel on a slice of the voxels in
   (2:nx,2:ny,2:nz).  This function is matched to the compute_div_e_err
   function in the parent directory.  Note that this pipeline is
   structurally the same as compute_div_e_err pipeline; updates to
   one likely should be reflected in the other. */

typedef struct compute_div_e_err_pipeline_args {
  field_t                      * ALIGNED f;
  const material_coefficient_t * ALIGNED m;
  const grid_t                 *         g;
} compute_div_e_err_pipeline_args_t;

DECLARE_PIPELINE(compute_div_e_err);

/* In compute_rms_div_e_err_pipeline.c */

/* This runs the compute_rms_div_e_err kernel on a slice of the voxels in
   (2:nx,2:ny,2:nz).  This function is matched to the compute_rms_div_e_err
   function in the parent directory.   Each pipeline gets its own arguments
   so it can pass the contribution to err back to the host.  */

typedef struct compute_rms_div_e_err_pipeline_args {
  field_t      * ALIGNED f;
  const grid_t *         g;
  double err[MAX_PIPELINE];
} compute_rms_div_e_err_pipeline_args_t;

DECLARE_PIPELINE(compute_rms_div_e_err);

/* In clean_div_e_pipeline.c */

/* This runs the clean_e_kernel on a slice of the voxels in
   (1:nx,1:ny,1:nz).  This function is matched to the clean_div_e
   function in the parent directory.  */

typedef struct clean_div_e_pipeline_args {
  field_t                      * ALIGNED f;
  const material_coefficient_t * ALIGNED m;
  const grid_t                 *         g;
} clean_div_e_pipeline_args_t;

DECLARE_PIPELINE(clean_div_e);

/* In compute_div_b_err_pipeline.c */

/* This runs the compute_div_b_err kernel on a slice of the voxels in
   (1:nx,1:ny,1:nz).  This function is matched to the compute_div_b_err
   function in the parent directory.  */

typedef struct compute_div_b_err_pipeline_args {
  field_t      * ALIGNED f;
  const grid_t *         g;
} compute_div_b_err_pipeline_args_t;

DECLARE_PIPELINE(compute_div_b_err);

/* In compute_rms_div_b_err_pipeline.c */

/* This runs the compute_rms_div_b_err kernel on a slice of the voxels in
   (1:nx,1:ny,1:nz).  This function is matched to the compute_rms_div_b_err
   function in the parent directory.   Each pipeline gets its own arguments
   so it can pass the err back to the host.  */

typedef struct compute_rms_div_b_err_pipeline_args {
  field_t      * ALIGNED f;
  const grid_t *         g;
  double err[MAX_PIPELINE];
} compute_rms_div_b_err_pipeline_args_t;

DECLARE_PIPELINE(compute_rms_div_b_err);

/* In clean_div_b_pipeline.c */

/* This runs the clean_b_kernel on a slice of the voxels in
   (2:nx,2:ny,2:nz).  This function is matched to the clean_div_b
   function in the parent directory.  */

typedef struct clean_div_b_pipeline_args {
  field_t      * ALIGNED f;
  const grid_t *         g;
} clean_div_b_pipeline_args_t;

DECLARE_PIPELINE(clean_div_b);

#undef DECLARE_PIPELINE

END_C_DECLS

#endif /* _field_pipelines_h_ */
