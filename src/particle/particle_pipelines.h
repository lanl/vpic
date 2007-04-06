#ifndef _particle_pipelines_h_
#define _particle_pipelines_h_

#include <particle.h>
#include <pipeline.h>

BEGIN_C_DECLS

#define DECLARE_PIPELINE(name)                                              \
void name##_pipeline(    name##_pipeline_args_t *args, int pipeline_rank ); \
void name##_pipeline_v4( name##_pipeline_args_t *args, int pipeline_rank ) 

/* In advance_p_pipeline.c */

typedef struct advance_p_pipeline_args {
  particle_t           * ALIGNED p;   /* Particle array */
  int                            n;   /* Number of particles */
  float                          q_m; /* Charge to mass ratio */
  particle_mover_t     * ALIGNED pm;  /* Particle mover array */
  int                            nm;  /* Number of movers */
  accumulator_t        * ALIGNED a;   /* Accumuator arrays */
  const interpolator_t * ALIGNED f;   /* Interpolator array */
  const grid_t         *         g;   /* Local domain grid parameters */
 
  /* Return values */

  struct {
    particle_mover_t * ALIGNED pm; /* First mover in segment */
    int nm;                        /* Number of used movers in segment */
  } seg[MAX_PIPELINE+1];           /* seg[n_pipeline] used by host */

} advance_p_pipeline_args_t;

DECLARE_PIPELINE(advance_p);

/* In center_p_pipeline.c */

typedef struct center_p_pipeline_args {
  particle_t           * ALIGNED p;   /* Particle array */
  int                            n;   /* Number of particles */
  float                          q_m; /* Charge to mass ratio */
  const interpolator_t * ALIGNED f;   /* Interpolator array */
  const grid_t         *         g;   /* Local domain grid parameters */
} center_p_pipeline_args_t;

DECLARE_PIPELINE(center_p);

/* In uncenter_p_pipeline.c */

typedef struct uncenter_p_pipeline_args {
  particle_t           * ALIGNED p;   /* Particle array */
  int                            n;   /* Number of particle */
  float                          q_m; /* Charge to mass ratio */
  const interpolator_t * ALIGNED f;   /* Interpolator array */
  const grid_t         *         g;   /* Local domain grid parameters */
} uncenter_p_pipeline_args_t;

DECLARE_PIPELINE(uncenter_p);

/* In energy_p_pipeline.c */

typedef struct energy_p_pipeline_args {
  const particle_t     * ALIGNED p;   /* Particle array */
  int                            n;   /* Number of particles */
  float                          q_m; /* Charge to mass ratio */
  const interpolator_t * ALIGNED f;   /* Interpolator array */
  const grid_t         *         g;   /* Local domain grid parameters */
  double en[MAX_PIPELINE];            /* Return values */
} energy_p_pipeline_args_t;

DECLARE_PIPELINE(energy_p);

#undef DECLARE_PIPELINES

END_C_DECLS

#endif /* _particle_pipelines_h_ */
