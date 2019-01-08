#define IN_spa

#define HAS_V4_PIPELINE
#define HAS_V8_PIPELINE
#define HAS_V16_PIPELINE

#include "spa_private.h"

#include "../../../util/pipelines/pipelines_exec.h"

//----------------------------------------------------------------------------//
// Reference implementation for an energy_p pipeline function which does not
// make use of explicit calls to vector intrinsic functions.  This function
// calculates kinetic energy, normalized by c^2.
//----------------------------------------------------------------------------//

void
energy_p_pipeline_scalar( energy_p_pipeline_args_t * RESTRICT args,
                          int pipeline_rank,
                          int n_pipeline )
{
  const interpolator_t * RESTRICT ALIGNED(128) f = args->f;
  const particle_t     * RESTRICT ALIGNED(32)  p = args->p;

  const float qdt_2mc = args->qdt_2mc;
  const float msp     = args->msp;
  const float one     = 1.0;

  float dx, dy, dz;
  float v0, v1, v2;

  double en = 0.0;

  int i, n, n0, n1;

  // Determine which particles this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, n1 );

  n1 += n0;

  // Process particles quads for this pipeline.

  for( n = n0; n < n1; n++ )
  {
    dx  = p[n].dx;
    dy  = p[n].dy;
    dz  = p[n].dz;
    i   = p[n].i;

    v0  = p[n].ux + qdt_2mc*(    ( f[i].ex    + dy*f[i].dexdy    ) +
                              dz*( f[i].dexdz + dy*f[i].d2exdydz ) );

    v1  = p[n].uy + qdt_2mc*(    ( f[i].ey    + dz*f[i].deydz    ) +
                              dx*( f[i].deydx + dz*f[i].d2eydzdx ) );

    v2  = p[n].uz + qdt_2mc*(    ( f[i].ez    + dx*f[i].dezdx    ) +
                              dy*( f[i].dezdy + dx*f[i].d2ezdxdy ) );

    v0  = v0*v0 + v1*v1 + v2*v2;

    v0  = (msp * p[n].w) * (v0 / (one + sqrtf(one + v0)));

    en += ( double ) v0;
  }

  args->en[pipeline_rank] = en;
}

//----------------------------------------------------------------------------//
// Top level function to select and call the proper energy_p pipeline
// function.
//----------------------------------------------------------------------------//

double
energy_p_pipeline( const species_t * RESTRICT sp,
                   const interpolator_array_t * RESTRICT ia )
{
  DECLARE_ALIGNED_ARRAY( energy_p_pipeline_args_t, 128, args, 1 );

  DECLARE_ALIGNED_ARRAY( double, 128, en, MAX_PIPELINE+1 );

  double local, global;
  int rank;

  if ( !sp || !ia || sp->g != ia->g )
  {
    ERROR( ( "Bad args" ) );
  }

  // Have the pipelines do the bulk of particles in blocks and have the
  // host do the final incomplete block.

  args->p       = sp->p;
  args->f       = ia->i;
  args->en      = en;
  args->qdt_2mc = (sp->q*sp->g->dt)/(2*sp->m*sp->g->cvac);
  args->msp     = sp->m;
  args->np      = sp->np;

  EXEC_PIPELINES( energy_p, args, 0 );

  WAIT_PIPELINES();

  local = 0.0;
  for( rank = 0; rank <= N_PIPELINE; rank++ )
  {
    local += en[rank];
  }

  mp_allsum_d( &local, &global, 1 );

  return global * ( ( double ) sp->g->cvac *
		    ( double ) sp->g->cvac );
}
