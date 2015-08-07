#define IN_spa
#define HAS_V4_PIPELINE
#include "spa_private.h"

void
energy_p_pipeline( energy_p_pipeline_args_t * RESTRICT args,
                   int pipeline_rank,
                   int n_pipeline ) {
  const interpolator_t * RESTRICT ALIGNED(128) f = args->f;
  const particle_t     * RESTRICT ALIGNED(32)  p = args->p;
  const float qdt_2mc = args->qdt_2mc;
  const float msp     = args->msp;
  const float one     = 1;

  float dx, dy, dz;
  float v0, v1, v2;

  double en = 0;

  int i, n, n0, n1;

  // Determine which particles this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, n1 ); n1 += n0;

  // Process particles quads for this pipeline

  for( n=n0; n<n1; n++ ) {
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
    en += (double)v0;
  }

  args->en[pipeline_rank] = en;
}

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

using namespace v4;

void
energy_p_pipeline_v4( energy_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline ) {
  const interpolator_t * RESTRICT ALIGNED(128) f = args->f;
  const particle_t     * RESTRICT ALIGNED(128) p = args->p;

  const float          * RESTRICT ALIGNED(16)  vp0;
  const float          * RESTRICT ALIGNED(16)  vp1;
  const float          * RESTRICT ALIGNED(16)  vp2;
  const float          * RESTRICT ALIGNED(16)  vp3;

  const v4float qdt_2mc(args->qdt_2mc);
  const v4float msp(args->msp);
  const v4float one(1.);

  v4float dx, dy, dz;
  v4float ex, ey, ez;
  v4float v0, v1, v2, w;
  v4int i;

  double en0 = 0, en1 = 0, en2 = 0, en3 = 0;

  int n0, nq;

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, nq );
  p += n0;
  nq >>= 2;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=4 ) {
    load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,i);

    // Interpolate fields

    vp0 = (float *)(f + i(0));
    vp1 = (float *)(f + i(1));
    vp2 = (float *)(f + i(2));
    vp3 = (float *)(f + i(3));
    load_4x4_tr(vp0,  vp1,  vp2,  vp3,  ex,v0,v1,v2); ex = fma( fma( dy, v2, v1 ), dz, fma( dy, v0, ex ) );
    load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,ey,v0,v1,v2); ey = fma( fma( dz, v2, v1 ), dx, fma( dz, v0, ey ) );
    load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,ez,v0,v1,v2); ez = fma( fma( dx, v2, v1 ), dy, fma( dx, v0, ez ) );

    // Update momentum to half step
    // (note Boris rotation does not change energy so it is unnecessary)

    load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,v0,v1,v2,w);
    v0  = fma( ex, qdt_2mc, v0 );
    v1  = fma( ey, qdt_2mc, v1 );
    v2  = fma( ez, qdt_2mc, v2 );

    // Accumulate energy

    v0 = fma( v0,v0, fma( v1,v1, v2*v2 ) );
    v0 = (msp * w) * (v0 / (one + sqrt(one + v0))); 
    en0 += (double)v0(0);
    en1 += (double)v0(1);
    en2 += (double)v0(2);
    en3 += (double)v0(3);
  }

  args->en[pipeline_rank] = en0 + en1 + en2 + en3;
}

#endif

double
energy_p( const species_t            * RESTRICT sp,
          const interpolator_array_t * RESTRICT ia ) {
  DECLARE_ALIGNED_ARRAY( energy_p_pipeline_args_t, 128, args, 1 );
  DECLARE_ALIGNED_ARRAY( double, 128, en, MAX_PIPELINE+1 );
  double local, global;
  int rank;

  if( !sp || !ia || sp->g!=ia->g ) ERROR(( "Bad args" ));

  // Have the pipelines do the bulk of particles in quads and have the
  // host do the final incomplete quad.

  args->p       = sp->p;
  args->f       = ia->i;
  args->en      = en;
  args->qdt_2mc = (sp->q*sp->g->dt)/(2*sp->m*sp->g->cvac);
  args->msp     = sp->m;
  args->np      = sp->np;

  EXEC_PIPELINES( energy_p, args, 0 );
  WAIT_PIPELINES();

  local = 0; for( rank=0; rank<=N_PIPELINE; rank++ ) local += en[rank];
  mp_allsum_d( &local, &global, 1 );
  return global*((double)sp->g->cvac*(double)sp->g->cvac);
}
