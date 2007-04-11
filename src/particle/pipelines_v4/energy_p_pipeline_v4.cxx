#include <particle_pipelines.h>

#ifdef V4_ACCELERATION
using namespace v4;

void
energy_p_pipeline_v4( energy_p_pipeline_args_t * args,
                      int pipeline_rank ) {
  const particle_t     * ALIGNED p   = args->p;
  int                            nq  = args->n >> 2;
  const float                    q_m = args->q_m;
  const interpolator_t * ALIGNED f0  = args->f;
  const grid_t *                 g   = args->g;
  double n_target;

  v4float dx, dy, dz; v4int ii;
  v4float ex, ey, ez;
  v4float v0, v1, v2, q;
  v4float qdt_2mc(0.5*q_m*g->dt/g->cvac);
  v4float one(1.);
  float *vp0, *vp1, *vp2, *vp3;
  double en0 = 0, en1 = 0, en2 = 0, en3 = 0;

  // Determine which particle quads to process in this pipeline

  n_target = (double)nq / (double)n_pipeline;
  nq  = (int)( n_target*(double) pipeline_rank    + 0.5 );
  p  += 4*nq;
  nq  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - nq;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=4 ) {
    load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

    // Interpolate fields
    vp0 = (float *)(f0 + ii(0));
    vp1 = (float *)(f0 + ii(1));
    vp2 = (float *)(f0 + ii(2));
    vp3 = (float *)(f0 + ii(3));
    load_4x4_tr(vp0,  vp1,  vp2,  vp3,  ex,v0,v1,v2); ex = fma( fma( dy, v2, v1 ), dz, fma( dy, v0, ex ) );
    load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,ey,v0,v1,v2); ey = fma( fma( dz, v2, v1 ), dx, fma( dz, v0, ey ) );
    load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,ez,v0,v1,v2); ez = fma( fma( dx, v2, v1 ), dy, fma( dx, v0, ez ) );

    // Update momentum to half step
    // (note Boris rotation does not change energy so it is unnecessary)
    load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,v0,v1,v2,q);
    v0 = fma( ex, qdt_2mc, v0 );
    v1 = fma( ey, qdt_2mc, v1 );
    v2 = fma( ez, qdt_2mc, v2 );

    // Accumulate energy
    v0  = fma( v0,v0, fma( v1,v1, v2*v2 ) );
    v0 /= sqrt(one+v0)+one; 
    en0 += (double)v0(0)*(double)q(0);
    en1 += (double)v0(1)*(double)q(1);
    en2 += (double)v0(2)*(double)q(2);
    en3 += (double)v0(3)*(double)q(3);
  }

  args->en[pipeline_rank] = en0 + en1 + en2 + en3;
}

#endif
