#include <particle_pipelines.h>

#ifdef V4_ACCELERATION
using namespace v4;

void
uncenter_p_pipeline_v4( uncenter_p_pipeline_args_t * args,
                        int pipeline_rank,
                        int n_pipeline ) {
  particle_t           * ALIGNED p   = args->p;
  int                            nq  = args->n >> 2;
  const float                    q_m = args->q_m;
  const interpolator_t * ALIGNED f0  = args->f;
  const grid_t         *         g   = args->g;
  double n_target;

  v4float dx, dy, dz; v4int ii;
  v4float ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v0, v1, v2, v3, v4, v5;
  v4float qdt_2mc(-0.5 *q_m*g->dt/g->cvac); // Negative for backward half advance
  v4float qdt_4mc(-0.25*q_m*g->dt/g->cvac); // Negative for backward half Boris rotate
  v4float one(1.);
  v4float one_third(1./3.);
  v4float two_fifteenths(2./15.);
  float *vp0, *vp1, *vp2, *vp3;

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
    load_4x4_tr(vp0,  vp1,  vp2,  vp3,  hax,v0,v1,v2); hax = qdt_2mc*fma( fma( dy, v2, v1 ), dz, fma( dy, v0, hax ) );
    load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,hay,v3,v4,v5); hay = qdt_2mc*fma( fma( dz, v5, v4 ), dx, fma( dz, v3, hay ) );
    load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,haz,v0,v1,v2); haz = qdt_2mc*fma( fma( dx, v2, v1 ), dy, fma( dx, v0, haz ) );
    load_4x4_tr(vp0+12,vp1+12,vp2+12,vp3+12,cbx,v3,cby,v4); cbx = fma( v3, dx, cbx );
    /**/                                                    cby = fma( v4, dy, cby );
    load_4x2_tr(vp0+16,vp1+16,vp2+16,vp3+16,cbz,v5);        cbz = fma( v5, dz, cbz );

    // Update momentum
    load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,q); // Could use load_4x3_tr
    v0  = qdt_4mc*rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
    v1  = fma( cbx,cbx, fma( cby,cby, cbz*cbz ) );
    v2  = (v0*v0)*v1;
    v3  = v0*fma( v2, fma( v2, two_fifteenths, one_third ), one );
    v4  = v3*rcp( fma( v3*v3, v1, one ) ); v4 += v4;
    v0  = fma( fms( uy,cbz, uz*cby ), v3, ux );
    v1  = fma( fms( uz,cbx, ux*cbz ), v3, uy );
    v2  = fma( fms( ux,cby, uy*cbx ), v3, uz );
    ux  = fma( fms( v1,cbz, v2*cby ), v4, ux );
    uy  = fma( fms( v2,cbx, v0*cbz ), v4, uy );
    uz  = fma( fms( v0,cby, v1*cbx ), v4, uz );
    ux += hax;
    uy += hay;
    uz += haz;
    store_4x4_tr(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux); // Could use store_4x3_tr
  }
}

#endif
