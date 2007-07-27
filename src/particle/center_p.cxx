#define IN_particle_pipeline
#define HAS_V4_PIPELINE
#include <particle_pipelines.h>

static void
center_p_pipeline( center_p_pipeline_args_t * args,
                   int pipeline_rank,
                   int n_pipeline ) {
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  particle_t           * ALIGNED(32)  p;
  const interpolator_t * ALIGNED(16)  f;

  const float qdt_2mc        =     args->qdt_2mc;
  const float qdt_4mc        = 0.5*args->qdt_2mc; // For half Boris rotate
  const float one            = 1.;
  const float one_third      = 1./3.;
  const float two_fifteenths = 2./15.;

  float dx, dy, dz, ux, uy, uz;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4;

  int first, ii, n;

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, first, n );
  p = args->p0 + first;

  // Process particles for this pipeline

  for(;n;n--,p++) {
    dx   = p->dx;                            // Load position
    dy   = p->dy;
    dz   = p->dz;
    ii   = p->i;
    f    = f0 + ii;                          // Interpolate E
    hax  = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                     dz*( f->dexdz + dy*f->d2exdydz ) );
    hay  = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                     dx*( f->deydx + dz*f->d2eydzdx ) );
    haz  = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                     dy*( f->dezdy + dx*f->d2ezdxdy ) );
    cbx  = f->cbx + dx*f->dcbxdx;            // Interpolate B
    cby  = f->cby + dy*f->dcbydy;
    cbz  = f->cbz + dz*f->dcbzdz;
    ux   = p->ux;                            // Load momentum
    uy   = p->uy;
    uz   = p->uz;
    ux  += hax;                              // Half advance E
    uy  += hay;
    uz  += haz;
    v0   = qdt_4mc/(float)sqrt(one + (ux*ux + (uy*uy + uz*uz)));
    /**/                                     // Boris - scalars
    v1   = cbx*cbx + (cby*cby + cbz*cbz);
    v2   = (v0*v0)*v1;
    v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
    v4   = v3/(one+v1*(v3*v3));
    v4  += v4;
    v0   = ux + v3*( uy*cbz - uz*cby );      // Boris - uprime
    v1   = uy + v3*( uz*cbx - ux*cbz );
    v2   = uz + v3*( ux*cby - uy*cbx );
    ux  += v4*( v1*cbz - v2*cby );           // Boris - rotation
    uy  += v4*( v2*cbx - v0*cbz );
    uz  += v4*( v0*cby - v1*cbx );
    p->ux = ux;                              // Store momentum
    p->uy = uy;
    p->uz = uz;
  }
}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

#error "SPU version not hooked up yet!"

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

using namespace v4;

static void
center_p_pipeline_v4( center_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline ) {
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  particle_t           * ALIGNED(128) p;
  const float          * ALIGNED(16)  vp0;
  const float          * ALIGNED(16)  vp1;
  const float          * ALIGNED(16)  vp2;
  const float          * ALIGNED(16)  vp3;

  const v4float qdt_2mc(    args->qdt_2mc);
  const v4float qdt_4mc(0.5*args->qdt_2mc); // For half Boris rotate
  const v4float one(1.);
  const v4float one_third(1./3.);
  const v4float two_fifteenths(2./15.);

  v4float dx, dy, dz, ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v0, v1, v2, v3, v4, v5;
  v4int   ii;

  int itmp, nq;

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );
  p = args->p0 + itmp;
  nq >>= 2;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=4 ) {
    load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

    // Interpolate fields
    vp0 = (const float * ALIGNED(16))(f0 + ii(0));
    vp1 = (const float * ALIGNED(16))(f0 + ii(1));
    vp2 = (const float * ALIGNED(16))(f0 + ii(2));
    vp3 = (const float * ALIGNED(16))(f0 + ii(3));
    load_4x4_tr(vp0,  vp1,  vp2,  vp3,  hax,v0,v1,v2); hax = qdt_2mc*fma( fma( dy, v2, v1 ), dz, fma( dy, v0, hax ) );
    load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,hay,v3,v4,v5); hay = qdt_2mc*fma( fma( dz, v5, v4 ), dx, fma( dz, v3, hay ) );
    load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,haz,v0,v1,v2); haz = qdt_2mc*fma( fma( dx, v2, v1 ), dy, fma( dx, v0, haz ) );
    load_4x4_tr(vp0+12,vp1+12,vp2+12,vp3+12,cbx,v3,cby,v4); cbx = fma( v3, dx, cbx );
    /**/                                                    cby = fma( v4, dy, cby );
    load_4x2_tr(vp0+16,vp1+16,vp2+16,vp3+16,cbz,v5);        cbz = fma( v5, dz, cbz );

    // Update momentum
    load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,q);
    /**/                                             // Could use load_4x3_tr
    ux += hax;
    uy += hay;
    uz += haz;
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
    store_4x4_tr(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux);
    /**/                                             // Could use store_4x3_tr
  }
}

#endif

void
center_p( particle_t           * ALIGNED(128) p0,
          const int                           np,
          const float                         q_m,
          const interpolator_t * ALIGNED(128) f0,
          const grid_t         *              g ) {
  DECLARE_ALIGNED_ARRAY( center_p_pipeline_args_t, 128, args, 1 );

  // FIXME: p0 NULL checking
  if( np<0     ) ERROR(("Bad number of particles"));
  if( f0==NULL ) ERROR(("Bad interpolator"));
  if( g==NULL  ) ERROR(("Bad grid"));

  // Have the pipelines do the bulk of particles in quads and have the
  // host do the final incomplete quad.

  args->p0      = p0;
  args->f0      = f0;
  args->qdt_2mc = 0.5*q_m*g->dt/g->cvac;
  args->np      = np;

  EXEC_PIPELINES( center_p, args, 0 );
  WAIT_PIPELINES();
}
