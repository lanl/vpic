#include <particle.h>

#ifndef V4_ACCELERATION
#define CENTER_P_PIPELINE (pipeline_func_t)center_p_pipeline
#else
#define CENTER_P_PIPELINE (pipeline_func_t)center_p_pipeline_v4
#endif

typedef struct center_p_pipeline_args {
  particle_t           * ALIGNED(128) p0;  // Particle array
  int                                 np;  // Number of particles
  float                               q_m; // Charge to mass ratio
  const interpolator_t * ALIGNED(16)  f;   // Interpolator array
  const grid_t         *              g;   // Local domain grid parameters
} center_p_pipeline_args_t;

static void
center_p_pipeline( center_p_pipeline_args_t * args,
                   int pipeline_rank,
                   int n_pipeline ) {
  particle_t           * ALIGNED(16) p   = args->p0;
  int                                n   = args->np;
  const float                        q_m = args->q_m;
  const interpolator_t * ALIGNED(16) f0  = args->f;
  const grid_t *                     g   = args->g;

  const float qdt_2mc        = 0.5 *q_m*g->dt/g->cvac;
  const float qdt_4mc        = 0.25*q_m*g->dt/g->cvac; // For half Boris rotate
  const float one            = 1.;
  const float one_third      = 1./3.;
  const float two_fifteenths = 2./15.;

  const interpolator_t * ALIGNED(16) f;

  int ii;
  float dx, dy, dz, ux, uy, uz;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4;

  if( pipeline_rank==n_pipeline ) { // Host does left over cleanup

    // Determine which particles the host processes

    p += n;
    n &= 7;
    p -= n;

  } else { // Pipelines do any rough equal number of particle quads

    // Determine which particles to process in this pipeline

    double n_target = (double)(n>>3)/(double)n_pipeline;
    n  = (int)( n_target*(double) pipeline_rank    + 0.5 );
    p += 8*n;
    n  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - n;
    n *= 8;

  }

  // Process particles for this pipeline

  for(;n;n--,p++) {
    dx = p->dx;                              // Load position
    dy = p->dy;
    dz = p->dz;
    ii = p->i;
    f = f0 + ii;                             // Interpolate E
    hax = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                    dz*( f->dexdz + dy*f->d2exdydz ) );
    hay = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                    dx*( f->deydx + dz*f->d2eydzdx ) );
    haz = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                    dy*( f->dezdy + dx*f->d2ezdxdy ) );
    cbx = f->cbx + dx*f->dcbxdx;             // Interpolate B
    cby = f->cby + dy*f->dcbydy;
    cbz = f->cbz + dz*f->dcbzdz;
    ux = p->ux;                              // Load momentum
    uy = p->uy;
    uz = p->uz;
    ux += hax;                               // Half advance E
    uy += hay;
    uz += haz;
    v0 = qdt_4mc/(float)sqrt(one + (ux*ux + (uy*uy + uz*uz))); // Boris - scalars
    v1 = cbx*cbx + (cby*cby + cbz*cbz);
    v2 = (v0*v0)*v1;
    v3 = v0*(one+v2*(one_third+v2*two_fifteenths));
    v4 = v3/(one+v1*(v3*v3));
    v4 += v4;
    v0 = ux + v3*( uy*cbz - uz*cby );        // Boris - uprime
    v1 = uy + v3*( uz*cbx - ux*cbz );
    v2 = uz + v3*( ux*cby - uy*cbx );
    ux += v4*( v1*cbz - v2*cby );            // Boris - rotation
    uy += v4*( v2*cbx - v0*cbz );
    uz += v4*( v0*cby - v1*cbx );
    p->ux = ux;                              // Store momentum
    p->uy = uy;
    p->uz = uz;
  }
}

#ifdef V4_ACCELERATION

using namespace v4;

static void
center_p_pipeline_v4( center_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline ) {
  particle_t           * ALIGNED(64) p   = args->p0;
  int                                nq  = args->np;
  const float                        q_m = args->q_m;
  const interpolator_t * ALIGNED(16) f0  = args->f;
  const grid_t         *             g   = args->g;
  double n_target;

  v4float dx, dy, dz; v4int ii;
  v4float ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v0, v1, v2, v3, v4, v5;
  v4float qdt_2mc(0.5 *q_m*g->dt/g->cvac);
  v4float qdt_4mc(0.25*q_m*g->dt/g->cvac); // For half Boris rotate
  v4float one(1.);
  v4float one_third(1./3.);
  v4float two_fifteenths(2./15.);
  float *vp0, *vp1, *vp2, *vp3;

  // Determine which particles quads to process in this pipeline

  n_target = (double)(nq>>3) / (double)n_pipeline;
  nq  = (int)( n_target*(double) pipeline_rank    + 0.5 );
  p  += 8*nq;
  nq  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - nq;
  nq *= 2;

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
    store_4x4_tr(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux); // Could use store_4x3_tr
  }
}

#endif

void
center_p( particle_t           * ALIGNED(128) p0,
          const int                           np,
          const float                         q_m,
          const interpolator_t * ALIGNED(16)  f,
          const grid_t         *              g ) {
  center_p_pipeline_args_t args[1];

  // FIXME: p0 NULL checking
  if( np<0    ) ERROR(("Bad number of particles"));
  if( f==NULL ) ERROR(("Bad interpolator"));
  if( g==NULL ) ERROR(("Bad grid"));

  // Have the pipelines do the bulk of particles in quads and have the
  // host do the final incomplete quad.

  args->p0  = p0;
  args->np  = np;
  args->q_m = q_m;
  args->f   = f;
  args->g   = g;

  PSTYLE.dispatch( CENTER_P_PIPELINE, args, 0 );
  center_p_pipeline( args, PSTYLE.n_pipeline, PSTYLE.n_pipeline );
  PSTYLE.wait();
}
