/* FIXME: PARTICLE MOVERS NEED TO BE OVERALLOCATED IN STRUCTORS TO
   ACCOUNT FOR SPLITTING THE MOVER ARRAY BETWEEN HOST AND PIPELINES */

#include <particle.h>

#ifndef V4_ACCELERATION
#define ADVANCE_P_PIPELINE (pipeline_func_t)advance_p_pipeline
#else
#define ADVANCE_P_PIPELINE (pipeline_func_t)advance_p_pipeline_v4
#endif

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

static void
advance_p_pipeline( advance_p_pipeline_args_t * args,
                    int pipeline_rank,
                    int n_pipeline ) {
  particle_t           * ALIGNED p0  = args->p;
  int                            n   = args->n;
  const float                    q_m = args->q_m;
  particle_mover_t     * ALIGNED pm  = args->pm;
  int                            nm  = args->nm;
  accumulator_t        * ALIGNED a0  = args->a;
  const interpolator_t * ALIGNED f0  = args->f;
  const grid_t *                 g   = args->g;

  const float qdt_2mc        = 0.5*q_m*g->dt/g->cvac;
  const float cdt_dx         = g->cvac*g->dt/g->dx;
  const float cdt_dy         = g->cvac*g->dt/g->dy;
  const float cdt_dz         = g->cvac*g->dt/g->dz;
  const float one            = 1.;
  const float one_third      = 1./3.;
  const float two_fifteenths = 2./15.;

  particle_t           * ALIGNED p;
  const interpolator_t * ALIGNED f;
  float                * ALIGNED a;

  int ii;
  float dx, dy, dz, ux, uy, uz, q;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4, v5;

  if( pipeline_rank==n_pipeline ) { /* Host does left over cleanup */

    /* Determine which particles the host processes, which movers are
       reserved for the host.  Note the host uses the first
       accumulator array. */

    p  = p0 + n; pm += nm;
    n &= 3;      nm  = n>nm ? nm : n;
    p -= n;      pm -= nm;

  } else { /* Pipelines do any rough equal number of particle quads */

    double n_target;

    /* Determine which particles to process in this pipeline */

    n_target = (double)(n>>2)/(double)n_pipeline;
    n  = (int)( n_target*(double) pipeline_rank    + 0.5 );
    p  = p0 + 4*n;
    n  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - n;
    n *= 4;

    /* Determine which movers are reserved for this pipeline */

    nm -= (args->n&3)>nm ? nm : (args->n&3); /* Reserve last movers for host */

    n_target = (double)nm / (double)n_pipeline; 
    nm  = (int)( n_target*(double) pipeline_rank    + 0.5 );
    pm += nm;
    nm  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - nm;

    /* Determine which accumulator array to use */

    a0 += (1+pipeline_rank)*(g->nx+2)*(g->ny+2)*(g->nz+2);

  }

  args->seg[pipeline_rank].pm = pm;
  args->seg[pipeline_rank].nm = nm;

  /* Process particles quads for this pipeline */

  for(;n;n--,p++) {
    dx = p->dx;                              /* Load position */
    dy = p->dy;
    dz = p->dz;
    ii = p->i;
    f = f0 + ii;                             /* Interpolate E */
    hax = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                    dz*( f->dexdz + dy*f->d2exdydz ) );
    hay = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                    dx*( f->deydx + dz*f->d2eydzdx ) );
    haz = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                    dy*( f->dezdy + dx*f->d2ezdxdy ) );
    cbx = f->cbx + dx*f->dcbxdx;             /* Interpolate B */
    cby = f->cby + dy*f->dcbydy;
    cbz = f->cbz + dz*f->dcbzdz;
    ux = p->ux;                              /* Load momentum */
    uy = p->uy;
    uz = p->uz;
    q  = p->q;
    ux += hax;                               /* Half advance E */
    uy += hay;
    uz += haz;
    v0 = qdt_2mc/(float)sqrt(one + (ux*ux + (uy*uy + uz*uz))); /* Boris - scalars */
    v1 = cbx*cbx + (cby*cby + cbz*cbz);
    v2 = (v0*v0)*v1;
    v3 = v0*(one+v2*(one_third+v2*two_fifteenths));
    v4 = v3/(one+v1*(v3*v3));
    v4 += v4;
    v0 = ux + v3*( uy*cbz - uz*cby );        /* Boris - uprime */
    v1 = uy + v3*( uz*cbx - ux*cbz );
    v2 = uz + v3*( ux*cby - uy*cbx );
    ux += v4*( v1*cbz - v2*cby );            /* Boris - rotation */
    uy += v4*( v2*cbx - v0*cbz );
    uz += v4*( v0*cby - v1*cbx );
    ux += hax;                               /* Half advance E */
    uy += hay;
    uz += haz;
    p->ux = ux;                              /* Store momentum */
    p->uy = uy;
    p->uz = uz;
    v0 = one/(float)sqrt(one + (ux*ux+ (uy*uy + uz*uz))); /* Get norm displacement */
    ux *= cdt_dx;
    uy *= cdt_dy;
    uz *= cdt_dz;
    ux *= v0;
    uy *= v0;
    uz *= v0;
    v0 = dx + ux;                            /* Streak midpoint (inbnds) */
    v1 = dy + uy;
    v2 = dz + uz;
    v3 = v0 + ux;                            /* New position */
    v4 = v1 + uy;
    v5 = v2 + uz;
    if(  v3<=one && v4<=one && v5<=one &&          /* Check if inbnds */
        -v3<=one && -v4<=one && -v5<=one ) {

      /* Common case (inbnds)
         Note: accumulator values are 4 times the total physical charge that
         passed through the appropriate current quadrant in a time-step */

      p->dx = v3;                            /* Store new position */
      p->dy = v4;
      p->dz = v5;
      dx = v0;                               /* Streak mid */
      dy = v1;
      dz = v2;
      v5 = q*ux*uy*uz*one_third;             /* Compute correction */
      a = (float *)( a0 + ii );              /* Get accumulator */
#     define accumulate_j(X,Y,Z)                                  \
      v4  = q*u##X;   /* v2 = q ux                            */  \
      v1  = v4*d##Y;  /* v1 = q ux dy                         */  \
      v0  = v4-v1;    /* v0 = q ux (1-dy)                     */  \
      v1 += v4;       /* v1 = q ux (1+dy)                     */  \
      v4  = one+d##Z; /* v4 = 1+dz                            */  \
      v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */  \
      v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */  \
      v4  = one-d##Z; /* v4 = 1-dz                            */  \
      v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */  \
      v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */  \
      v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */  \
      v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */  \
      v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */  \
      v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */  \
      a[0] += v0;                                                 \
      a[1] += v1;                                                 \
      a[2] += v2;                                                 \
      a[3] += v3
      accumulate_j(x,y,z); a += 4;
      accumulate_j(y,z,x); a += 4;
      accumulate_j(z,x,y);
#     undef accumulate_j

    } else if( nm>0 ) {

      pm->dispx = ux;
      pm->dispy = uy;
      pm->dispz = uz;
      pm->i     = p - p0;
      if( move_p( p0, pm, a0, g ) ) pm++, nm--;

    }
  }

  args->seg[pipeline_rank].nm -= nm;
}

#ifdef V4_ACCELERATION

using namespace v4;

// Note: The cell pipeline will need to have a SPU version of the
// move_p function for it as the normal move_p function will exist
// only on the PPU.

static void
advance_p_pipeline_v4( advance_p_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {
  double n_target;

  particle_t           * ALIGNED p   = args->p;
  particle_t           * ALIGNED p0  = p;
  int                            nq  = args->n >> 2;
  const float                    q_m = args->q_m;
  particle_mover_t     * ALIGNED pm  = args->pm;
  int                            nm  = args->nm;
  accumulator_t        * ALIGNED a0  = args->a;
  const interpolator_t * ALIGNED f0  = args->f;
  const grid_t *                 g   = args->g;

  v4float dx, dy, dz; v4int ii;
  v4float ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v0, v1, v2, v3, v4, v5;
  v4int outbnd;
  v4float qdt_2mc(0.5*q_m*g->dt/g->cvac);
  v4float cdt_dx(g->cvac*g->dt/g->dx);
  v4float cdt_dy(g->cvac*g->dt/g->dy);
  v4float cdt_dz(g->cvac*g->dt/g->dz);
  v4float one(1.);
  v4float one_third(1./3.);
  v4float two_fifteenths(2./15.);
  v4float neg_one(-1.);
  float *vp0, *vp1, *vp2, *vp3;

  // Determine which particles to process in this pipeline

  n_target = (double)nq / (double)n_pipeline;
  nq  = (int)( n_target*(double) pipeline_rank    + 0.5 );
  p  += 4*nq;
  nq  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - nq;

  // Determine which movers are reserved for this pipeline

  nm -= (args->n&3)>nm ? nm : (args->n&3); // Reserve last movers for host

  n_target = (double)nm / (double)n_pipeline; 
  nm  = (int)( n_target*(double) pipeline_rank    + 0.5 );
  pm += nm;
  nm  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - nm;

  args->seg[pipeline_rank].pm = pm;
  args->seg[pipeline_rank].nm = nm;

  // Determine which accumulator array is reserved for this pipeline

  a0 += (1+pipeline_rank)*(g->nx+2)*(g->ny+2)*(g->nz+2);

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=4 ) {
    load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

    // Interpolate fields
    vp0 = (float *)(f0 + ii(0));
    vp1 = (float *)(f0 + ii(1));
    vp2 = (float *)(f0 + ii(2));
    vp3 = (float *)(f0 + ii(3));
    load_4x4_tr(vp0,  vp1,  vp2,  vp3,  hax,v0,v1,v2); hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );
    load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,hay,v3,v4,v5); hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );
    load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,haz,v0,v1,v2); haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );
    load_4x4_tr(vp0+12,vp1+12,vp2+12,vp3+12,cbx,v3,cby,v4); cbx = fma( v3, dx, cbx );
    /**/                                                    cby = fma( v4, dy, cby );
    load_4x2_tr(vp0+16,vp1+16,vp2+16,vp3+16,cbz,v5);        cbz = fma( v5, dz, cbz );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,q);
    ux += hax;
    uy += hay;
    uz += haz;
    v0  = qdt_2mc*rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
    v1  = fma( cbx,cbx, fma( cby,cby, cbz*cbz ) );
    v2  = (v0*v0)*v1;
    v3  = v0*fma( fma( two_fifteenths, v2, one_third ), v2, one );
    v4  = v3*rcp(fma( v3*v3, v1, one ));
    v4 += v4;
    v0  = fma( fms( uy,cbz, uz*cby ), v3, ux );
    v1  = fma( fms( uz,cbx, ux*cbz ), v3, uy );
    v2  = fma( fms( ux,cby, uy*cbx ), v3, uz );
    ux  = fma( fms( v1,cbz, v2*cby ), v4, ux );
    uy  = fma( fms( v2,cbx, v0*cbz ), v4, uy );
    uz  = fma( fms( v0,cby, v1*cbx ), v4, uz );
    ux += hax;
    uy += hay;
    uz += haz;
    store_4x4_tr(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux);
    
    // Update the position of inbnd particles
    v0  = rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
    ux *= cdt_dx;
    uy *= cdt_dy;
    uz *= cdt_dz;
    ux *= v0;
    uy *= v0;
    uz *= v0;      // ux,uy,uz are normalized displ (relative to cell size)
    v0  = dx + ux;
    v1  = dy + uy;
    v2  = dz + uz; // New particle midpoint
    v3  = v0 + ux;
    v4  = v1 + uy;
    v5  = v2 + uz; // New particle position
    outbnd = (v3>one) | (v3<neg_one) |
             (v4>one) | (v4<neg_one) |
             (v5>one) | (v5<neg_one);
    v3  = merge(outbnd,dx,v3); // Do not update outbnd particles
    v4  = merge(outbnd,dy,v4);
    v5  = merge(outbnd,dz,v5);
    store_4x4_tr(v3,v4,v5,ii,&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx);
    
    // Accumulate current of inbnd particles
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step
    q  = czero(outbnd,q);          // Do not accumulate outbnd particles
    dx = v0;                       // Streak midpoint (valid for inbnd only)
    dy = v1;
    dz = v2;
    v5 = q*ux*uy*uz*one_third;     // Charge conservation correction
    vp0 = (float *)(a0 + ii(0));   // Accumulator pointers
    vp1 = (float *)(a0 + ii(1));
    vp2 = (float *)(a0 + ii(2));
    vp3 = (float *)(a0 + ii(3));
#   define accumulate_j(X,Y,Z)                                      \
    v4  = q*u##X;   /* v4 = q ux                            */      \
    v1  = v4*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v4-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v4;       /* v1 = q ux (1+dy)                     */      \
    v4  = one+d##Z; /* v4 = 1+dz                            */      \
    v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v4  = one-d##Z; /* v4 = 1-dz                            */      \
    v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */      \
    transpose(v0,v1,v2,v3);                                         \
    increment_4x1(vp0,v0); vp0+=4;                                  \
    increment_4x1(vp1,v1); vp1+=4;                                  \
    increment_4x1(vp2,v2); vp2+=4;                                  \
    increment_4x1(vp3,v3); vp3+=4
    accumulate_j(x,y,z);
    accumulate_j(y,z,x);
    accumulate_j(z,x,y);
#   undef accumulate_j

    // Update position and accumulate outbnd

#   define move_outbnd(N)                             \
    if( outbnd(N) && nm>0 ) {                         \
      pm->dispx = ux(N);                              \
      pm->dispy = uy(N);                              \
      pm->dispz = uz(N);                              \
      pm->i     = (p - p0) + N;                       \
      if( move_p( p0, pm, a0, g ) ) pm++, nm--;       \
    }
    move_outbnd(0);
    move_outbnd(1);
    move_outbnd(2);
    move_outbnd(3);
#   undef move_outbnd
  }

  args->seg[pipeline_rank].nm -= nm;
}

#endif

int
advance_p( particle_t           * ALIGNED p,
           const int                      n,
           const float                    q_m,
           particle_mover_t     * ALIGNED pm,
           int                            nm,       
           accumulator_t        * ALIGNED a,
           const interpolator_t * ALIGNED f,
           const grid_t         *         g ) {
  advance_p_pipeline_args_t args[1];
  int rank;

  if( p==NULL  ) ERROR(("Bad particle array"));
  if( n<0      ) ERROR(("Bad number of particles"));
  if( pm==NULL ) ERROR(("Bad particle mover"));
  if( nm<0     ) ERROR(("Bad number of movers"));
  if( a==NULL  ) ERROR(("Bad accumulator"));
  if( f==NULL  ) ERROR(("Bad interpolator"));
  if( g==NULL  ) ERROR(("Bad grid"));

  args->p   = p;
  args->n   = n;
  args->q_m = q_m;
  args->pm  = pm;
  args->nm  = nm;
  args->a   = a;
  args->f   = f;
  args->g   = g;

  PMETHOD.dispatch( ADVANCE_P_PIPELINE, args, 0 );

  /* Have the host processor do the incomplete quad if necessary.
     Note: This is overlapped with the pipelined processing.  As such,
     it uses an entire accumulator.  Reserving an entirely accumulator
     for the host processor to handle at most 3 particles is wasteful
     of memory.  It is anticipated that it may be useful at some point
     in the future have pipelines accumulating currents while the host
     processor is doing other more substantive work (e.g. accumulating
     currents from particles received from neighboring nodes).
     However, it is worth reconsidering this at some point in the
     future. */

  advance_p_pipeline( args, PMETHOD.n_pipeline, PMETHOD.n_pipeline );

  PMETHOD.wait();

  /* FIXME: HIDEOUS HACK UNTIL BETTER PARTICLE MOVER SEMANTICS
     INSTALLED FOR DEALING WITH PIPELINES.  COMPACT THE PARTICLE
     MOVERS TO ELIMINATE HOLES IN THE ALLOCATION. */

  nm = 0;
  for( rank=0; rank<=PMETHOD.n_pipeline; rank++ ) {
    if( pm+nm!=args->seg[rank].pm )
      memmove( pm+nm, args->seg[rank].pm, args->seg[rank].nm*sizeof(*pm) );
    nm += args->seg[rank].nm;
  }

  return nm;
}
