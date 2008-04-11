// FIXME: PARTICLE MOVERS NEED TO BE OVERALLOCATED IN STRUCTORS TO
// ACCOUNT FOR SPLITTING THE MOVER ARRAY BETWEEN HOST AND PIPELINES

#define IN_particle_pipeline
#define HAS_V4_PIPELINE
#define HAS_SPU_PIPELINE
#include <particle_pipelines.h>

static void
advance_p_pipeline( advance_p_pipeline_args_t * args,
                    int pipeline_rank,
                    int n_pipeline ) {
  particle_t           * ALIGNED(128) p0 = args->p0;
  accumulator_t        * ALIGNED(128) a0 = args->a0;
  const interpolator_t * ALIGNED(128) f0 = args->f0;
  const grid_t *                      g  = args->g;

  particle_t           * ALIGNED(32)  p;
  particle_mover_t     * ALIGNED(16)  pm;
  const interpolator_t * ALIGNED(16)  f;
  float                * ALIGNED(16)  a;

  const float qdt_2mc        = args->qdt_2mc;
  const float cdt_dx         = args->cdt_dx;
  const float cdt_dy         = args->cdt_dy;
  const float cdt_dz         = args->cdt_dz;
  const float one            = 1.;
  const float one_third      = 1./3.;
  const float two_fifteenths = 2./15.;

  float dx, dy, dz, ux, uy, uz, q;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4, v5;

  int itmp, ii, n, nm, max_nm;
  
  particle_mover_t local_pm[1];

  // Determine which quads of particles quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, n );
  p = args->p0 + itmp;

  // Determine which movers are reserved for this pipeline
  // Movers (16 bytes) should be reserved for pipelines in at least
  // multiples of 8 such that the set of particle movers reserved for
  // a pipeline is 128-byte aligned and a multiple of 128-byte in
  // size.  The host is guaranteed to get enough movers to process its
  // particles with this allocation.

  max_nm = args->max_nm - (args->np&15);
  if( max_nm<0 ) max_nm = 0;
  DISTRIBUTE( max_nm, 8, pipeline_rank, n_pipeline, itmp, max_nm );
  if( pipeline_rank==n_pipeline ) max_nm = args->max_nm - itmp;
  pm   = args->pm + itmp;
  nm   = 0;
  itmp = 0;

  // Determine which accumulator array to use
  // The host gets the first accumulator array

  if( pipeline_rank!=n_pipeline )
    a0 += (1+pipeline_rank)*
          POW2_CEIL((args->nx+2)*(args->ny+2)*(args->nz+2),2);

  // Process particles for this pipeline

  for(;n;n--,p++) {
    dx   = p->dx;                             // Load position
    dy   = p->dy;
    dz   = p->dz;
    ii   = p->i;
    f    = f0 + ii;                           // Interpolate E
    hax  = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                     dz*( f->dexdz + dy*f->d2exdydz ) );
    hay  = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                     dx*( f->deydx + dz*f->d2eydzdx ) );
    haz  = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                     dy*( f->dezdy + dx*f->d2ezdxdy ) );
    cbx  = f->cbx + dx*f->dcbxdx;             // Interpolate B
    cby  = f->cby + dy*f->dcbydy;
    cbz  = f->cbz + dz*f->dcbzdz;
    ux   = p->ux;                             // Load momentum
    uy   = p->uy;
    uz   = p->uz;
    q    = p->q;
    ux  += hax;                               // Half advance E
    uy  += hay;
    uz  += haz;
    v0   = qdt_2mc/sqrtf(one + (ux*ux + (uy*uy + uz*uz)));
    /**/                                      // Boris - scalars
    v1   = cbx*cbx + (cby*cby + cbz*cbz);
    v2   = (v0*v0)*v1;
    v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
    v4   = v3/(one+v1*(v3*v3));
    v4  += v4;
    v0   = ux + v3*( uy*cbz - uz*cby );       // Boris - uprime
    v1   = uy + v3*( uz*cbx - ux*cbz );
    v2   = uz + v3*( ux*cby - uy*cbx );
    ux  += v4*( v1*cbz - v2*cby );            // Boris - rotation
    uy  += v4*( v2*cbx - v0*cbz );
    uz  += v4*( v0*cby - v1*cbx );
    ux  += hax;                               // Half advance E
    uy  += hay;
    uz  += haz;
    p->ux = ux;                               // Store momentum
    p->uy = uy;
    p->uz = uz;
    v0   = one/sqrtf(one + (ux*ux+ (uy*uy + uz*uz)));
    /**/                                      // Get norm displacement
    ux  *= cdt_dx;
    uy  *= cdt_dy;
    uz  *= cdt_dz;
    ux  *= v0;
    uy  *= v0;
    uz  *= v0;
    v0   = dx + ux;                           // Streak midpoint (inbnds)
    v1   = dy + uy;
    v2   = dz + uz;
    v3   = v0 + ux;                           // New position
    v4   = v1 + uy;
    v5   = v2 + uz;

    if(  v3<=one &&  v4<=one &&  v5<=one &&   // Check if inbnds
        -v3<=one && -v4<=one && -v5<=one ) {

      // Common case (inbnds).  Note: accumulator values are 4 times
      // the total physical charge that passed through the appropriate
      // current quadrant in a time-step

      p->dx = v3;                             // Store new position
      p->dy = v4;
      p->dz = v5;
      dx = v0;                                // Streak midpoint
      dy = v1;
      dz = v2;
      v5 = q*ux*uy*uz*one_third;              // Compute correction
      a  = (float *)( a0 + ii );              // Get accumulator

#     define ACCUMULATE_J(X,Y,Z,offset)                                 \
      v4  = q*u##X;   /* v2 = q ux                            */        \
      v1  = v4*d##Y;  /* v1 = q ux dy                         */        \
      v0  = v4-v1;    /* v0 = q ux (1-dy)                     */        \
      v1 += v4;       /* v1 = q ux (1+dy)                     */        \
      v4  = one+d##Z; /* v4 = 1+dz                            */        \
      v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */        \
      v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */        \
      v4  = one-d##Z; /* v4 = 1-dz                            */        \
      v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */        \
      v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */        \
      v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */        \
      v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */        \
      v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */        \
      v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */        \
      a[offset+0] += v0;                                                \
      a[offset+1] += v1;                                                \
      a[offset+2] += v2;                                                \
      a[offset+3] += v3

      ACCUMULATE_J( x,y,z, 0 );
      ACCUMULATE_J( y,z,x, 4 );
      ACCUMULATE_J( z,x,y, 8 );

#     undef ACCUMULATE_J

    } else {                                    // Unlikely
      local_pm->dispx = ux;
      local_pm->dispy = uy;
      local_pm->dispz = uz;
      local_pm->i     = p - p0;
      if( move_p( p0, local_pm, a0, g ) ) {     // Unlikely
        if( nm<max_nm ) pm[nm++] = local_pm[0];
        else            itmp++;                 // Unlikely
      }
    }

  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

extern spe_program_handle_t advance_p_pipeline_spu;

#if 0
static void
advance_p_pipeline_signal() {
    spu.signal(READ_ARGS_AND_ADVANCE_PARTICLES);
}
#endif

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

using namespace v4;

// Note: The cell pipeline will need to have a SPU version of the
// move_p function for it as the normal move_p function will exist
// only on the PPU.

static void
advance_p_pipeline_v4( advance_p_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {
  particle_t           * ALIGNED(128) p0 = args->p0;
  accumulator_t        * ALIGNED(128) a0 = args->a0;
  const interpolator_t * ALIGNED(128) f0 = args->f0;
  const grid_t         *              g  = args->g;

  particle_t           * ALIGNED(128) p;
  particle_mover_t     * ALIGNED(16)  pm; 
  float                * ALIGNED(16)  vp0;
  float                * ALIGNED(16)  vp1;
  float                * ALIGNED(16)  vp2;
  float                * ALIGNED(16)  vp3;

  const v4float qdt_2mc(args->qdt_2mc);
  const v4float cdt_dx(args->cdt_dx);
  const v4float cdt_dy(args->cdt_dy);
  const v4float cdt_dz(args->cdt_dz);
  const v4float one(1.);
  const v4float one_third(1./3.);
  const v4float two_fifteenths(2./15.);
  const v4float neg_one(-1.);

  v4float dx, dy, dz, ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v0, v1, v2, v3, v4, v5;
  v4int   ii, outbnd;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );
  p = args->p0 + itmp;
  nq>>=2;

  // Determine which movers are reserved for this pipeline.
  // Movers (16 bytes) should be reserved for pipelines in at least
  // multiples of 8 such that the set of particle movers reserved for
  // a pipeline is 128-byte aligned and a multiple of 128-byte in
  // size.  The host is guaranteed to get enough movers to process its
  // particles with this allocation.

  max_nm = args->max_nm - (args->np&15);
  if( max_nm<0 ) max_nm = 0;
  DISTRIBUTE( max_nm, 8, pipeline_rank, n_pipeline, itmp, max_nm );
  if( pipeline_rank==n_pipeline ) max_nm = args->max_nm - itmp;
  pm   = args->pm + itmp;
  nm   = 0;
  itmp = 0;

  // Determine which accumulator array to use
  // The host gets the first accumulator array

  a0 += (1+pipeline_rank)*
        POW2_CEIL((args->nx+2)*(args->ny+2)*(args->nz+2),2);

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=4 ) {
    load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

    // Interpolate fields
    vp0 = (float * ALIGNED(16))(f0 + ii(0));
    vp1 = (float * ALIGNED(16))(f0 + ii(1));
    vp2 = (float * ALIGNED(16))(f0 + ii(2));
    vp3 = (float * ALIGNED(16))(f0 + ii(3));
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
    vp0 = (float * ALIGNED(16))(a0 + ii(0)); // Accumulator pointers
    vp1 = (float * ALIGNED(16))(a0 + ii(1));
    vp2 = (float * ALIGNED(16))(a0 + ii(2));
    vp3 = (float * ALIGNED(16))(a0 + ii(3));

#   define ACCUMULATE_J(X,Y,Z,offset)                               \
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
    increment_4x1(vp0+offset,v0);                                   \
    increment_4x1(vp1+offset,v1);                                   \
    increment_4x1(vp2+offset,v2);                                   \
    increment_4x1(vp3+offset,v3)

    ACCUMULATE_J( x,y,z, 0 );
    ACCUMULATE_J( y,z,x, 4 );
    ACCUMULATE_J( z,x,y, 8 );

#   undef ACCUMULATE_J

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
    if( outbnd(N) ) {                       /* Unlikely */              \
      local_pm->dispx = ux(N);                                          \
      local_pm->dispy = uy(N);                                          \
      local_pm->dispz = uz(N);                                          \
      local_pm->i     = (p - p0) + N;                                   \
      if( move_p( p0, local_pm, a0, g ) ) { /* Unlikely */              \
        if( nm<max_nm ) copy_4x1( &pm[nm++], local_pm );                \
        else            itmp++;             /* Unlikely */              \
      }                                                                 \
    }

    MOVE_OUTBND(0);
    MOVE_OUTBND(1);
    MOVE_OUTBND(2);
    MOVE_OUTBND(3);

#   undef MOVE_OUTBND

  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}

#endif

DECLARE_ALIGNED_ARRAY( advance_p_pipeline_args_t, 128, args, 1 );
DECLARE_ALIGNED_ARRAY( particle_mover_seg_t, 128, seg, MAX_PIPELINE+1 );

void advance_p_initialize() {
  INIT_PIPELINES( advance_p, args, 0 );
}

void advance_p_finalize() {
  FINALIZE_PIPELINES();
}

int
advance_p( particle_t           * ALIGNED(128) p0,
           const int                           np,
           const float                         q_m,
           particle_mover_t     * ALIGNED(128) pm,
           int                                 max_nm,       
           accumulator_t        * ALIGNED(128) a0,
           const interpolator_t * ALIGNED(128) f0,
           const grid_t         *              g ) {
  int nm, rank;

  if( p0==NULL ) ERROR(("Bad particle array"));
  if( np<0     ) ERROR(("Bad number of particles"));
  if( pm==NULL ) ERROR(("Bad particle mover"));
  if( max_nm<0 ) ERROR(("Bad number of movers"));
  if( a0==NULL ) ERROR(("Bad accumulator"));
  if( f0==NULL ) ERROR(("Bad interpolator"));
  if( g==NULL  ) ERROR(("Bad grid"));

  args->p0       = p0;
  args->pm       = pm;
  args->a0       = a0;
  args->f0       = f0;
  args->seg      = seg;
  args->g        = g;

  args->qdt_2mc  = 0.5*q_m*g->dt/g->cvac;
  args->cdt_dx   = g->cvac*g->dt/g->dx;
  args->cdt_dy   = g->cvac*g->dt/g->dy;
  args->cdt_dz   = g->cvac*g->dt/g->dz;

  args->np       = np;
  args->max_nm   = max_nm;
  args->nx       = g->nx;
  args->ny       = g->ny;
  args->nz       = g->nz;

# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
     defined(HAS_SPU_PIPELINE)
  args->neighbor = g->neighbor;
  args->rangel   = g->rangel;
  args->rangeh   = g->rangeh;
# endif

  // Have the host processor do the last incomplete bundle if necessary.
  // Note: This is overlapped with the pipelined processing.  As such,
  // it uses an entire accumulator.  Reserving an entirely accumulator
  // for the host processor to handle at most 15 particles is wasteful
  // of memory.  It is anticipated that it may be useful at some point
  // in the future have pipelines accumulating currents while the host
  // processor is doing other more substantive work (e.g. accumulating
  // currents from particles received from neighboring nodes).
  // However, it is worth reconsidering this at some point in the
  // future.

  EXEC_PIPELINES( advance_p, args, 0 );
  WAIT_PIPELINES();

  // FIXME: HIDEOUS HACK UNTIL BETTER PARTICLE MOVER SEMANTICS
  // INSTALLED FOR DEALING WITH PIPELINES.  COMPACT THE PARTICLE
  // MOVERS TO ELIMINATE HOLES FROM THE PIPELINING.

  nm = 0;
  for( rank=0; rank<=N_PIPELINE; rank++ ) {
    if( args->seg[rank].n_ignored )
      WARNING(( "Pipeline %i ran out of storage for %i movers",
                rank, args->seg[rank].n_ignored ));
    if( pm+nm!=args->seg[rank].pm )
      memmove( pm+nm, args->seg[rank].pm, args->seg[rank].nm*sizeof(*pm) );
    nm += args->seg[rank].nm;
  }
  
  return nm;
}
