using namespace v4;

// #if 0
// Method 1.
void
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
  const v4float qsp(args->qsp);
  const v4float one(1.0f);
  const v4float one_third(1.0f/3.0f);
  const v4float two_fifteenths(2.0f/15.0f);
  const v4float neg_one(-1.0f);

  const float _qsp = args->qsp;

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

    load_4x4_tr(vp0,  vp1,  vp2,  vp3,  hax,v0,v1,v2);
    hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );

    load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,hay,v3,v4,v5);
    hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );

    load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,haz,v0,v1,v2);
    haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );

    load_4x4_tr(vp0+12,vp1+12,vp2+12,vp3+12,cbx,v3,cby,v4);
    cbx = fma( v3, dx, cbx );
    cby = fma( v4, dy, cby );

    load_4x2_tr(vp0+16,vp1+16,vp2+16,vp3+16,cbz,v5);
    cbz = fma( v5, dz, cbz );

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
    q  = czero(outbnd,q*qsp);       // Do not accumulate outbnd particles
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
      if( move_p( p0, local_pm, a0, g, _qsp ) ) { /* Unlikely */        \
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
// #endif // Method 1.

#if 0
// Method 2.
void
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
  const v4float qsp(args->qsp);
  const v4float one(1.);
  const v4float one_third(1./3.);
  const v4float two_fifteenths(2./15.);
  const v4float neg_one(-1.);

  const float _qsp = args->qsp;

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
    _mm_prefetch( (const char *) &p[4].dx, _MM_HINT_T0 );  // added by Jim Kohn
    _mm_prefetch( (const char *) &p[5].dx, _MM_HINT_T0 );  // added by Jim Kohn
    _mm_prefetch( (const char *) &p[6].dx, _MM_HINT_T0 );  // added by Jim Kohn
    _mm_prefetch( (const char *) &p[7].dx, _MM_HINT_T0 );  // added by Jim Kohn

    load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

    // Interpolate fields
    vp0 = (float * ALIGNED(16))(f0 + ii(0));
    vp1 = (float * ALIGNED(16))(f0 + ii(1));
    vp2 = (float * ALIGNED(16))(f0 + ii(2));
    vp3 = (float * ALIGNED(16))(f0 + ii(3));

    load_4x4_tr(vp0,  vp1,  vp2,  vp3,  hax,v0,v1,v2);
    hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );

    load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,hay,v3,v4,v5);
    hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );

    load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,haz,v0,v1,v2);
    haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );

    load_4x4_tr(vp0+12,vp1+12,vp2+12,vp3+12,cbx,v3,cby,v4);
    cbx = fma( v3, dx, cbx );
    cby = fma( v4, dy, cby );

    load_4x2_tr(vp0+16,vp1+16,vp2+16,vp3+16,cbz,v5);
    cbz = fma( v5, dz, cbz );

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
    q  = czero(outbnd,q*qsp);       // Do not accumulate outbnd particles
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
      if( move_p( p0, local_pm, a0, g, _qsp ) ) { /* Unlikely */        \
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
#endif // Method 2.

#if 0
// Method 3.
void
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
  const v4float qsp(args->qsp);
  const v4float one(1.);
  const v4float one_third(1./3.);
  const v4float two_fifteenths(2./15.);
  const v4float neg_one(-1.);

  const float _qsp = args->qsp;

  v4float dx, dy, dz, ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v0, v1, v2, v3, v4, v5;
  v4int   ii, outbnd;

  __m128 t00, t01, t02, t03, t04, t05;

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

  for( ; nq; nq--, p += 4 )
  {
    //--------------------------------------------------------------------------//
    // load_4x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
    // 		    dx, dy, dz, ii );
    //--------------------------------------------------------------------------//

    t00 = _mm_load_ps( (const float *) &p[0].dx );
    t01 = _mm_load_ps( (const float *) &p[1].dx );
    t02 = _mm_load_ps( (const float *) &p[2].dx );
    t03 = _mm_load_ps( (const float *) &p[3].dx );

    t04 = _mm_unpackhi_ps( t00, t01 );
    t00 = _mm_unpacklo_ps( t00, t01 );
    t05 = _mm_unpackhi_ps( t02, t03 );
    t02 = _mm_unpacklo_ps( t02, t03 );

    dx.v = _mm_movelh_ps( t00, t02 );
    dy.v = _mm_movehl_ps( t02, t00 );
    dz.v = _mm_movelh_ps( t04, t05 );
    ii.v = _mm_movehl_ps( t05, t04 );

    //--------------------------------------------------------------------------//
    // Interpolate fields
    //--------------------------------------------------------------------------//

    vp0 = ( float * ALIGNED(16) )( f0 + ii(0) );
    vp1 = ( float * ALIGNED(16) )( f0 + ii(1) );
    vp2 = ( float * ALIGNED(16) )( f0 + ii(2) );
    vp3 = ( float * ALIGNED(16) )( f0 + ii(3) );

    //--------------------------------------------------------------------------//
    // load_4x4_tr( vp0, vp1, vp2, vp3,
    //              hax, v0, v1, v2 );
    //--------------------------------------------------------------------------//

    t00 = _mm_load_ps( (const float *) vp0 );
    t01 = _mm_load_ps( (const float *) vp1 );
    t02 = _mm_load_ps( (const float *) vp2 );
    t03 = _mm_load_ps( (const float *) vp3 );

    t04 = _mm_unpackhi_ps( t00, t01 );
    t00 = _mm_unpacklo_ps( t00, t01 );
    t05 = _mm_unpackhi_ps( t02, t03 );
    t02 = _mm_unpacklo_ps( t02, t03 );

    hax.v = _mm_movelh_ps( t00, t02 );
    v0.v  = _mm_movehl_ps( t02, t00 );
    v1.v  = _mm_movelh_ps( t04, t05 );
    v2.v  = _mm_movehl_ps( t05, t04 );

    //--------------------------------------------------------------------------//
    // hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );
    //--------------------------------------------------------------------------//

    hax.v = _mm_mul_ps( qdt_2mc.v,
			_mm_fmadd_ps( _mm_fmadd_ps( v2.v, dy.v, v1.v ),
				      dz.v,
				      _mm_fmadd_ps( v0.v, dy.v, hax.v ) ) );

    //--------------------------------------------------------------------------//
    // load_4x4_tr( vp0+4, vp1+4, vp2+4, vp3+4,
    //              hay, v3, v4, v5 );
    //--------------------------------------------------------------------------//

    t00 = _mm_load_ps( (const float *) vp0 + 4 );
    t01 = _mm_load_ps( (const float *) vp1 + 4 );
    t02 = _mm_load_ps( (const float *) vp2 + 4 );
    t03 = _mm_load_ps( (const float *) vp3 + 4 );

    t04 = _mm_unpackhi_ps( t00, t01 );
    t00 = _mm_unpacklo_ps( t00, t01 );
    t05 = _mm_unpackhi_ps( t02, t03 );
    t02 = _mm_unpacklo_ps( t02, t03 );

    hay.v = _mm_movelh_ps( t00, t02 );
    v3.v  = _mm_movehl_ps( t02, t00 );
    v4.v  = _mm_movelh_ps( t04, t05 );
    v5.v  = _mm_movehl_ps( t05, t04 );

    //--------------------------------------------------------------------------//
    // hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );
    //--------------------------------------------------------------------------//

    hay.v = _mm_mul_ps( qdt_2mc.v,
			_mm_fmadd_ps( _mm_fmadd_ps( v5.v, dz.v, v4.v ),
				      dx.v,
				      _mm_fmadd_ps( v3.v, dz.v, hay.v ) ) );

    //--------------------------------------------------------------------------//
    // load_4x4_tr( vp0+8, vp1+8, vp2+8, vp3+8,
    //              haz, v0, v1, v2 );
    //--------------------------------------------------------------------------//

    t00 = _mm_load_ps( (const float *) vp0 + 8 );
    t01 = _mm_load_ps( (const float *) vp1 + 8 );
    t02 = _mm_load_ps( (const float *) vp2 + 8 );
    t03 = _mm_load_ps( (const float *) vp3 + 8 );

    t04 = _mm_unpackhi_ps( t00, t01 );
    t00 = _mm_unpacklo_ps( t00, t01 );
    t05 = _mm_unpackhi_ps( t02, t03 );
    t02 = _mm_unpacklo_ps( t02, t03 );

    haz.v = _mm_movelh_ps( t00, t02 );
    v0.v  = _mm_movehl_ps( t02, t00 );
    v1.v  = _mm_movelh_ps( t04, t05 );
    v2.v  = _mm_movehl_ps( t05, t04 );

    //--------------------------------------------------------------------------//
    // haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );
    //--------------------------------------------------------------------------//

    haz.v = _mm_mul_ps( qdt_2mc.v,
			_mm_fmadd_ps( _mm_fmadd_ps( v2.v, dx.v, v1.v ),
				      dy.v,
				      _mm_fmadd_ps( v0.v, dx.v, haz.v ) ) );

    //--------------------------------------------------------------------------//
    // load_4x4_tr( vp0+12, vp1+12, vp2+12, vp3+12,
    //              cbx, v3, cby, v4 );
    //--------------------------------------------------------------------------//

    t00 = _mm_load_ps( (const float *) vp0 + 12 );
    t01 = _mm_load_ps( (const float *) vp1 + 12 );
    t02 = _mm_load_ps( (const float *) vp2 + 12 );
    t03 = _mm_load_ps( (const float *) vp3 + 12 );

    t04 = _mm_unpackhi_ps( t00, t01 );
    t00 = _mm_unpacklo_ps( t00, t01 );
    t05 = _mm_unpackhi_ps( t02, t03 );
    t02 = _mm_unpacklo_ps( t02, t03 );

    cbx.v = _mm_movelh_ps( t00, t02 );
    v3.v  = _mm_movehl_ps( t02, t00 );
    cby.v = _mm_movelh_ps( t04, t05 );
    v4.v  = _mm_movehl_ps( t05, t04 );

    //--------------------------------------------------------------------------//
    // cbx = fma( v3, dx, cbx );
    //--------------------------------------------------------------------------//

    cbx.v = _mm_fmadd_ps( v3.v, dx.v, cbx.v );

    //--------------------------------------------------------------------------//
    // cby = fma( v4, dy, cby );
    //--------------------------------------------------------------------------//

    cby.v = _mm_fmadd_ps( v4.v, dy.v, cby.v );

    //--------------------------------------------------------------------------//
    // load_4x2_tr( vp0+16, vp1+16, vp2+16, vp3+16,
    //              cbz, v5 );
    //--------------------------------------------------------------------------//

    t00 = _mm_setzero_ps();

    t01 = _mm_loadh_pi( _mm_loadl_pi( t00, (__m64 *)vp0+16 ), (__m64 *)vp1+16 );
    t00 = _mm_loadh_pi( _mm_loadl_pi( t00, (__m64 *)vp2+16 ), (__m64 *)vp3+16 );

    cbz.v = _mm_shuffle_ps( t01, t00, 0x88 );
    v5.v  = _mm_shuffle_ps( t01, t00, 0xdd );

    //--------------------------------------------------------------------------//
    // cbz = fma( v5, dz, cbz );
    //--------------------------------------------------------------------------//

    cbz.v = _mm_fmadd_ps( v5.v, dz.v, cbz.v );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    //--------------------------------------------------------------------------//
    // load_4x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
    //              ux, uy, uz, q );
    //--------------------------------------------------------------------------//

    t00 = _mm_load_ps( (const float *) &p[0].ux );
    t01 = _mm_load_ps( (const float *) &p[1].ux );
    t02 = _mm_load_ps( (const float *) &p[2].ux );
    t03 = _mm_load_ps( (const float *) &p[3].ux );

    t04 = _mm_unpackhi_ps( t00, t01 );
    t00 = _mm_unpacklo_ps( t00, t01 );
    t05 = _mm_unpackhi_ps( t02, t03 );
    t02 = _mm_unpacklo_ps( t02, t03 );

    ux.v = _mm_movelh_ps( t00, t02 );
    uy.v = _mm_movehl_ps( t02, t00 );
    uz.v = _mm_movelh_ps( t04, t05 );
    q.v  = _mm_movehl_ps( t05, t04 );

    //--------------------------------------------------------------------------//
    // ux += hax;
    //--------------------------------------------------------------------------//

    ux.v = _mm_add_ps( ux.v, hax.v );

    //--------------------------------------------------------------------------//
    // uy += hay;
    //--------------------------------------------------------------------------//

    uy.v = _mm_add_ps( uy.v, hay.v );

    //--------------------------------------------------------------------------//
    // uz += haz;
    //--------------------------------------------------------------------------//

    uz.v = _mm_add_ps( uz.v, haz.v );

    //--------------------------------------------------------------------------//
    // v0  = qdt_2mc*rsqrt( one + fma( ux, ux, fma( uy, uy, uz*uz ) ) );
    //--------------------------------------------------------------------------//

    t00 = _mm_add_ps( one.v,
		      _mm_fmadd_ps( ux.v,
				    ux.v,
				    _mm_fmadd_ps( uy.v,
						  uy.v,
						  _mm_mul_ps( uz.v,
							      uz.v )
						  )
				    )
		      );

    t01 = _mm_rsqrt_ps( t00 );

    t02 = _mm_fmadd_ps( _mm_set1_ps( 0.5f ),
			_mm_fnmadd_ps( t00,
				       _mm_mul_ps( t01,
						   _mm_mul_ps( t01, t01 ) ),
				       t01 ),
			t01 );

    v0.v = _mm_mul_ps( qdt_2mc.v, t02 );

    //--------------------------------------------------------------------------//
    // v1  = fma( cbx, cbx, fma( cby, cby, cbz*cbz ) );
    //--------------------------------------------------------------------------//

    v1.v = _mm_fmadd_ps( cbx.v,
			 cbx.v,
			 _mm_fmadd_ps( cby.v,
				       cby.v,
				       _mm_mul_ps( cbz.v,
						   cbz.v ) ) );

    //--------------------------------------------------------------------------//
    // v2  = (v0*v0)*v1;
    //--------------------------------------------------------------------------//

    v2.v = _mm_mul_ps( _mm_mul_ps( v0.v, v0.v ), v1.v );

    //--------------------------------------------------------------------------//
    // v3  = v0*fma( fma( two_fifteenths, v2, one_third ), v2, one );
    //--------------------------------------------------------------------------//

    v3.v = _mm_mul_ps( v0.v,
		       _mm_fmadd_ps( _mm_fmadd_ps( two_fifteenths.v,
						   v2.v,
						   one_third.v ),
				     v2.v,
				     one.v ) );

    //--------------------------------------------------------------------------//
    // v4  = v3*rcp( fma( v3*v3, v1, one ) );
    //--------------------------------------------------------------------------//

    t00 = _mm_fmadd_ps( _mm_mul_ps( v3.v, v3.v ), v1.v, one.v );

    t01 = _mm_rcp_ps( t00 );

    v4.v = _mm_mul_ps( v3.v,
		       _mm_fnmadd_ps( t00,
				      _mm_mul_ps( t01, t01 ),
				      _mm_add_ps( t01, t01 ) ) );

    //--------------------------------------------------------------------------//
    // v4 += v4;
    //--------------------------------------------------------------------------//

    v4.v = _mm_add_ps( v4.v, v4.v );

    //--------------------------------------------------------------------------//
    // v0  = fma( fms( uy, cbz, uz*cby ), v3, ux );
    //--------------------------------------------------------------------------//

    v0.v = _mm_fmadd_ps( _mm_fmsub_ps( uy.v,
				       cbz.v,
				       _mm_mul_ps( uz.v, cby.v ) ),
			 v3.v,
			 ux.v );

    //--------------------------------------------------------------------------//
    // v1  = fma( fms( uz, cbx, ux*cbz ), v3, uy );
    //--------------------------------------------------------------------------//

    v1.v = _mm_fmadd_ps( _mm_fmsub_ps( uz.v,
				       cbx.v,
				       _mm_mul_ps( ux.v, cbz.v ) ),
			 v3.v,
			 uy.v );

    //--------------------------------------------------------------------------//
    // v2  = fma( fms( ux, cby, uy*cbx ), v3, uz );
    //--------------------------------------------------------------------------//

    v2.v = _mm_fmadd_ps( _mm_fmsub_ps( ux.v,
				       cby.v,
				       _mm_mul_ps( uy.v, cbx.v ) ),
			 v3.v,
			 uz.v );

    //--------------------------------------------------------------------------//
    // ux  = fma( fms( v1, cbz, v2*cby ), v4, ux );
    //--------------------------------------------------------------------------//

    ux.v = _mm_fmadd_ps( _mm_fmsub_ps( v1.v,
				       cbz.v,
				       _mm_mul_ps( v2.v, cby.v ) ),
			 v4.v,
			 ux.v );

    //--------------------------------------------------------------------------//
    // uy  = fma( fms( v2, cbx, v0*cbz ), v4, uy );
    //--------------------------------------------------------------------------//

    uy.v = _mm_fmadd_ps( _mm_fmsub_ps( v2.v,
				       cbx.v,
				       _mm_mul_ps( v0.v, cbz.v ) ),
			 v4.v,
			 uy.v );

    //--------------------------------------------------------------------------//
    // uz  = fma( fms( v0, cby, v1*cbx ), v4, uz );
    //--------------------------------------------------------------------------//

    uz.v = _mm_fmadd_ps( _mm_fmsub_ps( v0.v,
				       cby.v,
				       _mm_mul_ps( v1.v, cbx.v ) ),
			 v4.v,
			 uz.v );

    //--------------------------------------------------------------------------//
    // ux += hax;
    //--------------------------------------------------------------------------//

    ux.v = _mm_add_ps( ux.v, hax.v );

    //--------------------------------------------------------------------------//
    // uy += hay;
    //--------------------------------------------------------------------------//

    uy.v = _mm_add_ps( uy.v, hay.v );

    //--------------------------------------------------------------------------//
    // uz += haz;
    //--------------------------------------------------------------------------//

    uz.v = _mm_add_ps( uz.v, haz.v );

    //--------------------------------------------------------------------------//
    // store_4x4_tr( ux, uy, uz, q,
    //               &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux );
    //--------------------------------------------------------------------------//

    t04 = _mm_unpackhi_ps( ux.v, uy.v );
    t00 = _mm_unpacklo_ps( ux.v, uy.v );

    t05 = _mm_unpackhi_ps( uz.v, q.v );
    t02 = _mm_unpacklo_ps( uz.v, q.v );

    t01 = _mm_movehl_ps( t02, t00 );
    t00 = _mm_movelh_ps( t00, t02 );
    t02 = _mm_movelh_ps( t04, t05 );
    t03 = _mm_movehl_ps( t05, t04 );

    _mm_store_ps( (float *) &p[0].ux, t00 );
    _mm_store_ps( (float *) &p[1].ux, t01 );
    _mm_store_ps( (float *) &p[2].ux, t02 );
    _mm_store_ps( (float *) &p[3].ux, t03 );

    //--------------------------------------------------------------------------//
    // Update the position of inbnd particles
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v0  = rsqrt( one + fma( ux, ux, fma( uy, uy, uz*uz ) ) );
    //--------------------------------------------------------------------------//

    t00 = _mm_add_ps( one.v,
		      _mm_fmadd_ps( ux.v,
				    ux.v,
				    _mm_fmadd_ps( uy.v,
						  uy.v,
						  _mm_mul_ps( uz.v,
							      uz.v )
						  )
				    )
		      );

    t01 = _mm_rsqrt_ps( t00 );

    v0.v = _mm_fmadd_ps( _mm_set1_ps( 0.5f ),
			 _mm_fnmadd_ps( t00,
					_mm_mul_ps( t01,
						    _mm_mul_ps( t01, t01 ) ),
					t01 ),
			 t01 );

    //--------------------------------------------------------------------------//
    // ux *= cdt_dx;
    //--------------------------------------------------------------------------//

    ux.v = _mm_mul_ps( ux.v, cdt_dx.v );

    //--------------------------------------------------------------------------//
    // uy *= cdt_dy;
    //--------------------------------------------------------------------------//

    uy.v = _mm_mul_ps( uy.v, cdt_dy.v );

    //--------------------------------------------------------------------------//
    // uz *= cdt_dz;
    //--------------------------------------------------------------------------//

    uz.v = _mm_mul_ps( uz.v, cdt_dz.v );

    //--------------------------------------------------------------------------//
    // ux *= v0;
    //--------------------------------------------------------------------------//

    ux.v = _mm_mul_ps( ux.v, v0.v );

    //--------------------------------------------------------------------------//
    // uy *= v0;
    //--------------------------------------------------------------------------//

    uy.v = _mm_mul_ps( uy.v, v0.v );

    //--------------------------------------------------------------------------//
    // uz *= v0;      // ux,uy,uz are normalized displ (relative to cell size)
    //--------------------------------------------------------------------------//

    uz.v = _mm_mul_ps( uz.v, v0.v );

    //--------------------------------------------------------------------------//
    // v0  = dx + ux;
    //--------------------------------------------------------------------------//

    v0.v = _mm_add_ps( dx.v, ux.v );

    //--------------------------------------------------------------------------//
    // v1  = dy + uy;
    //--------------------------------------------------------------------------//

    v1.v = _mm_add_ps( dy.v, uy.v );

    //--------------------------------------------------------------------------//
    // v2  = dz + uz; // New particle midpoint
    //--------------------------------------------------------------------------//

    v2.v = _mm_add_ps( dz.v, uz.v );

    //--------------------------------------------------------------------------//
    // v3  = v0 + ux;
    //--------------------------------------------------------------------------//

    v3.v = _mm_add_ps( v0.v, ux.v );

    //--------------------------------------------------------------------------//
    // v4  = v1 + uy;
    //--------------------------------------------------------------------------//

    v4.v = _mm_add_ps( v1.v, uy.v );

    //--------------------------------------------------------------------------//
    // v5  = v2 + uz; // New particle position
    //--------------------------------------------------------------------------//

    v5.v = _mm_add_ps( v2.v, uz.v );

    //--------------------------------------------------------------------------//
    outbnd = ( v3 > one ) | ( v3 < neg_one ) |
             ( v4 > one ) | ( v4 < neg_one ) |
             ( v5 > one ) | ( v5 < neg_one );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v3  = merge( outbnd, dx, v3 ); // Do not update outbnd particles
    //--------------------------------------------------------------------------//

    v3.v = _mm_or_ps( _mm_andnot_ps( outbnd.v,
				     v3.v ),
		      _mm_and_ps( outbnd.v,
				  dx.v ) );

    //--------------------------------------------------------------------------//
    // v4  = merge( outbnd, dy, v4 );
    //--------------------------------------------------------------------------//

    v4.v = _mm_or_ps( _mm_andnot_ps( outbnd.v,
				     v4.v ),
		      _mm_and_ps( outbnd.v,
				  dy.v ) );

    //--------------------------------------------------------------------------//
    // v5  = merge( outbnd, dz, v5 );
    //--------------------------------------------------------------------------//

    v5.v = _mm_or_ps( _mm_andnot_ps( outbnd.v,
				     v5.v ),
		      _mm_and_ps( outbnd.v,
				  dz.v ) );

    //--------------------------------------------------------------------------//
    // store_4x4_tr( v3, v4, v5, ii,
    // 		     &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx );
    //--------------------------------------------------------------------------//

    t04 = _mm_unpackhi_ps( v3.v, v4.v );
    t00 = _mm_unpacklo_ps( v3.v, v4.v );

    t05 = _mm_unpackhi_ps( v5.v, ii.v );
    t02 = _mm_unpacklo_ps( v5.v, ii.v );

    t01 = _mm_movehl_ps( t02, t00 );
    t00 = _mm_movelh_ps( t00, t02 );
    t02 = _mm_movelh_ps( t04, t05 );
    t03 = _mm_movehl_ps( t05, t04 );

    _mm_store_ps( (float *) &p[0].dx, t00 );
    _mm_store_ps( (float *) &p[1].dx, t01 );
    _mm_store_ps( (float *) &p[2].dx, t02 );
    _mm_store_ps( (float *) &p[3].dx, t03 );

    //--------------------------------------------------------------------------//
    // Accumulate current of inbnd particles
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // q  = czero( outbnd, q*qsp );    // Do not accumulate outbnd particles
    //--------------------------------------------------------------------------//

    q.v = _mm_andnot_ps( outbnd.v,
			 _mm_mul_ps( q.v, qsp.v ) );

    //--------------------------------------------------------------------------//
    // dx = v0;                        // Streak midpoint (valid for inbnd only)
    //--------------------------------------------------------------------------//

    dx.v = v0.v;

    //--------------------------------------------------------------------------//
    // dy = v1;
    //--------------------------------------------------------------------------//

    dy.v = v1.v;

    //--------------------------------------------------------------------------//
    // dz = v2;
    //--------------------------------------------------------------------------//

    dz.v = v2.v;

    //--------------------------------------------------------------------------//
    // v5 = q*ux*uy*uz*one_third;      // Charge conservation correction
    //--------------------------------------------------------------------------//

    v5.v = _mm_mul_ps( q.v,
		       _mm_mul_ps( ux.v,
				   _mm_mul_ps( uy.v,
					       _mm_mul_ps( uz.v,
							   one_third.v )
					       )
				   )
		       );

    //--------------------------------------------------------------------------//
    // Accumulator pointers.
    //--------------------------------------------------------------------------//

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

    //--------------------------------------------------------------------------//
    // ACCUMULATE_J( x, y, z, 0 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v4  = q*ux;
    //--------------------------------------------------------------------------//

    v4.v = _mm_mul_ps( q.v, ux.v );

    //--------------------------------------------------------------------------//
    // v1  = v4*dy;
    //--------------------------------------------------------------------------//

    v1.v = _mm_mul_ps( v4.v, dy.v );

    //--------------------------------------------------------------------------//
    // v0  = v4-v1;
    //--------------------------------------------------------------------------//

    v0.v = _mm_sub_ps( v4.v, v1.v );

    //--------------------------------------------------------------------------//
    // v1 += v4;
    //--------------------------------------------------------------------------//

    v1.v = _mm_add_ps( v1.v, v4.v );

    //--------------------------------------------------------------------------//
    // v4  = one+dz;
    //--------------------------------------------------------------------------//

    v4.v = _mm_add_ps( one.v, dz.v );

    //--------------------------------------------------------------------------//
    // v2  = v0*v4;
    //--------------------------------------------------------------------------//

    v2.v = _mm_mul_ps( v0.v, v4.v );

    //--------------------------------------------------------------------------//
    // v3  = v1*v4;
    //--------------------------------------------------------------------------//

    v3.v = _mm_mul_ps( v1.v, v4.v );

    //--------------------------------------------------------------------------//
    // v4  = one-dz;
    //--------------------------------------------------------------------------//

    v4.v = _mm_sub_ps( one.v, dz.v );

    //--------------------------------------------------------------------------//
    // v0 *= v4;
    //--------------------------------------------------------------------------//

    v0.v = _mm_mul_ps( v0.v, v4.v );

    //--------------------------------------------------------------------------//
    // v1 *= v4;
    //--------------------------------------------------------------------------//

    v1.v = _mm_mul_ps( v1.v, v4.v );

    //--------------------------------------------------------------------------//
    // v0 += v5;
    //--------------------------------------------------------------------------//

    v0.v = _mm_add_ps( v0.v, v5.v );

    //--------------------------------------------------------------------------//
    // v1 -= v5;
    //--------------------------------------------------------------------------//

    v1.v = _mm_sub_ps( v1.v, v5.v );

    //--------------------------------------------------------------------------//
    // v2 -= v5;
    //--------------------------------------------------------------------------//

    v2.v = _mm_sub_ps( v2.v, v5.v );

    //--------------------------------------------------------------------------//
    // v3 += v5;
    //--------------------------------------------------------------------------//

    v3.v = _mm_add_ps( v3.v, v5.v );

    //--------------------------------------------------------------------------//
    // transpose( v0, v1, v2, v3 );
    //--------------------------------------------------------------------------//

    t00  = _mm_unpackhi_ps( v0.v, v1.v );
    v0.v = _mm_unpacklo_ps( v0.v, v1.v );
    t01  = _mm_unpackhi_ps( v2.v, v3.v );
    v2.v = _mm_unpacklo_ps( v2.v, v3.v );

    v1.v = _mm_movehl_ps( v2.v, v0.v );
    v0.v = _mm_movelh_ps( v0.v, v2.v );
    v2.v = _mm_movelh_ps( t00, t01 );
    v3.v = _mm_movehl_ps( t01, t00 );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp0, v0 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp0, _mm_add_ps( _mm_load_ps( vp0 ), v0.v ) );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp1, v1 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp1, _mm_add_ps( _mm_load_ps( vp1 ), v1.v ) );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp2, v2 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp2, _mm_add_ps( _mm_load_ps( vp2 ), v2.v ) );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp3, v3 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp3, _mm_add_ps( _mm_load_ps( vp3 ), v3.v ) );

    //--------------------------------------------------------------------------//
    // ACCUMULATE_J( y, z, x, 4 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v4  = q*uy;
    //--------------------------------------------------------------------------//

    v4.v = _mm_mul_ps( q.v, uy.v );

    //--------------------------------------------------------------------------//
    // v1  = v4*dz;
    //--------------------------------------------------------------------------//

    v1.v = _mm_mul_ps( v4.v, dz.v );

    //--------------------------------------------------------------------------//
    // v0  = v4-v1;
    //--------------------------------------------------------------------------//

    v0.v = _mm_sub_ps( v4.v, v1.v );

    //--------------------------------------------------------------------------//
    // v1 += v4;
    //--------------------------------------------------------------------------//

    v1.v = _mm_add_ps( v1.v, v4.v );

    //--------------------------------------------------------------------------//
    // v4  = one+dx;
    //--------------------------------------------------------------------------//

    v4.v = _mm_add_ps( one.v, dx.v );

    //--------------------------------------------------------------------------//
    // v2  = v0*v4;
    //--------------------------------------------------------------------------//

    v2.v = _mm_mul_ps( v0.v, v4.v );

    //--------------------------------------------------------------------------//
    // v3  = v1*v4;
    //--------------------------------------------------------------------------//

    v3.v = _mm_mul_ps( v1.v, v4.v );

    //--------------------------------------------------------------------------//
    // v4  = one-dx;
    //--------------------------------------------------------------------------//

    v4.v = _mm_sub_ps( one.v, dx.v );

    //--------------------------------------------------------------------------//
    // v0 *= v4;
    //--------------------------------------------------------------------------//

    v0.v = _mm_mul_ps( v0.v, v4.v );

    //--------------------------------------------------------------------------//
    // v1 *= v4;
    //--------------------------------------------------------------------------//

    v1.v = _mm_mul_ps( v1.v, v4.v );

    //--------------------------------------------------------------------------//
    // v0 += v5;
    //--------------------------------------------------------------------------//

    v0.v = _mm_add_ps( v0.v, v5.v );

    //--------------------------------------------------------------------------//
    // v1 -= v5;
    //--------------------------------------------------------------------------//

    v1.v = _mm_sub_ps( v1.v, v5.v );

    //--------------------------------------------------------------------------//
    // v2 -= v5;
    //--------------------------------------------------------------------------//

    v2.v = _mm_sub_ps( v2.v, v5.v );

    //--------------------------------------------------------------------------//
    // v3 += v5;
    //--------------------------------------------------------------------------//

    v3.v = _mm_add_ps( v3.v, v5.v );

    //--------------------------------------------------------------------------//
    // transpose( v0, v1, v2, v3 );
    //--------------------------------------------------------------------------//

    t00  = _mm_unpackhi_ps( v0.v, v1.v );
    v0.v = _mm_unpacklo_ps( v0.v, v1.v );
    t01  = _mm_unpackhi_ps( v2.v, v3.v );
    v2.v = _mm_unpacklo_ps( v2.v, v3.v );

    v1.v = _mm_movehl_ps( v2.v, v0.v );
    v0.v = _mm_movelh_ps( v0.v, v2.v );
    v2.v = _mm_movelh_ps( t00, t01 );
    v3.v = _mm_movehl_ps( t01, t00 );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp0 + 4, v0 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp0 + 4, _mm_add_ps( _mm_load_ps( vp0 + 4 ), v0.v ) );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp1 + 4, v1 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp1 + 4, _mm_add_ps( _mm_load_ps( vp1 + 4 ), v1.v ) );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp2 + 4, v2 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp2 + 4, _mm_add_ps( _mm_load_ps( vp2 + 4 ), v2.v ) );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp3 + 4, v3 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp3 + 4, _mm_add_ps( _mm_load_ps( vp3 + 4 ), v3.v ) );

    //--------------------------------------------------------------------------//
    // ACCUMULATE_J( z, x, y, 8 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v4  = q*uz;
    //--------------------------------------------------------------------------//

    v4.v = _mm_mul_ps( q.v, uz.v );

    //--------------------------------------------------------------------------//
    // v1  = v4*dx;
    //--------------------------------------------------------------------------//

    v1.v = _mm_mul_ps( v4.v, dx.v );

    //--------------------------------------------------------------------------//
    // v0  = v4-v1;
    //--------------------------------------------------------------------------//

    v0.v = _mm_sub_ps( v4.v, v1.v );

    //--------------------------------------------------------------------------//
    // v1 += v4;
    //--------------------------------------------------------------------------//

    v1.v = _mm_add_ps( v1.v, v4.v );

    //--------------------------------------------------------------------------//
    // v4  = one+dy;
    //--------------------------------------------------------------------------//

    v4.v = _mm_add_ps( one.v, dy.v );

    //--------------------------------------------------------------------------//
    // v2  = v0*v4;
    //--------------------------------------------------------------------------//

    v2.v = _mm_mul_ps( v0.v, v4.v );

    //--------------------------------------------------------------------------//
    // v3  = v1*v4;
    //--------------------------------------------------------------------------//

    v3.v = _mm_mul_ps( v1.v, v4.v );

    //--------------------------------------------------------------------------//
    // v4  = one-dy;
    //--------------------------------------------------------------------------//

    v4.v = _mm_sub_ps( one.v, dy.v );

    //--------------------------------------------------------------------------//
    // v0 *= v4;
    //--------------------------------------------------------------------------//

    v0.v = _mm_mul_ps( v0.v, v4.v );

    //--------------------------------------------------------------------------//
    // v1 *= v4;
    //--------------------------------------------------------------------------//

    v1.v = _mm_mul_ps( v1.v, v4.v );

    //--------------------------------------------------------------------------//
    // v0 += v5;
    //--------------------------------------------------------------------------//

    v0.v = _mm_add_ps( v0.v, v5.v );

    //--------------------------------------------------------------------------//
    // v1 -= v5;
    //--------------------------------------------------------------------------//

    v1.v = _mm_sub_ps( v1.v, v5.v );

    //--------------------------------------------------------------------------//
    // v2 -= v5;
    //--------------------------------------------------------------------------//

    v2.v = _mm_sub_ps( v2.v, v5.v );

    //--------------------------------------------------------------------------//
    // v3 += v5;
    //--------------------------------------------------------------------------//

    v3.v = _mm_add_ps( v3.v, v5.v );

    //--------------------------------------------------------------------------//
    // transpose( v0, v1, v2, v3 );
    //--------------------------------------------------------------------------//

    t00  = _mm_unpackhi_ps( v0.v, v1.v );
    v0.v = _mm_unpacklo_ps( v0.v, v1.v );
    t01  = _mm_unpackhi_ps( v2.v, v3.v );
    v2.v = _mm_unpacklo_ps( v2.v, v3.v );

    v1.v = _mm_movehl_ps( v2.v, v0.v );
    v0.v = _mm_movelh_ps( v0.v, v2.v );
    v2.v = _mm_movelh_ps( t00, t01 );
    v3.v = _mm_movehl_ps( t01, t00 );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp0 + 8, v0 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp0 + 8, _mm_add_ps( _mm_load_ps( vp0 + 8 ), v0.v ) );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp1 + 8, v1 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp1 + 8, _mm_add_ps( _mm_load_ps( vp1 + 8 ), v1.v ) );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp2 + 8, v2 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp2 + 8, _mm_add_ps( _mm_load_ps( vp2 + 8 ), v2.v ) );

    //--------------------------------------------------------------------------//
    // increment_4x1( vp3 + 8, v3 );
    //--------------------------------------------------------------------------//

    _mm_store_ps( vp3 + 8, _mm_add_ps( _mm_load_ps( vp3 + 8 ), v3.v ) );

#   undef ACCUMULATE_J

    //--------------------------------------------------------------------------//
    // Update position and accumulate outbnd
    //--------------------------------------------------------------------------//

#   define MOVE_OUTBND(N)                                               \
    if( outbnd(N) ) {                       /* Unlikely */              \
      local_pm->dispx = ux(N);                                          \
      local_pm->dispy = uy(N);                                          \
      local_pm->dispz = uz(N);                                          \
      local_pm->i     = (p - p0) + N;                                   \
      if( move_p( p0, local_pm, a0, g, _qsp ) ) { /* Unlikely */        \
        if( nm<max_nm ) copy_4x1( &pm[nm++], local_pm );                \
        else            itmp++;             /* Unlikely */              \
      }                                                                 \
    }

    //--------------------------------------------------------------------------//
    // MOVE_OUTBND(0);
    //--------------------------------------------------------------------------//

    if( outbnd(0) )                                           /* Unlikely */
    {
      local_pm->dispx = ux(0);
      local_pm->dispy = uy(0);
      local_pm->dispz = uz(0);
      local_pm->i     = (p - p0) + 0;
      if( move_p( p0, local_pm, a0, g, _qsp ) )               /* Unlikely */
      {
        if( nm < max_nm )
	{
	  copy_4x1( &pm[nm++], local_pm );
	}
        else                                                  /* Unlikely */
	{
	  itmp++;
	}
      }
    }

    //--------------------------------------------------------------------------//
    // MOVE_OUTBND(1);
    //--------------------------------------------------------------------------//

    if( outbnd(1) )                                           /* Unlikely */
    {
      local_pm->dispx = ux(1);
      local_pm->dispy = uy(1);
      local_pm->dispz = uz(1);
      local_pm->i     = (p - p0) + 1;
      if( move_p( p0, local_pm, a0, g, _qsp ) )               /* Unlikely */
      {
        if( nm < max_nm )
	{
	  copy_4x1( &pm[nm++], local_pm );
	}
        else                                                  /* Unlikely */
	{
	  itmp++;
	}
      }
    }

    //--------------------------------------------------------------------------//
    // MOVE_OUTBND(2);
    //--------------------------------------------------------------------------//

    if( outbnd(2) )                                           /* Unlikely */
    {
      local_pm->dispx = ux(2);
      local_pm->dispy = uy(2);
      local_pm->dispz = uz(2);
      local_pm->i     = (p - p0) + 2;
      if( move_p( p0, local_pm, a0, g, _qsp ) )               /* Unlikely */
      {
        if( nm < max_nm )
	{
	  copy_4x1( &pm[nm++], local_pm );
	}
        else                                                  /* Unlikely */
	{
	  itmp++;
	}
      }
    }

    //--------------------------------------------------------------------------//
    // MOVE_OUTBND(3);
    //--------------------------------------------------------------------------//

    if( outbnd(3) )                                           /* Unlikely */
    {
      local_pm->dispx = ux(3);
      local_pm->dispy = uy(3);
      local_pm->dispz = uz(3);
      local_pm->i     = (p - p0) + 3;
      if( move_p( p0, local_pm, a0, g, _qsp ) )               /* Unlikely */
      {
        if( nm < max_nm )
	{
	  copy_4x1( &pm[nm++], local_pm );
	}
        else                                                  /* Unlikely */
	{
	  itmp++;
	}
      }
    }

#   undef MOVE_OUTBND

  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
#endif // Method 3.

#if 0
// Method 4.
void
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
  const v4float qsp(args->qsp);
  const v4float one(1.);
  const v4float one_third(1./3.);
  const v4float two_fifteenths(2./15.);
  const v4float neg_one(-1.);

  const float _qsp = args->qsp;

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
    load_4x4_tr_v0(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);
    transpose_v0(dx,dy,dz,ii);

    // Interpolate fields
    vp0 = (float * ALIGNED(16))(f0 + ii(0));
    vp1 = (float * ALIGNED(16))(f0 + ii(1));
    vp2 = (float * ALIGNED(16))(f0 + ii(2));
    vp3 = (float * ALIGNED(16))(f0 + ii(3));

    load_4x4_tr_v1(vp0,  vp1,  vp2,  vp3,  hax,v0,v1,v2);
    transpose_v1(hax,v0,v1,v2);
    hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );

    load_4x4_tr_v2(vp0+4,vp1+4,vp2+4,vp3+4,hay,v3,v4,v5);
    transpose_v2(hay,v3,v4,v5);
    hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );

    load_4x4_tr_v3(vp0+8,vp1+8,vp2+8,vp3+8,haz,v0,v1,v2);
    transpose_v3(haz,v0,v1,v2);
    haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );

    load_4x4_tr_v4(vp0+12,vp1+12,vp2+12,vp3+12,cbx,v3,cby,v4);
    transpose_v4(cbx,v3,cby,v4);
    cbx = fma( v3, dx, cbx );
    cby = fma( v4, dy, cby );

    load_4x2_tr_v0(vp0+16,vp1+16,vp2+16,vp3+16,cbz,v5);
    cbz = fma( v5, dz, cbz );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    load_4x4_tr_v5(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,q);
    transpose_v5(ux,uy,uz,q);
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
    transpose_v6(ux,uy,uz,q);
    store_4x4_tr_v0(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux);
    
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
    transpose_v7(v3,v4,v5,ii);
    store_4x4_tr_v1(v3,v4,v5,ii,&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx);
    
    // Accumulate current of inbnd particles
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step
    q  = czero(outbnd,q*qsp);       // Do not accumulate outbnd particles
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
      if( move_p( p0, local_pm, a0, g, _qsp ) ) { /* Unlikely */        \
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
#endif // Method 4.
