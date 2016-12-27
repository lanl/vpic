using namespace v4;

void
advance_p_pipeline_v4( advance_p_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline )
{
  particle_t           * ALIGNED(128) p0 = args->p0;
  accumulator_t        * ALIGNED(128) a0 = args->a0;
  const interpolator_t * ALIGNED(128) f0 = args->f0;
  const grid_t         *              g  = args->g;

  particle_t           * ALIGNED(128) p;
  particle_mover_t     * ALIGNED(16)  pm;

  float                * ALIGNED(16)  vp00;
  float                * ALIGNED(16)  vp01;
  float                * ALIGNED(16)  vp02;
  float                * ALIGNED(16)  vp03;

  // Basic constants.
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
  v4float v00, v01, v02, v03, v04, v05;
  v4int   ii, outbnd;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 2;

  // Determine which movers are reserved for this pipeline.
  // Movers (16 bytes) should be reserved for pipelines in at least
  // multiples of 8 such that the set of particle movers reserved for
  // a pipeline is 128-byte aligned and a multiple of 128-byte in
  // size.  The host is guaranteed to get enough movers to process its
  // particles with this allocation.

  max_nm = args->max_nm - ( args->np&15 );

  if ( max_nm < 0 ) max_nm = 0;

  DISTRIBUTE( max_nm, 8, pipeline_rank, n_pipeline, itmp, max_nm );

  if ( pipeline_rank == n_pipeline ) max_nm = args->max_nm - itmp;

  pm   = args->pm + itmp;
  nm   = 0;
  itmp = 0;

  // Determine which accumulator array to use.
  // The host gets the first accumulator array.

  a0 += ( 1 + pipeline_rank ) *
        POW2_CEIL( (args->nx+2)*(args->ny+2)*(args->nz+2), 2 );

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=4 )
  {
    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------
    load_4x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		 dx, dy, dz, ii );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( float * ALIGNED(16) ) ( f0 + ii( 0) );
    vp01 = ( float * ALIGNED(16) ) ( f0 + ii( 1) );
    vp02 = ( float * ALIGNED(16) ) ( f0 + ii( 2) );
    vp03 = ( float * ALIGNED(16) ) ( f0 + ii( 3) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00, vp01, vp02, vp03,
		 hax, v00, v01, v02 );

    //--------------------------------------------------------------------------
    // Experiment to understand cost of transposing data.
    //--------------------------------------------------------------------------
    // #if 0
    transpose( hax, v00, v01, v02 );
    transpose( hax, v00, v01, v02 );
    // #endif

    hax = qdt_2mc*fma( fma( v02, dy, v01 ), dz, fma( v00, dy, hax ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00+4, vp01+4, vp02+4, vp03+4,
		 hay, v03, v04, v05 );

    hay = qdt_2mc*fma( fma( v05, dz, v04 ), dx, fma( v03, dz, hay ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00+8, vp01+8, vp02+8, vp03+8,
		 haz, v00, v01, v02 );

    haz = qdt_2mc*fma( fma( v02, dx, v01 ), dy, fma( v00, dx, haz ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00+12, vp01+12, vp02+12, vp03+12,
		 cbx, v03, cby, v04 );

    cbx = fma( v03, dx, cbx );

    cby = fma( v04, dy, cby );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles, final.
    //--------------------------------------------------------------------------
    load_4x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
		 cbz, v05 );

    cbz = fma( v05, dz, cbz );

    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------
    load_4x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		 ux, uy, uz, q );

    //--------------------------------------------------------------------------
    // Update momentum.
    //--------------------------------------------------------------------------
    // For a 5-10% performance hit, v00 = qdt_2mc/sqrt(blah) is a few ulps more
    // accurate (but still quite in the noise numerically) for cyclotron
    // frequencies approaching the nyquist frequency.
    //--------------------------------------------------------------------------

    ux  += hax;
    uy  += hay;
    uz  += haz;

    v00  = qdt_2mc*rsqrt( one + fma( ux, ux, fma( uy, uy, uz*uz ) ) );
    v01  = fma( cbx, cbx, fma( cby, cby, cbz*cbz ) );
    v02  = (v00*v00)*v01;
    v03  = v00*fma( fma( two_fifteenths, v02, one_third ), v02, one );
    v04  = v03*rcp( fma( v03*v03, v01, one ) );
    v04 += v04;

    v00  = fma( fms(  uy, cbz,  uz*cby ), v03, ux );
    v01  = fma( fms(  uz, cbx,  ux*cbz ), v03, uy );
    v02  = fma( fms(  ux, cby,  uy*cbx ), v03, uz );

    ux   = fma( fms( v01, cbz, v02*cby ), v04, ux );
    uy   = fma( fms( v02, cbx, v00*cbz ), v04, uy );
    uz   = fma( fms( v00, cby, v01*cbx ), v04, uz );

    ux  += hax;
    uy  += hay;
    uz  += haz;

    //--------------------------------------------------------------------------
    // Store particle data.
    //--------------------------------------------------------------------------
    store_4x4_tr( ux, uy, uz, q,
		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux );

    //--------------------------------------------------------------------------
    // Update the position of in bound particles.
    //--------------------------------------------------------------------------
    v00 = rsqrt( one + fma( ux, ux, fma( uy, uy, uz*uz ) ) );

    ux *= cdt_dx;
    uy *= cdt_dy;
    uz *= cdt_dz;

    ux *= v00;
    uy *= v00;
    uz *= v00;      // ux,uy,uz are normalized displ (relative to cell size)

    v00 =  dx + ux;
    v01 =  dy + uy;
    v02 =  dz + uz; // New particle midpoint

    v03 = v00 + ux;
    v04 = v01 + uy;
    v05 = v02 + uz; // New particle position

    //--------------------------------------------------------------------------
    // Determine which particles are out of bounds.
    //--------------------------------------------------------------------------
    outbnd = ( v03 > one ) | ( v03 < neg_one ) |
             ( v04 > one ) | ( v04 < neg_one ) |
             ( v05 > one ) | ( v05 < neg_one );

    v03 = merge( outbnd, dx, v03 ); // Do not update outbnd particles
    v04 = merge( outbnd, dy, v04 );
    v05 = merge( outbnd, dz, v05 );

    //--------------------------------------------------------------------------
    // Store particle data, final.
    //--------------------------------------------------------------------------
    store_4x4_tr( v03, v04, v05, ii,
		  &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx );

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q  = czero( outbnd, q*qsp );   // Do not accumulate outbnd particles

    dx = v00;                      // Streak midpoint (valid for inbnd only)
    dy = v01;
    dz = v02;

    v05 = q*ux*uy*uz*one_third;    // Charge conservation correction

    //--------------------------------------------------------------------------
    // Set current density accumulation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( float * ALIGNED(16) ) ( a0 + ii( 0) );
    vp01 = ( float * ALIGNED(16) ) ( a0 + ii( 1) );
    vp02 = ( float * ALIGNED(16) ) ( a0 + ii( 2) );
    vp03 = ( float * ALIGNED(16) ) ( a0 + ii( 3) );

    //--------------------------------------------------------------------------
    // Accumulate current density.
    //--------------------------------------------------------------------------
#   define ACCUMULATE_J(X,Y,Z,offset)                                  \
    v04  = q*u##X;    /* v04 = q ux                            */      \
    v01  = v04*d##Y;  /* v01 = q ux dy                         */      \
    v00  = v04-v01;   /* v00 = q ux (1-dy)                     */      \
    v01 += v04;       /* v01 = q ux (1+dy)                     */      \
    v04  = one+d##Z;  /* v04 = 1+dz                            */      \
    v02  = v00*v04;   /* v02 = q ux (1-dy)(1+dz)               */      \
    v03  = v01*v04;   /* v03 = q ux (1+dy)(1+dz)               */      \
    v04  = one-d##Z;  /* v04 = 1-dz                            */      \
    v00 *= v04;       /* v00 = q ux (1-dy)(1-dz)               */      \
    v01 *= v04;       /* v01 = q ux (1+dy)(1-dz)               */      \
    v00 += v05;       /* v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v01 -= v05;       /* v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v02 -= v05;       /* v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v03 += v05;       /* v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */      \
    transpose( v00, v01, v02, v03 );                                   \
    increment_4x1( vp00+offset, v00 );                                 \
    increment_4x1( vp01+offset, v01 );                                 \
    increment_4x1( vp02+offset, v02 );                                 \
    increment_4x1( vp03+offset, v03 )

    ACCUMULATE_J( x, y, z, 0 );
    ACCUMULATE_J( y, z, x, 4 );
    ACCUMULATE_J( z, x, y, 8 );

#   undef ACCUMULATE_J

    //--------------------------------------------------------------------------
    // Update position and accumulate current density for out of bounds
    // particles.
    //--------------------------------------------------------------------------

#   define MOVE_OUTBND(N)                                               \
    if ( outbnd(N) )                                /* Unlikely */      \
    {                                                                   \
      local_pm->dispx = ux(N);                                          \
      local_pm->dispy = uy(N);                                          \
      local_pm->dispz = uz(N);                                          \
      local_pm->i     = ( p - p0 ) + N;                                 \
      if ( move_p( p0, local_pm, a0, g, _qsp ) )    /* Unlikely */      \
      {                                                                 \
        if ( nm < max_nm )                                              \
        {                                                               \
	  copy_4x1( &pm[nm++], local_pm );                              \
	}                                                               \
        else                                        /* Unlikely */      \
	{                                                               \
	  itmp++;                                                       \
	}                                                               \
      }                                                                 \
    }

    MOVE_OUTBND( 0);
    MOVE_OUTBND( 1);
    MOVE_OUTBND( 2);
    MOVE_OUTBND( 3);

#   undef MOVE_OUTBND
  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
