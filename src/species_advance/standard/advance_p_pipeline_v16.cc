using namespace v16;

//----------------------------------------------------------------------------//
// Method 4
//----------------------------------------------------------------------------//
// This method processes 16 particles at a time instead of 32.
//----------------------------------------------------------------------------//
// This method processes the particles in the same order as the reference
// implementation and gives good reproducibility. This is achieved using
// modified load_16x8_tr_p and store_16x8_tr_p functions which load or store
// the particle data in the correct order in a single step instead of using
// two steps.
//----------------------------------------------------------------------------//

void
advance_p_pipeline_v16( advance_p_pipeline_args_t * args,
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
  float                * ALIGNED(16)  vp04;
  float                * ALIGNED(16)  vp05;
  float                * ALIGNED(16)  vp06;
  float                * ALIGNED(16)  vp07;
  float                * ALIGNED(16)  vp08;
  float                * ALIGNED(16)  vp09;
  float                * ALIGNED(16)  vp10;
  float                * ALIGNED(16)  vp11;
  float                * ALIGNED(16)  vp12;
  float                * ALIGNED(16)  vp13;
  float                * ALIGNED(16)  vp14;
  float                * ALIGNED(16)  vp15;

  // Basic constants.
  const v16float qdt_2mc(args->qdt_2mc);
  const v16float cdt_dx(args->cdt_dx);
  const v16float cdt_dy(args->cdt_dy);
  const v16float cdt_dz(args->cdt_dz);
  const v16float qsp(args->qsp);
  const v16float one(1.0f);
  const v16float one_third(1.0f/3.0f);
  const v16float two_fifteenths(2.0f/15.0f);
  const v16float neg_one(-1.0f);

  const float _qsp = args->qsp;

  v16float dx, dy, dz, ux, uy, uz, q;
  v16float hax, hay, haz, cbx, cby, cbz;
  v16float v00, v01, v02, v03, v04, v05, v06, v07;
  v16float v08, v09, v10, v11, v12, v13, v14, v15;
  v16int   ii, outbnd;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 4;

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

  int num_load_16x16_tr = 0;
  int num_load_16x16_bc = 0;
  for( ; nq; nq--, p+=16 )
  {
    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------
    load_16x8_tr_p( &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
                    &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		    dx, dy, dz, ii, ux, uy, uz, q );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( float * ALIGNED(16) ) ( f0 + ii( 0) );
    vp01 = ( float * ALIGNED(16) ) ( f0 + ii( 1) );
    vp02 = ( float * ALIGNED(16) ) ( f0 + ii( 2) );
    vp03 = ( float * ALIGNED(16) ) ( f0 + ii( 3) );
    vp04 = ( float * ALIGNED(16) ) ( f0 + ii( 4) );
    vp05 = ( float * ALIGNED(16) ) ( f0 + ii( 5) );
    vp06 = ( float * ALIGNED(16) ) ( f0 + ii( 6) );
    vp07 = ( float * ALIGNED(16) ) ( f0 + ii( 7) );
    vp08 = ( float * ALIGNED(16) ) ( f0 + ii( 8) );
    vp09 = ( float * ALIGNED(16) ) ( f0 + ii( 9) );
    vp10 = ( float * ALIGNED(16) ) ( f0 + ii(10) );
    vp11 = ( float * ALIGNED(16) ) ( f0 + ii(11) );
    vp12 = ( float * ALIGNED(16) ) ( f0 + ii(12) );
    vp13 = ( float * ALIGNED(16) ) ( f0 + ii(13) );
    vp14 = ( float * ALIGNED(16) ) ( f0 + ii(14) );
    vp15 = ( float * ALIGNED(16) ) ( f0 + ii(15) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00, vp01, vp02, vp03,
                   vp04, vp05, vp06, vp07,
                   vp08, vp09, vp10, vp11,
                   vp12, vp13, vp14, vp15,
                   hax, v00, v01, v02, hay, v03, v04, v05,
                   haz, v06, v07, v08, cbx, v09, cby, v10 );

    //--------------------------------------------------------------------------
    // Experiment to explore broadcast of data elements to vectors when all 16
    // particles are in the same cell.
    //--------------------------------------------------------------------------
    // Perhaps a better way would be to broadcast ii(0) into all the elements
    // of an ii_tmp simd vector, do a simd compare and then a simd reduction
    // to a single int value.  But for now, use the brute force approach.
    //--------------------------------------------------------------------------
    #if 0
    if ( ii( 1) == ii(0) && ii( 2) == ii(0) && ii( 3) == ii(0) &&
	 ii( 4) == ii(0) && ii( 5) == ii(0) && ii( 6) == ii(0) &&
	 ii( 7) == ii(0) && ii( 8) == ii(0) && ii( 9) == ii(0) &&
	 ii(10) == ii(0) && ii(11) == ii(0) && ii(12) == ii(0) &&
	 ii(13) == ii(0) && ii(14) == ii(0) && ii(15) == ii(0) )
    {
      num_load_16x16_bc++;
      load_16x16_bc( vp00,
                     hax, v00, v01, v02, hay, v03, v04, v05,
                     haz, v06, v07, v08, cbx, v09, cby, v10 );
    }
    else
    {
      num_load_16x16_tr++;
      load_16x16_tr( vp00, vp01, vp02, vp03,
                     vp04, vp05, vp06, vp07,
                     vp08, vp09, vp10, vp11,
                     vp12, vp13, vp14, vp15,
                     hax, v00, v01, v02, hay, v03, v04, v05,
                     haz, v06, v07, v08, cbx, v09, cby, v10 );
    }
    #endif

    //--------------------------------------------------------------------------
    // Experiment to understand cost of transposing data.
    //--------------------------------------------------------------------------
    #if 0
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    transpose( hax, v00, v01, v02, hay, v03, v04, v05,
	       haz, v06, v07, v08, cbx, v09, cby, v10 );
    #endif

    hax = qdt_2mc*fma( fma( v02, dy, v01 ), dz, fma( v00, dy, hax ) );

    hay = qdt_2mc*fma( fma( v05, dz, v04 ), dx, fma( v03, dz, hay ) );

    haz = qdt_2mc*fma( fma( v08, dx, v07 ), dy, fma( v06, dx, haz ) );

    cbx = fma( v09, dx, cbx );

    cby = fma( v10, dy, cby );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles, final.
    //--------------------------------------------------------------------------
    // load_16x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
    //               vp04+16, vp05+16, vp06+16, vp07+16,
    //               vp08+16, vp09+16, vp10+16, vp11+16,
    //               vp12+16, vp13+16, vp14+16, vp15+16,
    //               cbz, v05 );

    load_16x16_tr( vp00+16, vp01+16, vp02+16, vp03+16,
                   vp04+16, vp05+16, vp06+16, vp07+16,
                   vp08+16, vp09+16, vp10+16, vp11+16,
                   vp12+16, vp13+16, vp14+16, vp15+16,
                   cbz, v05, v00, v01, v02, v03, v04, v06,
                   v07, v08, v09, v10, v11, v12, v13, v14 );

    //--------------------------------------------------------------------------
    // Experiment to explore broadcast of data elements to vectors when all 16
    // particles are in the same cell.
    //--------------------------------------------------------------------------
    // Perhaps a better way would be to broadcast ii(0) into all the elements
    // of an ii_tmp simd vector, do a simd compare and then a simd reduction
    // to a single int value.  But for now, use the brute force approach.
    //--------------------------------------------------------------------------
    #if 0
    if ( ii( 1) == ii(0) && ii( 2) == ii(0) && ii( 3) == ii(0) &&
	 ii( 4) == ii(0) && ii( 5) == ii(0) && ii( 6) == ii(0) &&
	 ii( 7) == ii(0) && ii( 8) == ii(0) && ii( 9) == ii(0) &&
	 ii(10) == ii(0) && ii(11) == ii(0) && ii(12) == ii(0) &&
	 ii(13) == ii(0) && ii(14) == ii(0) && ii(15) == ii(0) )
    {
      load_16x2_bc( vp00+16,
                    cbz, v05 );
    }
    else
    {
      load_16x16_tr( vp00+16, vp01+16, vp02+16, vp03+16,
                     vp04+16, vp05+16, vp06+16, vp07+16,
                     vp08+16, vp09+16, vp10+16, vp11+16,
                     vp12+16, vp13+16, vp14+16, vp15+16,
                     cbz, v05, v00, v01, v02, v03, v04, v06,
                     v07, v08, v09, v10, v11, v12, v13, v14 );
    }
    #endif

    cbz = fma( v05, dz, cbz );

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

    // Store ux, uy, uz in v06, v07, v08 so particle velocity store can be done
    // later with the particle positions.
    v06  = ux;
    v07  = uy;
    v08  = uz;

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
    store_16x8_tr_p( v03, v04, v05, ii, v06, v07, v08, q,
                     &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
                     &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx );

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q  = czero( outbnd, q*qsp );   // Do not accumulate outbnd particles

    dx = v00;                      // Streak midpoint (valid for inbnd only)
    dy = v01;
    dz = v02;

    v13 = q*ux*uy*uz*one_third;    // Charge conservation correction

    //--------------------------------------------------------------------------
    // Set current density accumulation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( float * ALIGNED(16) ) ( a0 + ii( 0) );
    vp01 = ( float * ALIGNED(16) ) ( a0 + ii( 1) );
    vp02 = ( float * ALIGNED(16) ) ( a0 + ii( 2) );
    vp03 = ( float * ALIGNED(16) ) ( a0 + ii( 3) );
    vp04 = ( float * ALIGNED(16) ) ( a0 + ii( 4) );
    vp05 = ( float * ALIGNED(16) ) ( a0 + ii( 5) );
    vp06 = ( float * ALIGNED(16) ) ( a0 + ii( 6) );
    vp07 = ( float * ALIGNED(16) ) ( a0 + ii( 7) );
    vp08 = ( float * ALIGNED(16) ) ( a0 + ii( 8) );
    vp09 = ( float * ALIGNED(16) ) ( a0 + ii( 9) );
    vp10 = ( float * ALIGNED(16) ) ( a0 + ii(10) );
    vp11 = ( float * ALIGNED(16) ) ( a0 + ii(11) );
    vp12 = ( float * ALIGNED(16) ) ( a0 + ii(12) );
    vp13 = ( float * ALIGNED(16) ) ( a0 + ii(13) );
    vp14 = ( float * ALIGNED(16) ) ( a0 + ii(14) );
    vp15 = ( float * ALIGNED(16) ) ( a0 + ii(15) );

    //--------------------------------------------------------------------------
    // Accumulate current density.
    //--------------------------------------------------------------------------
    // Accumulate Jx for 16 particles into the v0-v3 vectors.
    v12  = q*ux;     // v12 = q ux
    v01  = v12*dy;   // v01 = q ux dy
    v00  = v12-v01;  // v00 = q ux (1-dy)
    v01 += v12;      // v01 = q ux (1+dy)
    v12  = one+dz;   // v12 = 1+dz
    v02  = v00*v12;  // v02 = q ux (1-dy)(1+dz)
    v03  = v01*v12;  // v03 = q ux (1+dy)(1+dz)
    v12  = one-dz;   // v12 = 1-dz
    v00 *= v12;      // v00 = q ux (1-dy)(1-dz)
    v01 *= v12;      // v01 = q ux (1+dy)(1-dz)
    v00 += v13;      // v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ]
    v01 -= v13;      // v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ]
    v02 -= v13;      // v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ]
    v03 += v13;      // v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ]

    // Accumulate Jy for 16 particles into the v4-v7 vectors.
    v12  = q*uy;     // v12 = q uy
    v05  = v12*dz;   // v05 = q uy dz
    v04  = v12-v05;  // v04 = q uy (1-dz)
    v05 += v12;      // v05 = q uy (1+dz)
    v12  = one+dx;   // v12 = 1+dx
    v06  = v04*v12;  // v06 = q uy (1-dz)(1+dx)
    v07  = v05*v12;  // v07 = q uy (1+dz)(1+dx)
    v12  = one-dx;   // v12 = 1-dx
    v04 *= v12;      // v04 = q uy (1-dz)(1-dx)
    v05 *= v12;      // v05 = q uy (1+dz)(1-dx)
    v04 += v13;      // v04 = q uy [ (1-dz)(1-dx) + ux*uz/3 ]
    v05 -= v13;      // v05 = q uy [ (1+dz)(1-dx) - ux*uz/3 ]
    v06 -= v13;      // v06 = q uy [ (1-dz)(1+dx) - ux*uz/3 ]
    v07 += v13;      // v07 = q uy [ (1+dz)(1+dx) + ux*uz/3 ]

    // Accumulate Jz for 16 particles into the v8-v11 vectors.
    v12  = q*uz;     // v12 = q uz
    v09  = v12*dx;   // v09 = q uz dx
    v08  = v12-v09;  // v08 = q uz (1-dx)
    v09 += v12;      // v09 = q uz (1+dx)
    v12  = one+dy;   // v12 = 1+dy
    v10  = v08*v12;  // v10 = q uz (1-dx)(1+dy)
    v11  = v09*v12;  // v11 = q uz (1+dx)(1+dy)
    v12  = one-dy;   // v12 = 1-dy
    v08 *= v12;      // v08 = q uz (1-dx)(1-dy)
    v09 *= v12;      // v09 = q uz (1+dx)(1-dy)
    v08 += v13;      // v08 = q uz [ (1-dx)(1-dy) + ux*uy/3 ]
    v09 -= v13;      // v09 = q uz [ (1+dx)(1-dy) - ux*uy/3 ]
    v10 -= v13;      // v10 = q uz [ (1-dx)(1+dy) - ux*uy/3 ]
    v11 += v13;      // v11 = q uz [ (1+dx)(1+dy) + ux*uy/3 ]

    // Zero the v12-v15 vectors prior to transposing the data.
    v12 = 0.0;
    v13 = 0.0;
    v14 = 0.0;
    v15 = 0.0;

    // Transpose the data in vectors v0-v15 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v00, v01, v02, v03, v04, v05, v06, v07,
	       v08, v09, v10, v11, v12, v13, v14, v15 );

    // Add the contributions to Jx, Jy and Jz from 16 particles into the
    // accumulator arrays for Jx, Jy and Jz.
    increment_16x1( vp00, v00 );
    increment_16x1( vp01, v01 );
    increment_16x1( vp02, v02 );
    increment_16x1( vp03, v03 );
    increment_16x1( vp04, v04 );
    increment_16x1( vp05, v05 );
    increment_16x1( vp06, v06 );
    increment_16x1( vp07, v07 );
    increment_16x1( vp08, v08 );
    increment_16x1( vp09, v09 );
    increment_16x1( vp10, v10 );
    increment_16x1( vp11, v11 );
    increment_16x1( vp12, v12 );
    increment_16x1( vp13, v13 );
    increment_16x1( vp14, v14 );
    increment_16x1( vp15, v15 );

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
	  v4::copy_4x1( &pm[nm++], local_pm ); 	                        \
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
    MOVE_OUTBND( 4);
    MOVE_OUTBND( 5);
    MOVE_OUTBND( 6);
    MOVE_OUTBND( 7);
    MOVE_OUTBND( 8);
    MOVE_OUTBND( 9);
    MOVE_OUTBND(10);
    MOVE_OUTBND(11);
    MOVE_OUTBND(12);
    MOVE_OUTBND(13);
    MOVE_OUTBND(14);
    MOVE_OUTBND(15);

#   undef MOVE_OUTBND
  }

  if ( world_rank == 17 && pipeline_rank == 0 )
  {
    std::cout << "world_rank: " << world_rank
	      << "   pipeline_rank: " << pipeline_rank
	      << "   num_load_16x16_tr: " << num_load_16x16_tr
	      << "   num_load_16x16_bc: " << num_load_16x16_bc
	      << std::flush << std::endl;
  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
