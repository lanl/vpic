using namespace v16;

#if 0
//----------------------------------------------------------------------------//
// Method 1
//----------------------------------------------------------------------------//
// This method processes the particles in the same order as the reference
// implementation and gives good reproducibility. This is achieved using
// modified load_16x16_tr_p and store_16x16_tr_p functions which load or
// store the particle data in the correct order in a single step instead
// of using two steps.
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

  float                * ALIGNED(16)  vp00_a;
  float                * ALIGNED(16)  vp01_a;
  float                * ALIGNED(16)  vp02_a;
  float                * ALIGNED(16)  vp03_a;
  float                * ALIGNED(16)  vp04_a;
  float                * ALIGNED(16)  vp05_a;
  float                * ALIGNED(16)  vp06_a;
  float                * ALIGNED(16)  vp07_a;
  float                * ALIGNED(16)  vp08_a;
  float                * ALIGNED(16)  vp09_a;
  float                * ALIGNED(16)  vp10_a;
  float                * ALIGNED(16)  vp11_a;
  float                * ALIGNED(16)  vp12_a;
  float                * ALIGNED(16)  vp13_a;
  float                * ALIGNED(16)  vp14_a;
  float                * ALIGNED(16)  vp15_a;

  float                * ALIGNED(16)  vp00_b;
  float                * ALIGNED(16)  vp01_b;
  float                * ALIGNED(16)  vp02_b;
  float                * ALIGNED(16)  vp03_b;
  float                * ALIGNED(16)  vp04_b;
  float                * ALIGNED(16)  vp05_b;
  float                * ALIGNED(16)  vp06_b;
  float                * ALIGNED(16)  vp07_b;
  float                * ALIGNED(16)  vp08_b;
  float                * ALIGNED(16)  vp09_b;
  float                * ALIGNED(16)  vp10_b;
  float                * ALIGNED(16)  vp11_b;
  float                * ALIGNED(16)  vp12_b;
  float                * ALIGNED(16)  vp13_b;
  float                * ALIGNED(16)  vp14_b;
  float                * ALIGNED(16)  vp15_b;

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

  v16float dx_a, dy_a, dz_a, ux_a, uy_a, uz_a, q_a;
  v16float hax_a, hay_a, haz_a, cbx_a, cby_a, cbz_a;
  v16float v00_a, v01_a, v02_a, v03_a, v04_a, v05_a, v06_a, v07_a;
  v16float v08_a, v09_a, v10_a, v11_a, v12_a, v13_a, v14_a, v15_a;
  v16int   ii_a, outbnd_a;

  v16float dx_b, dy_b, dz_b, ux_b, uy_b, uz_b, q_b;
  v16float hax_b, hay_b, haz_b, cbx_b, cby_b, cbz_b;
  v16float v00_b, v01_b, v02_b, v03_b, v04_b, v05_b, v06_b, v07_b;
  v16float v08_b, v09_b, v10_b, v11_b, v12_b, v13_b, v14_b, v15_b;
  v16int   ii_b, outbnd_b;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 32, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 5;

  // Determine which movers are reserved for this pipeline.
  // Movers (16 bytes) should be reserved for pipelines in at least
  // multiples of 8 such that the set of particle movers reserved for
  // a pipeline is 128-byte aligned and a multiple of 128-byte in
  // size.  The host is guaranteed to get enough movers to process its
  // particles with this allocation.

  max_nm = args->max_nm - ( args->np&15 );

  if ( max_nm < 0 ) max_nm = 0;

  DISTRIBUTE( max_nm, 8, pipeline_rank, n_pipeline, itmp, max_nm );

  if ( pipeline_rank==n_pipeline ) max_nm = args->max_nm - itmp;

  pm   = args->pm + itmp;
  nm   = 0;
  itmp = 0;

  // Determine which accumulator array to use.
  // The host gets the first accumulator array.

  a0 += ( 1 + pipeline_rank ) *
        POW2_CEIL( (args->nx+2)*(args->ny+2)*(args->nz+2), 2 );

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=32 )
  {
    load_16x16_tr_p( &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
		     &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		     &p[16].dx, &p[18].dx, &p[20].dx, &p[22].dx,
		     &p[24].dx, &p[26].dx, &p[28].dx, &p[30].dx,
		     dx_a, dy_a, dz_a, ii_a, ux_a, uy_a, uz_a, q_a,
		     dx_b, dy_b, dz_b, ii_b, ux_b, uy_b, uz_b, q_b );

    //--------------------------------------------------------------------------
    // Interpolate fields for particles of block a.
    //--------------------------------------------------------------------------
    vp00_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 0) );
    vp01_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 1) );
    vp02_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 2) );
    vp03_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 3) );
    vp04_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 4) );
    vp05_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 5) );
    vp06_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 6) );
    vp07_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 7) );
    vp08_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 8) );
    vp09_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 9) );
    vp10_a = ( float * ALIGNED(16) ) ( f0 + ii_a(10) );
    vp11_a = ( float * ALIGNED(16) ) ( f0 + ii_a(11) );
    vp12_a = ( float * ALIGNED(16) ) ( f0 + ii_a(12) );
    vp13_a = ( float * ALIGNED(16) ) ( f0 + ii_a(13) );
    vp14_a = ( float * ALIGNED(16) ) ( f0 + ii_a(14) );
    vp15_a = ( float * ALIGNED(16) ) ( f0 + ii_a(15) );

    //--------------------------------------------------------------------------
    // Load data for 16 particles for block a.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00_a, vp01_a, vp02_a, vp03_a,
		   vp04_a, vp05_a, vp06_a, vp07_a,
		   vp08_a, vp09_a, vp10_a, vp11_a,
		   vp12_a, vp13_a, vp14_a, vp15_a,
		   hax_a, v00_a, v01_a, v02_a, hay_a, v03_a, v04_a, v05_a,
		   haz_a, v06_a, v07_a, v08_a, cbx_a, v09_a, cby_a, v10_a );

    //--------------------------------------------------------------------------
    // Process particle block a.
    //--------------------------------------------------------------------------
    hax_a = qdt_2mc*
            fma( fma( v02_a, dy_a, v01_a ), dz_a, fma( v00_a, dy_a, hax_a ) );

    hay_a = qdt_2mc*
            fma( fma( v05_a, dz_a, v04_a ), dx_a, fma( v03_a, dz_a, hay_a ) );

    haz_a = qdt_2mc*
            fma( fma( v08_a, dx_a, v07_a ), dy_a, fma( v06_a, dx_a, haz_a ) );

    cbx_a = fma( v09_a, dx_a, cbx_a );

    cby_a = fma( v10_a, dy_a, cby_a );

    load_16x2_tr( vp00_a+16, vp01_a+16, vp02_a+16, vp03_a+16,
		  vp04_a+16, vp05_a+16, vp06_a+16, vp07_a+16,
		  vp08_a+16, vp09_a+16, vp10_a+16, vp11_a+16,
		  vp12_a+16, vp13_a+16, vp14_a+16, vp15_a+16,
		  cbz_a, v05_a );

    cbz_a = fma( v05_a, dz_a, cbz_a );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    ux_a += hax_a;
    uy_a += hay_a;
    uz_a += haz_a;
    v00_a  = qdt_2mc*rsqrt( one + fma( ux_a, ux_a,
                                       fma( uy_a, uy_a, uz_a*uz_a ) ) );
    v01_a  = fma( cbx_a, cbx_a, fma( cby_a, cby_a, cbz_a*cbz_a ) );
    v02_a  = (v00_a*v00_a)*v01_a;
    v03_a  = v00_a*fma( fma( two_fifteenths, v02_a, one_third ), v02_a, one );
    v04_a  = v03_a*rcp( fma( v03_a*v03_a, v01_a, one ) );
    v04_a += v04_a;
    v00_a  = fma( fms(  uy_a, cbz_a,  uz_a*cby_a ), v03_a, ux_a );
    v01_a  = fma( fms(  uz_a, cbx_a,  ux_a*cbz_a ), v03_a, uy_a );
    v02_a  = fma( fms(  ux_a, cby_a,  uy_a*cbx_a ), v03_a, uz_a );
    ux_a   = fma( fms( v01_a, cbz_a, v02_a*cby_a ), v04_a, ux_a );
    uy_a   = fma( fms( v02_a, cbx_a, v00_a*cbz_a ), v04_a, uy_a );
    uz_a   = fma( fms( v00_a, cby_a, v01_a*cbx_a ), v04_a, uz_a );
    ux_a  += hax_a;
    uy_a  += hay_a;
    uz_a  += haz_a;

    // Store ux, uy, uz in v6, v7, v8 so store_16x16_tr can be used below.
    v06_a  = ux_a;
    v07_a  = uy_a;
    v08_a  = uz_a;

    // Update the position of inbnd particles
    v00_a  = rsqrt( one + fma( ux_a, ux_a, fma( uy_a, uy_a, uz_a*uz_a ) ) );
    ux_a  *= cdt_dx;
    uy_a  *= cdt_dy;
    uz_a  *= cdt_dz;
    ux_a  *= v00_a;
    uy_a  *= v00_a;
    uz_a  *= v00_a;        // ux,uy,uz are normalized displ (relative to cell size)
    v00_a  =  dx_a + ux_a;
    v01_a  =  dy_a + uy_a;
    v02_a  =  dz_a + uz_a; // New particle midpoint
    v03_a  = v00_a + ux_a;
    v04_a  = v01_a + uy_a;
    v05_a  = v02_a + uz_a; // New particle position

    outbnd_a = ( v03_a > one ) | ( v03_a < neg_one ) |
               ( v04_a > one ) | ( v04_a < neg_one ) |
               ( v05_a > one ) | ( v05_a < neg_one );

    v03_a  = merge( outbnd_a, dx_a, v03_a ); // Do not update outbnd particles
    v04_a  = merge( outbnd_a, dy_a, v04_a );
    v05_a  = merge( outbnd_a, dz_a, v05_a );

    //--------------------------------------------------------------------------
    // Interpolate fields for particles of block b.
    //--------------------------------------------------------------------------
    vp00_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 0) );
    vp01_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 1) );
    vp02_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 2) );
    vp03_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 3) );
    vp04_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 4) );
    vp05_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 5) );
    vp06_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 6) );
    vp07_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 7) );
    vp08_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 8) );
    vp09_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 9) );
    vp10_b = ( float * ALIGNED(16) ) ( f0 + ii_b(10) );
    vp11_b = ( float * ALIGNED(16) ) ( f0 + ii_b(11) );
    vp12_b = ( float * ALIGNED(16) ) ( f0 + ii_b(12) );
    vp13_b = ( float * ALIGNED(16) ) ( f0 + ii_b(13) );
    vp14_b = ( float * ALIGNED(16) ) ( f0 + ii_b(14) );
    vp15_b = ( float * ALIGNED(16) ) ( f0 + ii_b(15) );

    //--------------------------------------------------------------------------
    // Load data for 16 particles for block b.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00_b, vp01_b, vp02_b, vp03_b,
		   vp04_b, vp05_b, vp06_b, vp07_b,
		   vp08_b, vp09_b, vp10_b, vp11_b,
		   vp12_b, vp13_b, vp14_b, vp15_b,
		   hax_b, v00_b, v01_b, v02_b, hay_b, v03_b, v04_b, v05_b,
		   haz_b, v06_b, v07_b, v08_b, cbx_b, v09_b, cby_b, v10_b );

    //--------------------------------------------------------------------------
    // Process particle block b.
    //--------------------------------------------------------------------------
    hax_b = qdt_2mc*
            fma( fma( v02_b, dy_b, v01_b ), dz_b, fma( v00_b, dy_b, hax_b ) );

    hay_b = qdt_2mc*
            fma( fma( v05_b, dz_b, v04_b ), dx_b, fma( v03_b, dz_b, hay_b ) );

    haz_b = qdt_2mc*
            fma( fma( v08_b, dx_b, v07_b ), dy_b, fma( v06_b, dx_b, haz_b ) );

    cbx_b = fma( v09_b, dx_b, cbx_b );

    cby_b = fma( v10_b, dy_b, cby_b );

    load_16x2_tr( vp00_b+16, vp01_b+16, vp02_b+16, vp03_b+16,
		  vp04_b+16, vp05_b+16, vp06_b+16, vp07_b+16,
		  vp08_b+16, vp09_b+16, vp10_b+16, vp11_b+16,
		  vp12_b+16, vp13_b+16, vp14_b+16, vp15_b+16,
		  cbz_b, v05_b );

    cbz_b = fma( v05_b, dz_b, cbz_b );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    ux_b  += hax_b;
    uy_b  += hay_b;
    uz_b  += haz_b;
    v00_b  = qdt_2mc*rsqrt( one + fma( ux_b, ux_b,
                                       fma( uy_b, uy_b, uz_b*uz_b ) ) );
    v01_b  = fma( cbx_b, cbx_b, fma( cby_b, cby_b, cbz_b*cbz_b ) );
    v02_b  = (v00_b*v00_b)*v01_b;
    v03_b  = v00_b*fma( fma( two_fifteenths, v02_b, one_third ), v02_b, one );
    v04_b  = v03_b*rcp( fma( v03_b*v03_b, v01_b, one ) );
    v04_b += v04_b;
    v00_b  = fma( fms(  uy_b, cbz_b,  uz_b*cby_b ), v03_b, ux_b );
    v01_b  = fma( fms(  uz_b, cbx_b,  ux_b*cbz_b ), v03_b, uy_b );
    v02_b  = fma( fms(  ux_b, cby_b,  uy_b*cbx_b ), v03_b, uz_b );
    ux_b   = fma( fms( v01_b, cbz_b, v02_b*cby_b ), v04_b, ux_b );
    uy_b   = fma( fms( v02_b, cbx_b, v00_b*cbz_b ), v04_b, uy_b );
    uz_b   = fma( fms( v00_b, cby_b, v01_b*cbx_b ), v04_b, uz_b );
    ux_b  += hax_b;
    uy_b  += hay_b;
    uz_b  += haz_b;

    // Store ux, uy, uz in v6, v7, v8 so store_16x16_tr can be used below.
    v06_b  = ux_b;
    v07_b  = uy_b;
    v08_b  = uz_b;

    // Update the position of inbnd particles
    v00_b  = rsqrt( one + fma( ux_b, ux_b, fma( uy_b, uy_b, uz_b*uz_b ) ) );
    ux_b  *= cdt_dx;
    uy_b  *= cdt_dy;
    uz_b  *= cdt_dz;
    ux_b  *= v00_b;
    uy_b  *= v00_b;
    uz_b  *= v00_b;        // ux,uy,uz are normalized displ (relative to cell size)
    v00_b  =  dx_b + ux_b;
    v01_b  =  dy_b + uy_b;
    v02_b  =  dz_b + uz_b; // New particle midpoint
    v03_b  = v00_b + ux_b;
    v04_b  = v01_b + uy_b;
    v05_b  = v02_b + uz_b; // New particle position

    outbnd_b = ( v03_b > one ) | ( v03_b < neg_one ) |
               ( v04_b > one ) | ( v04_b < neg_one ) |
               ( v05_b > one ) | ( v05_b < neg_one );

    v03_b  = merge( outbnd_b, dx_b, v03_b ); // Do not update outbnd particles
    v04_b  = merge( outbnd_b, dy_b, v04_b );
    v05_b  = merge( outbnd_b, dz_b, v05_b );

    //--------------------------------------------------------------------------
    // Store results for 2 blocks of particles.
    //--------------------------------------------------------------------------
    store_16x16_tr_p( v03_a, v04_a, v05_a, ii_a, v06_a, v07_a, v08_a, q_a,
		      v03_b, v04_b, v05_b, ii_b, v06_b, v07_b, v08_b, q_b,
		      &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
		      &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		      &p[16].dx, &p[18].dx, &p[20].dx, &p[22].dx,
		      &p[24].dx, &p[26].dx, &p[28].dx, &p[30].dx );

    //--------------------------------------------------------------------------
    // Process particle block a.
    //--------------------------------------------------------------------------

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q_a  = czero( outbnd_a, q_a*qsp );   // Do not accumulate outbnd particles
    dx_a = v00_a;                        // Streak midpoint (valid for inbnd only)
    dy_a = v01_a;
    dz_a = v02_a;

    v13_a = q_a*ux_a*uy_a*uz_a*one_third;     // Charge conservation correction

    vp00_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 0) ); // Accumulator pointers
    vp01_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 1) );
    vp02_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 2) );
    vp03_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 3) );
    vp04_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 4) );
    vp05_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 5) );
    vp06_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 6) );
    vp07_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 7) );
    vp08_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 8) );
    vp09_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 9) );
    vp10_a = ( float * ALIGNED(16) ) ( a0 + ii_a(10) );
    vp11_a = ( float * ALIGNED(16) ) ( a0 + ii_a(11) );
    vp12_a = ( float * ALIGNED(16) ) ( a0 + ii_a(12) );
    vp13_a = ( float * ALIGNED(16) ) ( a0 + ii_a(13) );
    vp14_a = ( float * ALIGNED(16) ) ( a0 + ii_a(14) );
    vp15_a = ( float * ALIGNED(16) ) ( a0 + ii_a(15) );

    // Accumulate Jx for 16 particles into the v0-v3 vectors.
    v12_a  = q_a*ux_a;     // v12 = q ux
    v01_a  = v12_a*dy_a;   // v01 = q ux dy
    v00_a  = v12_a-v01_a;  // v00 = q ux (1-dy)
    v01_a += v12_a;        // v01 = q ux (1+dy)
    v12_a  = one+dz_a;     // v12 = 1+dz
    v02_a  = v00_a*v12_a;  // v02 = q ux (1-dy)(1+dz)
    v03_a  = v01_a*v12_a;  // v03 = q ux (1+dy)(1+dz)
    v12_a  = one-dz_a;     // v12 = 1-dz
    v00_a *= v12_a;        // v00 = q ux (1-dy)(1-dz)
    v01_a *= v12_a;        // v01 = q ux (1+dy)(1-dz)
    v00_a += v13_a;        // v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ]
    v01_a -= v13_a;        // v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ]
    v02_a -= v13_a;        // v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ]
    v03_a += v13_a;        // v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ]

    // Accumulate Jy for 16 particles into the v4-v7 vectors.
    v12_a  = q_a*uy_a;     // v12 = q uy
    v05_a  = v12_a*dz_a;   // v05 = q uy dz
    v04_a  = v12_a-v05_a;  // v04 = q uy (1-dz)
    v05_a += v12_a;        // v05 = q uy (1+dz)
    v12_a  = one+dx_a;     // v12 = 1+dx
    v06_a  = v04_a*v12_a;  // v06 = q uy (1-dz)(1+dx)
    v07_a  = v05_a*v12_a;  // v07 = q uy (1+dz)(1+dx)
    v12_a  = one-dx_a;     // v12 = 1-dx
    v04_a *= v12_a;        // v04 = q uy (1-dz)(1-dx)
    v05_a *= v12_a;        // v05 = q uy (1+dz)(1-dx)
    v04_a += v13_a;        // v04 = q uy [ (1-dz)(1-dx) + ux*uz/3 ]
    v05_a -= v13_a;        // v05 = q uy [ (1+dz)(1-dx) - ux*uz/3 ]
    v06_a -= v13_a;        // v06 = q uy [ (1-dz)(1+dx) - ux*uz/3 ]
    v07_a += v13_a;        // v07 = q uy [ (1+dz)(1+dx) + ux*uz/3 ]

    // Accumulate Jz for 16 particles into the v8-v11 vectors.
    v12_a  = q_a*uz_a;     // v12 = q uz
    v09_a  = v12_a*dx_a;   // v09 = q uz dx
    v08_a  = v12_a-v09_a;  // v08 = q uz (1-dx)
    v09_a += v12_a;        // v09 = q uz (1+dx)
    v12_a  = one+dy_a;     // v12 = 1+dy
    v10_a  = v08_a*v12_a;  // v10 = q uz (1-dx)(1+dy)
    v11_a  = v09_a*v12_a;  // v11 = q uz (1+dx)(1+dy)
    v12_a  = one-dy_a;     // v12 = 1-dy
    v08_a *= v12_a;        // v08 = q uz (1-dx)(1-dy)
    v09_a *= v12_a;        // v09 = q uz (1+dx)(1-dy)
    v08_a += v13_a;        // v08 = q uz [ (1-dx)(1-dy) + ux*uy/3 ]
    v09_a -= v13_a;        // v09 = q uz [ (1+dx)(1-dy) - ux*uy/3 ]
    v10_a -= v13_a;        // v10 = q uz [ (1-dx)(1+dy) - ux*uy/3 ]
    v11_a += v13_a;        // v11 = q uz [ (1+dx)(1+dy) + ux*uy/3 ]

    // Zero the v12-v15 vectors prior to transposing the data.
    v12_a = 0.0;
    v13_a = 0.0;
    v14_a = 0.0;
    v15_a = 0.0;

    // Transpose the data in vectors v0-v15 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v00_a, v01_a, v02_a, v03_a, v04_a, v05_a, v06_a, v07_a,
	       v08_a, v09_a, v10_a, v11_a, v12_a, v13_a, v14_a, v15_a );

    // Add the contributions to Jx, Jy and Jz from 16 particles into the
    // accumulator arrays for Jx, Jy and Jz.
    increment_16x1( vp00_a, v00_a );
    increment_16x1( vp01_a, v01_a );
    increment_16x1( vp02_a, v02_a );
    increment_16x1( vp03_a, v03_a );
    increment_16x1( vp04_a, v04_a );
    increment_16x1( vp05_a, v05_a );
    increment_16x1( vp06_a, v06_a );
    increment_16x1( vp07_a, v07_a );
    increment_16x1( vp08_a, v08_a );
    increment_16x1( vp09_a, v09_a );
    increment_16x1( vp10_a, v10_a );
    increment_16x1( vp11_a, v11_a );
    increment_16x1( vp12_a, v12_a );
    increment_16x1( vp13_a, v13_a );
    increment_16x1( vp14_a, v14_a );
    increment_16x1( vp15_a, v15_a );

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
    if ( outbnd_a(N) )                              /* Unlikely */      \
    {                                                                   \
      local_pm->dispx = ux_a(N);                                        \
      local_pm->dispy = uy_a(N);                                        \
      local_pm->dispz = uz_a(N);                                        \
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

    //--------------------------------------------------------------------------
    // Process particle block b.
    //--------------------------------------------------------------------------

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q_b  = czero( outbnd_b, q_b*qsp );   // Do not accumulate outbnd particles
    dx_b = v00_b;                        // Streak midpoint (valid for inbnd only)
    dy_b = v01_b;
    dz_b = v02_b;

    v13_b = q_b*ux_b*uy_b*uz_b*one_third;     // Charge conservation correction

    vp00_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 0) ); // Accumulator pointers
    vp01_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 1) );
    vp02_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 2) );
    vp03_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 3) );
    vp04_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 4) );
    vp05_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 5) );
    vp06_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 6) );
    vp07_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 7) );
    vp08_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 8) );
    vp09_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 9) );
    vp10_b = ( float * ALIGNED(16) ) ( a0 + ii_b(10) );
    vp11_b = ( float * ALIGNED(16) ) ( a0 + ii_b(11) );
    vp12_b = ( float * ALIGNED(16) ) ( a0 + ii_b(12) );
    vp13_b = ( float * ALIGNED(16) ) ( a0 + ii_b(13) );
    vp14_b = ( float * ALIGNED(16) ) ( a0 + ii_b(14) );
    vp15_b = ( float * ALIGNED(16) ) ( a0 + ii_b(15) );

    // Accumulate Jx for 16 particles into the v0-v3 vectors.
    v12_b  = q_b*ux_b;     // v12 = q ux
    v01_b  = v12_b*dy_b;   // v01 = q ux dy
    v00_b  = v12_b-v01_b;  // v00 = q ux (1-dy)
    v01_b += v12_b;        // v01 = q ux (1+dy)
    v12_b  = one+dz_b;     // v12 = 1+dz
    v02_b  = v00_b*v12_b;  // v02 = q ux (1-dy)(1+dz)
    v03_b  = v01_b*v12_b;  // v03 = q ux (1+dy)(1+dz)
    v12_b  = one-dz_b;     // v12 = 1-dz
    v00_b *= v12_b;        // v00 = q ux (1-dy)(1-dz)
    v01_b *= v12_b;        // v01 = q ux (1+dy)(1-dz)
    v00_b += v13_b;        // v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ]
    v01_b -= v13_b;        // v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ]
    v02_b -= v13_b;        // v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ]
    v03_b += v13_b;        // v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ]

    // Accumulate Jy for 16 particles into the v4-v7 vectors.
    v12_b  = q_b*uy_b;     // v12 = q uy
    v05_b  = v12_b*dz_b;   // v05 = q uy dz
    v04_b  = v12_b-v05_b;  // v04 = q uy (1-dz)
    v05_b += v12_b;        // v05 = q uy (1+dz)
    v12_b  = one+dx_b;     // v12 = 1+dx
    v06_b  = v04_b*v12_b;  // v06 = q uy (1-dz)(1+dx)
    v07_b  = v05_b*v12_b;  // v07 = q uy (1+dz)(1+dx)
    v12_b  = one-dx_b;     // v12 = 1-dx
    v04_b *= v12_b;        // v04 = q uy (1-dz)(1-dx)
    v05_b *= v12_b;        // v05 = q uy (1+dz)(1-dx)
    v04_b += v13_b;        // v04 = q uy [ (1-dz)(1-dx) + ux*uz/3 ]
    v05_b -= v13_b;        // v05 = q uy [ (1+dz)(1-dx) - ux*uz/3 ]
    v06_b -= v13_b;        // v06 = q uy [ (1-dz)(1+dx) - ux*uz/3 ]
    v07_b += v13_b;        // v07 = q uy [ (1+dz)(1+dx) + ux*uz/3 ]

    // Accumulate Jz for 16 particles into the v8-v11 vectors.
    v12_b  = q_b*uz_b;     // v12 = q uz
    v09_b  = v12_b*dx_b;   // v09 = q uz dx
    v08_b  = v12_b-v09_b;  // v08 = q uz (1-dx)
    v09_b += v12_b;        // v09 = q uz (1+dx)
    v12_b  = one+dy_b;     // v12 = 1+dy
    v10_b  = v08_b*v12_b;  // v10 = q uz (1-dx)(1+dy)
    v11_b  = v09_b*v12_b;  // v11 = q uz (1+dx)(1+dy)
    v12_b  = one-dy_b;     // v12 = 1-dy
    v08_b *= v12_b;        // v08 = q uz (1-dx)(1-dy)
    v09_b *= v12_b;        // v09 = q uz (1+dx)(1-dy)
    v08_b += v13_b;        // v08 = q uz [ (1-dx)(1-dy) + ux*uy/3 ]
    v09_b -= v13_b;        // v09 = q uz [ (1+dx)(1-dy) - ux*uy/3 ]
    v10_b -= v13_b;        // v10 = q uz [ (1-dx)(1+dy) - ux*uy/3 ]
    v11_b += v13_b;        // v11 = q uz [ (1+dx)(1+dy) + ux*uy/3 ]

    // Zero the v12-v15 vectors prior to transposing the data.
    v12_b = 0.0;
    v13_b = 0.0;
    v14_b = 0.0;
    v15_b = 0.0;

    // Transpose the data in vectors v0-v15 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v00_b, v01_b, v02_b, v03_b, v04_b, v05_b, v06_b, v07_b,
	       v08_b, v09_b, v10_b, v11_b, v12_b, v13_b, v14_b, v15_b );

    // Add the contributions to Jx, Jy and Jz from 16 particles into the
    // accumulator arrays for Jx, Jy and Jz.
    increment_16x1( vp00_b, v00_b );
    increment_16x1( vp01_b, v01_b );
    increment_16x1( vp02_b, v02_b );
    increment_16x1( vp03_b, v03_b );
    increment_16x1( vp04_b, v04_b );
    increment_16x1( vp05_b, v05_b );
    increment_16x1( vp06_b, v06_b );
    increment_16x1( vp07_b, v07_b );
    increment_16x1( vp08_b, v08_b );
    increment_16x1( vp09_b, v09_b );
    increment_16x1( vp10_b, v10_b );
    increment_16x1( vp11_b, v11_b );
    increment_16x1( vp12_b, v12_b );
    increment_16x1( vp13_b, v13_b );
    increment_16x1( vp14_b, v14_b );
    increment_16x1( vp15_b, v15_b );

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
    if ( outbnd_b(N) )                              /* Unlikely */      \
    {                                                                   \
      local_pm->dispx = ux_b(N);                                        \
      local_pm->dispy = uy_b(N);                                        \
      local_pm->dispz = uz_b(N);                                        \
      local_pm->i     = ( p - p0 ) + N + 16;                            \
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

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
#endif // Method 1

#if 0
//----------------------------------------------------------------------------//
// Method 2
//----------------------------------------------------------------------------//
// This method processes the particles in the same order as the reference
// implementation and gives good reproducibility. This is achieved by an
// extra reordering step after the load and transpose step and before the
// transpose and store step.  A better approach is to modify the transpose
// step and put the particles in the correct order in a single step instead
// of two steps.
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

  float                * ALIGNED(16)  vp00_a;
  float                * ALIGNED(16)  vp01_a;
  float                * ALIGNED(16)  vp02_a;
  float                * ALIGNED(16)  vp03_a;
  float                * ALIGNED(16)  vp04_a;
  float                * ALIGNED(16)  vp05_a;
  float                * ALIGNED(16)  vp06_a;
  float                * ALIGNED(16)  vp07_a;
  float                * ALIGNED(16)  vp08_a;
  float                * ALIGNED(16)  vp09_a;
  float                * ALIGNED(16)  vp10_a;
  float                * ALIGNED(16)  vp11_a;
  float                * ALIGNED(16)  vp12_a;
  float                * ALIGNED(16)  vp13_a;
  float                * ALIGNED(16)  vp14_a;
  float                * ALIGNED(16)  vp15_a;

  float                * ALIGNED(16)  vp00_b;
  float                * ALIGNED(16)  vp01_b;
  float                * ALIGNED(16)  vp02_b;
  float                * ALIGNED(16)  vp03_b;
  float                * ALIGNED(16)  vp04_b;
  float                * ALIGNED(16)  vp05_b;
  float                * ALIGNED(16)  vp06_b;
  float                * ALIGNED(16)  vp07_b;
  float                * ALIGNED(16)  vp08_b;
  float                * ALIGNED(16)  vp09_b;
  float                * ALIGNED(16)  vp10_b;
  float                * ALIGNED(16)  vp11_b;
  float                * ALIGNED(16)  vp12_b;
  float                * ALIGNED(16)  vp13_b;
  float                * ALIGNED(16)  vp14_b;
  float                * ALIGNED(16)  vp15_b;

  const v16float qdt_2mc(args->qdt_2mc);
  const v16float cdt_dx(args->cdt_dx);
  const v16float cdt_dy(args->cdt_dy);
  const v16float cdt_dz(args->cdt_dz);
  const v16float qsp(args->qsp);
  const v16float one(1.);
  const v16float one_third(1./3.);
  const v16float two_fifteenths(2./15.);
  const v16float neg_one(-1.);

  const float _qsp = args->qsp;

  v16float dx_aa, dy_aa, dz_aa, ux_aa, uy_aa, uz_aa, q_aa;
  v16float dx_a, dy_a, dz_a, ux_a, uy_a, uz_a, q_a;
  v16float hax_a, hay_a, haz_a, cbx_a, cby_a, cbz_a;
  v16float v00_a, v01_a, v02_a, v03_a, v04_a, v05_a, v06_a, v07_a;
  v16float v08_a, v09_a, v10_a, v11_a, v12_a, v13_a, v14_a, v15_a;
  v16int   ii_aa, ii_a, outbnd_a;

  v16float dx_bb, dy_bb, dz_bb, ux_bb, uy_bb, uz_bb, q_bb;
  v16float dx_b, dy_b, dz_b, ux_b, uy_b, uz_b, q_b;
  v16float hax_b, hay_b, haz_b, cbx_b, cby_b, cbz_b;
  v16float v00_b, v01_b, v02_b, v03_b, v04_b, v05_b, v06_b, v07_b;
  v16float v08_b, v09_b, v10_b, v11_b, v12_b, v13_b, v14_b, v15_b;
  v16int   ii_bb, ii_b, outbnd_b;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 32, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 5;

  // Determine which movers are reserved for this pipeline.
  // Movers (16 bytes) should be reserved for pipelines in at least
  // multiples of 8 such that the set of particle movers reserved for
  // a pipeline is 128-byte aligned and a multiple of 128-byte in
  // size.  The host is guaranteed to get enough movers to process its
  // particles with this allocation.

  max_nm = args->max_nm - ( args->np&15 );

  if ( max_nm < 0 ) max_nm = 0;

  DISTRIBUTE( max_nm, 8, pipeline_rank, n_pipeline, itmp, max_nm );

  if ( pipeline_rank==n_pipeline ) max_nm = args->max_nm - itmp;

  pm   = args->pm + itmp;
  nm   = 0;
  itmp = 0;

  // Determine which accumulator array to use.
  // The host gets the first accumulator array.

  a0 += ( 1 + pipeline_rank ) *
        POW2_CEIL( (args->nx+2)*(args->ny+2)*(args->nz+2), 2 );

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=32 )
  {
    load_16x16_tr( &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
		   &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		   &p[16].dx, &p[18].dx, &p[20].dx, &p[22].dx,
		   &p[24].dx, &p[26].dx, &p[28].dx, &p[30].dx,
		   dx_aa, dy_aa, dz_aa, ii_aa, ux_aa, uy_aa, uz_aa, q_aa,
                   dx_bb, dy_bb, dz_bb, ii_bb, ux_bb, uy_bb, uz_bb, q_bb );

    // Properly order the particles in the vector arrays.
    for( int j = 0; j < 8; j++ )
    {
      dx_a[2*j  ] = dx_aa[j];
      dx_a[2*j+1] = dx_bb[j];

      dx_b[2*j  ] = dx_aa[j+8];
      dx_b[2*j+1] = dx_bb[j+8];

      dy_a[2*j  ] = dy_aa[j];
      dy_a[2*j+1] = dy_bb[j];

      dy_b[2*j  ] = dy_aa[j+8];
      dy_b[2*j+1] = dy_bb[j+8];

      dz_a[2*j  ] = dz_aa[j];
      dz_a[2*j+1] = dz_bb[j];

      dz_b[2*j  ] = dz_aa[j+8];
      dz_b[2*j+1] = dz_bb[j+8];

      ii_a[2*j  ] = ii_aa[j];
      ii_a[2*j+1] = ii_bb[j];

      ii_b[2*j  ] = ii_aa[j+8];
      ii_b[2*j+1] = ii_bb[j+8];

      ux_a[2*j  ] = ux_aa[j];
      ux_a[2*j+1] = ux_bb[j];

      ux_b[2*j  ] = ux_aa[j+8];
      ux_b[2*j+1] = ux_bb[j+8];

      uy_a[2*j  ] = uy_aa[j];
      uy_a[2*j+1] = uy_bb[j];

      uy_b[2*j  ] = uy_aa[j+8];
      uy_b[2*j+1] = uy_bb[j+8];

      uz_a[2*j  ] = uz_aa[j];
      uz_a[2*j+1] = uz_bb[j];

      uz_b[2*j  ] = uz_aa[j+8];
      uz_b[2*j+1] = uz_bb[j+8];

       q_a[2*j  ] =  q_aa[j];
       q_a[2*j+1] =  q_bb[j];

       q_b[2*j  ] =  q_aa[j+8];
       q_b[2*j+1] =  q_bb[j+8];
    }

    //--------------------------------------------------------------------------
    // Interpolate fields for particles of block a.
    //--------------------------------------------------------------------------
    vp00_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 0) );
    vp01_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 1) );
    vp02_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 2) );
    vp03_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 3) );
    vp04_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 4) );
    vp05_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 5) );
    vp06_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 6) );
    vp07_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 7) );
    vp08_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 8) );
    vp09_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 9) );
    vp10_a = ( float * ALIGNED(16) ) ( f0 + ii_a(10) );
    vp11_a = ( float * ALIGNED(16) ) ( f0 + ii_a(11) );
    vp12_a = ( float * ALIGNED(16) ) ( f0 + ii_a(12) );
    vp13_a = ( float * ALIGNED(16) ) ( f0 + ii_a(13) );
    vp14_a = ( float * ALIGNED(16) ) ( f0 + ii_a(14) );
    vp15_a = ( float * ALIGNED(16) ) ( f0 + ii_a(15) );

    //--------------------------------------------------------------------------
    // Load data for 16 particles for block a.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00_a, vp01_a, vp02_a, vp03_a,
		   vp04_a, vp05_a, vp06_a, vp07_a,
		   vp08_a, vp09_a, vp10_a, vp11_a,
		   vp12_a, vp13_a, vp14_a, vp15_a,
		   hax_a, v00_a, v01_a, v02_a, hay_a, v03_a, v04_a, v05_a,
		   haz_a, v06_a, v07_a, v08_a, cbx_a, v09_a, cby_a, v10_a );

    //--------------------------------------------------------------------------
    // Process particle block a.
    //--------------------------------------------------------------------------
    hax_a = qdt_2mc*
            fma( fma( v02_a, dy_a, v01_a ), dz_a, fma( v00_a, dy_a, hax_a ) );

    hay_a = qdt_2mc*
            fma( fma( v05_a, dz_a, v04_a ), dx_a, fma( v03_a, dz_a, hay_a ) );

    haz_a = qdt_2mc*
            fma( fma( v08_a, dx_a, v07_a ), dy_a, fma( v06_a, dx_a, haz_a ) );

    cbx_a = fma( v09_a, dx_a, cbx_a );

    cby_a = fma( v10_a, dy_a, cby_a );

    load_16x2_tr( vp00_a+16, vp01_a+16, vp02_a+16, vp03_a+16,
		  vp04_a+16, vp05_a+16, vp06_a+16, vp07_a+16,
		  vp08_a+16, vp09_a+16, vp10_a+16, vp11_a+16,
		  vp12_a+16, vp13_a+16, vp14_a+16, vp15_a+16,
		  cbz_a, v05_a );

    cbz_a = fma( v05_a, dz_a, cbz_a );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    ux_a += hax_a;
    uy_a += hay_a;
    uz_a += haz_a;
    v00_a  = qdt_2mc*rsqrt( one + fma( ux_a, ux_a,
                                       fma( uy_a, uy_a, uz_a*uz_a ) ) );
    v01_a  = fma( cbx_a, cbx_a, fma( cby_a, cby_a, cbz_a*cbz_a ) );
    v02_a  = (v00_a*v00_a)*v01_a;
    v03_a  = v00_a*fma( fma( two_fifteenths, v02_a, one_third ), v02_a, one );
    v04_a  = v03_a*rcp( fma( v03_a*v03_a, v01_a, one ) );
    v04_a += v04_a;
    v00_a  = fma( fms(  uy_a, cbz_a,  uz_a*cby_a ), v03_a, ux_a );
    v01_a  = fma( fms(  uz_a, cbx_a,  ux_a*cbz_a ), v03_a, uy_a );
    v02_a  = fma( fms(  ux_a, cby_a,  uy_a*cbx_a ), v03_a, uz_a );
    ux_a   = fma( fms( v01_a, cbz_a, v02_a*cby_a ), v04_a, ux_a );
    uy_a   = fma( fms( v02_a, cbx_a, v00_a*cbz_a ), v04_a, uy_a );
    uz_a   = fma( fms( v00_a, cby_a, v01_a*cbx_a ), v04_a, uz_a );
    ux_a  += hax_a;
    uy_a  += hay_a;
    uz_a  += haz_a;

    // Store ux, uy, uz in v6, v7, v8 so store_16x16_tr can be used below.
    v06_a  = ux_a;
    v07_a  = uy_a;
    v08_a  = uz_a;

    // Update the position of inbnd particles
    v00_a  = rsqrt( one + fma( ux_a, ux_a, fma( uy_a, uy_a, uz_a*uz_a ) ) );
    ux_a  *= cdt_dx;
    uy_a  *= cdt_dy;
    uz_a  *= cdt_dz;
    ux_a  *= v00_a;
    uy_a  *= v00_a;
    uz_a  *= v00_a;        // ux,uy,uz are normalized displ (relative to cell size)
    v00_a  =  dx_a + ux_a;
    v01_a  =  dy_a + uy_a;
    v02_a  =  dz_a + uz_a; // New particle midpoint
    v03_a  = v00_a + ux_a;
    v04_a  = v01_a + uy_a;
    v05_a  = v02_a + uz_a; // New particle position

    outbnd_a = ( v03_a > one ) | ( v03_a < neg_one ) |
               ( v04_a > one ) | ( v04_a < neg_one ) |
               ( v05_a > one ) | ( v05_a < neg_one );

    v03_a  = merge( outbnd_a, dx_a, v03_a ); // Do not update outbnd particles
    v04_a  = merge( outbnd_a, dy_a, v04_a );
    v05_a  = merge( outbnd_a, dz_a, v05_a );

    //--------------------------------------------------------------------------
    // Interpolate fields for particles of block b.
    //--------------------------------------------------------------------------
    vp00_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 0) );
    vp01_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 1) );
    vp02_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 2) );
    vp03_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 3) );
    vp04_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 4) );
    vp05_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 5) );
    vp06_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 6) );
    vp07_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 7) );
    vp08_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 8) );
    vp09_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 9) );
    vp10_b = ( float * ALIGNED(16) ) ( f0 + ii_b(10) );
    vp11_b = ( float * ALIGNED(16) ) ( f0 + ii_b(11) );
    vp12_b = ( float * ALIGNED(16) ) ( f0 + ii_b(12) );
    vp13_b = ( float * ALIGNED(16) ) ( f0 + ii_b(13) );
    vp14_b = ( float * ALIGNED(16) ) ( f0 + ii_b(14) );
    vp15_b = ( float * ALIGNED(16) ) ( f0 + ii_b(15) );

    //--------------------------------------------------------------------------
    // Load data for 16 particles for block b.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00_b, vp01_b, vp02_b, vp03_b,
		   vp04_b, vp05_b, vp06_b, vp07_b,
		   vp08_b, vp09_b, vp10_b, vp11_b,
		   vp12_b, vp13_b, vp14_b, vp15_b,
		   hax_b, v00_b, v01_b, v02_b, hay_b, v03_b, v04_b, v05_b,
		   haz_b, v06_b, v07_b, v08_b, cbx_b, v09_b, cby_b, v10_b );

    //--------------------------------------------------------------------------
    // Process particle block b.
    //--------------------------------------------------------------------------
    hax_b = qdt_2mc*
            fma( fma( v02_b, dy_b, v01_b ), dz_b, fma( v00_b, dy_b, hax_b ) );

    hay_b = qdt_2mc*
            fma( fma( v05_b, dz_b, v04_b ), dx_b, fma( v03_b, dz_b, hay_b ) );

    haz_b = qdt_2mc*
            fma( fma( v08_b, dx_b, v07_b ), dy_b, fma( v06_b, dx_b, haz_b ) );

    cbx_b = fma( v09_b, dx_b, cbx_b );

    cby_b = fma( v10_b, dy_b, cby_b );

    load_16x2_tr( vp00_b+16, vp01_b+16, vp02_b+16, vp03_b+16,
		  vp04_b+16, vp05_b+16, vp06_b+16, vp07_b+16,
		  vp08_b+16, vp09_b+16, vp10_b+16, vp11_b+16,
		  vp12_b+16, vp13_b+16, vp14_b+16, vp15_b+16,
		  cbz_b, v05_b );

    cbz_b = fma( v05_b, dz_b, cbz_b );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    ux_b  += hax_b;
    uy_b  += hay_b;
    uz_b  += haz_b;
    v00_b  = qdt_2mc*rsqrt( one + fma( ux_b, ux_b,
                                       fma( uy_b, uy_b, uz_b*uz_b ) ) );
    v01_b  = fma( cbx_b, cbx_b, fma( cby_b, cby_b, cbz_b*cbz_b ) );
    v02_b  = (v00_b*v00_b)*v01_b;
    v03_b  = v00_b*fma( fma( two_fifteenths, v02_b, one_third ), v02_b, one );
    v04_b  = v03_b*rcp( fma( v03_b*v03_b, v01_b, one ) );
    v04_b += v04_b;
    v00_b  = fma( fms(  uy_b, cbz_b,  uz_b*cby_b ), v03_b, ux_b );
    v01_b  = fma( fms(  uz_b, cbx_b,  ux_b*cbz_b ), v03_b, uy_b );
    v02_b  = fma( fms(  ux_b, cby_b,  uy_b*cbx_b ), v03_b, uz_b );
    ux_b   = fma( fms( v01_b, cbz_b, v02_b*cby_b ), v04_b, ux_b );
    uy_b   = fma( fms( v02_b, cbx_b, v00_b*cbz_b ), v04_b, uy_b );
    uz_b   = fma( fms( v00_b, cby_b, v01_b*cbx_b ), v04_b, uz_b );
    ux_b  += hax_b;
    uy_b  += hay_b;
    uz_b  += haz_b;

    // Store ux, uy, uz in v6, v7, v8 so store_16x16_tr can be used below.
    v06_b  = ux_b;
    v07_b  = uy_b;
    v08_b  = uz_b;

    // Update the position of inbnd particles
    v00_b  = rsqrt( one + fma( ux_b, ux_b, fma( uy_b, uy_b, uz_b*uz_b ) ) );
    ux_b  *= cdt_dx;
    uy_b  *= cdt_dy;
    uz_b  *= cdt_dz;
    ux_b  *= v00_b;
    uy_b  *= v00_b;
    uz_b  *= v00_b;        // ux,uy,uz are normalized displ (relative to cell size)
    v00_b  =  dx_b + ux_b;
    v01_b  =  dy_b + uy_b;
    v02_b  =  dz_b + uz_b; // New particle midpoint
    v03_b  = v00_b + ux_b;
    v04_b  = v01_b + uy_b;
    v05_b  = v02_b + uz_b; // New particle position

    outbnd_b = ( v03_b > one ) | ( v03_b < neg_one ) |
               ( v04_b > one ) | ( v04_b < neg_one ) |
               ( v05_b > one ) | ( v05_b < neg_one );

    v03_b  = merge( outbnd_b, dx_b, v03_b ); // Do not update outbnd particles
    v04_b  = merge( outbnd_b, dy_b, v04_b );
    v05_b  = merge( outbnd_b, dz_b, v05_b );

    // Properly order the particles in the vector arrays.
    for( int j = 0; j < 8; j++ )
    {
      dx_aa[j  ] = v03_a[2*j  ];
      dx_aa[j+8] = v03_b[2*j  ];

      dx_bb[j  ] = v03_a[2*j+1];
      dx_bb[j+8] = v03_b[2*j+1];

      dy_aa[j  ] = v04_a[2*j  ];
      dy_aa[j+8] = v04_b[2*j  ];

      dy_bb[j  ] = v04_a[2*j+1];
      dy_bb[j+8] = v04_b[2*j+1];

      dz_aa[j  ] = v05_a[2*j  ];
      dz_aa[j+8] = v05_b[2*j  ];

      dz_bb[j  ] = v05_a[2*j+1];
      dz_bb[j+8] = v05_b[2*j+1];

      ii_aa[j  ] = ii_a[2*j  ];
      ii_aa[j+8] = ii_b[2*j  ];

      ii_bb[j  ] = ii_a[2*j+1];
      ii_bb[j+8] = ii_b[2*j+1];

      ux_aa[j  ] = v06_a[2*j  ];
      ux_aa[j+8] = v06_b[2*j  ];

      ux_bb[j  ] = v06_a[2*j+1];
      ux_bb[j+8] = v06_b[2*j+1];

      uy_aa[j  ] = v07_a[2*j  ];
      uy_aa[j+8] = v07_b[2*j  ];

      uy_bb[j  ] = v07_a[2*j+1];
      uy_bb[j+8] = v07_b[2*j+1];

      uz_aa[j  ] = v08_a[2*j  ];
      uz_aa[j+8] = v08_b[2*j  ];

      uz_bb[j  ] = v08_a[2*j+1];
      uz_bb[j+8] = v08_b[2*j+1];

       q_aa[j  ] =  q_a[2*j  ];
       q_aa[j+8] =  q_b[2*j  ];

       q_bb[j  ] =  q_a[2*j+1];
       q_bb[j+8] =  q_b[2*j+1];
    }

    //--------------------------------------------------------------------------
    // Store results for 2 blocks of particles.
    //--------------------------------------------------------------------------
    store_16x16_tr( dx_aa, dy_aa, dz_aa, ii_aa, ux_aa, uy_aa, uz_aa, q_aa,
		    dx_bb, dy_bb, dz_bb, ii_bb, ux_bb, uy_bb, uz_bb, q_bb,
		    &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
		    &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		    &p[16].dx, &p[18].dx, &p[20].dx, &p[22].dx,
		    &p[24].dx, &p[26].dx, &p[28].dx, &p[30].dx );

    //--------------------------------------------------------------------------
    // Process particle block a.
    //--------------------------------------------------------------------------

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q_a  = czero( outbnd_a, q_a*qsp );   // Do not accumulate outbnd particles
    dx_a = v00_a;                        // Streak midpoint (valid for inbnd only)
    dy_a = v01_a;
    dz_a = v02_a;

    v13_a = q_a*ux_a*uy_a*uz_a*one_third;     // Charge conservation correction

    vp00_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 0) ); // Accumulator pointers
    vp01_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 1) );
    vp02_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 2) );
    vp03_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 3) );
    vp04_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 4) );
    vp05_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 5) );
    vp06_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 6) );
    vp07_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 7) );
    vp08_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 8) );
    vp09_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 9) );
    vp10_a = ( float * ALIGNED(16) ) ( a0 + ii_a(10) );
    vp11_a = ( float * ALIGNED(16) ) ( a0 + ii_a(11) );
    vp12_a = ( float * ALIGNED(16) ) ( a0 + ii_a(12) );
    vp13_a = ( float * ALIGNED(16) ) ( a0 + ii_a(13) );
    vp14_a = ( float * ALIGNED(16) ) ( a0 + ii_a(14) );
    vp15_a = ( float * ALIGNED(16) ) ( a0 + ii_a(15) );

    // Accumulate Jx for 16 particles into the v0-v3 vectors.
    v12_a  = q_a*ux_a;     // v12 = q ux
    v01_a  = v12_a*dy_a;   // v01 = q ux dy
    v00_a  = v12_a-v01_a;  // v00 = q ux (1-dy)
    v01_a += v12_a;        // v01 = q ux (1+dy)
    v12_a  = one+dz_a;     // v12 = 1+dz
    v02_a  = v00_a*v12_a;  // v02 = q ux (1-dy)(1+dz)
    v03_a  = v01_a*v12_a;  // v03 = q ux (1+dy)(1+dz)
    v12_a  = one-dz_a;     // v12 = 1-dz
    v00_a *= v12_a;        // v00 = q ux (1-dy)(1-dz)
    v01_a *= v12_a;        // v01 = q ux (1+dy)(1-dz)
    v00_a += v13_a;        // v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ]
    v01_a -= v13_a;        // v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ]
    v02_a -= v13_a;        // v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ]
    v03_a += v13_a;        // v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ]

    // Accumulate Jy for 16 particles into the v4-v7 vectors.
    v12_a  = q_a*uy_a;     // v12 = q uy
    v05_a  = v12_a*dz_a;   // v05 = q uy dz
    v04_a  = v12_a-v05_a;  // v04 = q uy (1-dz)
    v05_a += v12_a;        // v05 = q uy (1+dz)
    v12_a  = one+dx_a;     // v12 = 1+dx
    v06_a  = v04_a*v12_a;  // v06 = q uy (1-dz)(1+dx)
    v07_a  = v05_a*v12_a;  // v07 = q uy (1+dz)(1+dx)
    v12_a  = one-dx_a;     // v12 = 1-dx
    v04_a *= v12_a;        // v04 = q uy (1-dz)(1-dx)
    v05_a *= v12_a;        // v05 = q uy (1+dz)(1-dx)
    v04_a += v13_a;        // v04 = q uy [ (1-dz)(1-dx) + ux*uz/3 ]
    v05_a -= v13_a;        // v05 = q uy [ (1+dz)(1-dx) - ux*uz/3 ]
    v06_a -= v13_a;        // v06 = q uy [ (1-dz)(1+dx) - ux*uz/3 ]
    v07_a += v13_a;        // v07 = q uy [ (1+dz)(1+dx) + ux*uz/3 ]

    // Accumulate Jz for 16 particles into the v8-v11 vectors.
    v12_a  = q_a*uz_a;     // v12 = q uz
    v09_a  = v12_a*dx_a;   // v09 = q uz dx
    v08_a  = v12_a-v09_a;  // v08 = q uz (1-dx)
    v09_a += v12_a;        // v09 = q uz (1+dx)
    v12_a  = one+dy_a;     // v12 = 1+dy
    v10_a  = v08_a*v12_a;  // v10 = q uz (1-dx)(1+dy)
    v11_a  = v09_a*v12_a;  // v11 = q uz (1+dx)(1+dy)
    v12_a  = one-dy_a;     // v12 = 1-dy
    v08_a *= v12_a;        // v08 = q uz (1-dx)(1-dy)
    v09_a *= v12_a;        // v09 = q uz (1+dx)(1-dy)
    v08_a += v13_a;        // v08 = q uz [ (1-dx)(1-dy) + ux*uy/3 ]
    v09_a -= v13_a;        // v09 = q uz [ (1+dx)(1-dy) - ux*uy/3 ]
    v10_a -= v13_a;        // v10 = q uz [ (1-dx)(1+dy) - ux*uy/3 ]
    v11_a += v13_a;        // v11 = q uz [ (1+dx)(1+dy) + ux*uy/3 ]

    // Zero the v12-v15 vectors prior to transposing the data.
    v12_a = 0.0;
    v13_a = 0.0;
    v14_a = 0.0;
    v15_a = 0.0;

    // Transpose the data in vectors v0-v15 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v00_a, v01_a, v02_a, v03_a, v04_a, v05_a, v06_a, v07_a,
	       v08_a, v09_a, v10_a, v11_a, v12_a, v13_a, v14_a, v15_a );

    // Add the contributions to Jx, Jy and Jz from 16 particles into the
    // accumulator arrays for Jx, Jy and Jz.
    increment_16x1( vp00_a, v00_a );
    increment_16x1( vp01_a, v01_a );
    increment_16x1( vp02_a, v02_a );
    increment_16x1( vp03_a, v03_a );
    increment_16x1( vp04_a, v04_a );
    increment_16x1( vp05_a, v05_a );
    increment_16x1( vp06_a, v06_a );
    increment_16x1( vp07_a, v07_a );
    increment_16x1( vp08_a, v08_a );
    increment_16x1( vp09_a, v09_a );
    increment_16x1( vp10_a, v10_a );
    increment_16x1( vp11_a, v11_a );
    increment_16x1( vp12_a, v12_a );
    increment_16x1( vp13_a, v13_a );
    increment_16x1( vp14_a, v14_a );
    increment_16x1( vp15_a, v15_a );

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
    if ( outbnd_a(N) )                              /* Unlikely */      \
    {                                                                   \
      local_pm->dispx = ux_a(N);                                        \
      local_pm->dispy = uy_a(N);                                        \
      local_pm->dispz = uz_a(N);                                        \
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

    //--------------------------------------------------------------------------
    // Process particle block b.
    //--------------------------------------------------------------------------

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q_b  = czero( outbnd_b, q_b*qsp );   // Do not accumulate outbnd particles
    dx_b = v00_b;                        // Streak midpoint (valid for inbnd only)
    dy_b = v01_b;
    dz_b = v02_b;

    v13_b = q_b*ux_b*uy_b*uz_b*one_third;     // Charge conservation correction

    vp00_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 0) ); // Accumulator pointers
    vp01_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 1) );
    vp02_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 2) );
    vp03_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 3) );
    vp04_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 4) );
    vp05_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 5) );
    vp06_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 6) );
    vp07_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 7) );
    vp08_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 8) );
    vp09_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 9) );
    vp10_b = ( float * ALIGNED(16) ) ( a0 + ii_b(10) );
    vp11_b = ( float * ALIGNED(16) ) ( a0 + ii_b(11) );
    vp12_b = ( float * ALIGNED(16) ) ( a0 + ii_b(12) );
    vp13_b = ( float * ALIGNED(16) ) ( a0 + ii_b(13) );
    vp14_b = ( float * ALIGNED(16) ) ( a0 + ii_b(14) );
    vp15_b = ( float * ALIGNED(16) ) ( a0 + ii_b(15) );

    // Accumulate Jx for 16 particles into the v0-v3 vectors.
    v12_b  = q_b*ux_b;     // v12 = q ux
    v01_b  = v12_b*dy_b;   // v01 = q ux dy
    v00_b  = v12_b-v01_b;  // v00 = q ux (1-dy)
    v01_b += v12_b;        // v01 = q ux (1+dy)
    v12_b  = one+dz_b;     // v12 = 1+dz
    v02_b  = v00_b*v12_b;  // v02 = q ux (1-dy)(1+dz)
    v03_b  = v01_b*v12_b;  // v03 = q ux (1+dy)(1+dz)
    v12_b  = one-dz_b;     // v12 = 1-dz
    v00_b *= v12_b;        // v00 = q ux (1-dy)(1-dz)
    v01_b *= v12_b;        // v01 = q ux (1+dy)(1-dz)
    v00_b += v13_b;        // v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ]
    v01_b -= v13_b;        // v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ]
    v02_b -= v13_b;        // v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ]
    v03_b += v13_b;        // v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ]

    // Accumulate Jy for 16 particles into the v4-v7 vectors.
    v12_b  = q_b*uy_b;     // v12 = q uy
    v05_b  = v12_b*dz_b;   // v05 = q uy dz
    v04_b  = v12_b-v05_b;  // v04 = q uy (1-dz)
    v05_b += v12_b;        // v05 = q uy (1+dz)
    v12_b  = one+dx_b;     // v12 = 1+dx
    v06_b  = v04_b*v12_b;  // v06 = q uy (1-dz)(1+dx)
    v07_b  = v05_b*v12_b;  // v07 = q uy (1+dz)(1+dx)
    v12_b  = one-dx_b;     // v12 = 1-dx
    v04_b *= v12_b;        // v04 = q uy (1-dz)(1-dx)
    v05_b *= v12_b;        // v05 = q uy (1+dz)(1-dx)
    v04_b += v13_b;        // v04 = q uy [ (1-dz)(1-dx) + ux*uz/3 ]
    v05_b -= v13_b;        // v05 = q uy [ (1+dz)(1-dx) - ux*uz/3 ]
    v06_b -= v13_b;        // v06 = q uy [ (1-dz)(1+dx) - ux*uz/3 ]
    v07_b += v13_b;        // v07 = q uy [ (1+dz)(1+dx) + ux*uz/3 ]

    // Accumulate Jz for 16 particles into the v8-v11 vectors.
    v12_b  = q_b*uz_b;     // v12 = q uz
    v09_b  = v12_b*dx_b;   // v09 = q uz dx
    v08_b  = v12_b-v09_b;  // v08 = q uz (1-dx)
    v09_b += v12_b;        // v09 = q uz (1+dx)
    v12_b  = one+dy_b;     // v12 = 1+dy
    v10_b  = v08_b*v12_b;  // v10 = q uz (1-dx)(1+dy)
    v11_b  = v09_b*v12_b;  // v11 = q uz (1+dx)(1+dy)
    v12_b  = one-dy_b;     // v12 = 1-dy
    v08_b *= v12_b;        // v08 = q uz (1-dx)(1-dy)
    v09_b *= v12_b;        // v09 = q uz (1+dx)(1-dy)
    v08_b += v13_b;        // v08 = q uz [ (1-dx)(1-dy) + ux*uy/3 ]
    v09_b -= v13_b;        // v09 = q uz [ (1+dx)(1-dy) - ux*uy/3 ]
    v10_b -= v13_b;        // v10 = q uz [ (1-dx)(1+dy) - ux*uy/3 ]
    v11_b += v13_b;        // v11 = q uz [ (1+dx)(1+dy) + ux*uy/3 ]

    // Zero the v12-v15 vectors prior to transposing the data.
    v12_b = 0.0;
    v13_b = 0.0;
    v14_b = 0.0;
    v15_b = 0.0;

    // Transpose the data in vectors v0-v15 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v00_b, v01_b, v02_b, v03_b, v04_b, v05_b, v06_b, v07_b,
	       v08_b, v09_b, v10_b, v11_b, v12_b, v13_b, v14_b, v15_b );

    // Add the contributions to Jx, Jy and Jz from 16 particles into the
    // accumulator arrays for Jx, Jy and Jz.
    increment_16x1( vp00_b, v00_b );
    increment_16x1( vp01_b, v01_b );
    increment_16x1( vp02_b, v02_b );
    increment_16x1( vp03_b, v03_b );
    increment_16x1( vp04_b, v04_b );
    increment_16x1( vp05_b, v05_b );
    increment_16x1( vp06_b, v06_b );
    increment_16x1( vp07_b, v07_b );
    increment_16x1( vp08_b, v08_b );
    increment_16x1( vp09_b, v09_b );
    increment_16x1( vp10_b, v10_b );
    increment_16x1( vp11_b, v11_b );
    increment_16x1( vp12_b, v12_b );
    increment_16x1( vp13_b, v13_b );
    increment_16x1( vp14_b, v14_b );
    increment_16x1( vp15_b, v15_b );

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
    if ( outbnd_b(N) )                              /* Unlikely */      \
    {                                                                   \
      local_pm->dispx = ux_b(N);                                        \
      local_pm->dispy = uy_b(N);                                        \
      local_pm->dispz = uz_b(N);                                        \
      local_pm->i     = ( p - p0 ) + N + 16;                            \
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

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
#endif // Method 2

#if 0
//----------------------------------------------------------------------------//
// Method 3
//----------------------------------------------------------------------------//
// This approach does not process the particles in the same order as that of
// the reference implementation.  As a result, for certain configurations such
// as those that use Maxwellian reflux boundary conditions, reproducibility
// is lost because of random numbers being assigned to particles in a
// different order.
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

  float                * ALIGNED(16)  vp00_a;
  float                * ALIGNED(16)  vp01_a;
  float                * ALIGNED(16)  vp02_a;
  float                * ALIGNED(16)  vp03_a;
  float                * ALIGNED(16)  vp04_a;
  float                * ALIGNED(16)  vp05_a;
  float                * ALIGNED(16)  vp06_a;
  float                * ALIGNED(16)  vp07_a;
  float                * ALIGNED(16)  vp08_a;
  float                * ALIGNED(16)  vp09_a;
  float                * ALIGNED(16)  vp10_a;
  float                * ALIGNED(16)  vp11_a;
  float                * ALIGNED(16)  vp12_a;
  float                * ALIGNED(16)  vp13_a;
  float                * ALIGNED(16)  vp14_a;
  float                * ALIGNED(16)  vp15_a;

  float                * ALIGNED(16)  vp00_b;
  float                * ALIGNED(16)  vp01_b;
  float                * ALIGNED(16)  vp02_b;
  float                * ALIGNED(16)  vp03_b;
  float                * ALIGNED(16)  vp04_b;
  float                * ALIGNED(16)  vp05_b;
  float                * ALIGNED(16)  vp06_b;
  float                * ALIGNED(16)  vp07_b;
  float                * ALIGNED(16)  vp08_b;
  float                * ALIGNED(16)  vp09_b;
  float                * ALIGNED(16)  vp10_b;
  float                * ALIGNED(16)  vp11_b;
  float                * ALIGNED(16)  vp12_b;
  float                * ALIGNED(16)  vp13_b;
  float                * ALIGNED(16)  vp14_b;
  float                * ALIGNED(16)  vp15_b;

  const v16float qdt_2mc(args->qdt_2mc);
  const v16float cdt_dx(args->cdt_dx);
  const v16float cdt_dy(args->cdt_dy);
  const v16float cdt_dz(args->cdt_dz);
  const v16float qsp(args->qsp);
  const v16float one(1.);
  const v16float one_third(1./3.);
  const v16float two_fifteenths(2./15.);
  const v16float neg_one(-1.);

  const float _qsp = args->qsp;

  v16float dx_a, dy_a, dz_a, ux_a, uy_a, uz_a, q_a;
  v16float hax_a, hay_a, haz_a, cbx_a, cby_a, cbz_a;
  v16float v00_a, v01_a, v02_a, v03_a, v04_a, v05_a, v06_a, v07_a;
  v16float v08_a, v09_a, v10_a, v11_a, v12_a, v13_a, v14_a, v15_a;
  v16int   ii_a, outbnd_a;

  v16float dx_b, dy_b, dz_b, ux_b, uy_b, uz_b, q_b;
  v16float hax_b, hay_b, haz_b, cbx_b, cby_b, cbz_b;
  v16float v00_b, v01_b, v02_b, v03_b, v04_b, v05_b, v06_b, v07_b;
  v16float v08_b, v09_b, v10_b, v11_b, v12_b, v13_b, v14_b, v15_b;
  v16int   ii_b, outbnd_b;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 32, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 5;

  // Determine which movers are reserved for this pipeline.
  // Movers (16 bytes) should be reserved for pipelines in at least
  // multiples of 8 such that the set of particle movers reserved for
  // a pipeline is 128-byte aligned and a multiple of 128-byte in
  // size.  The host is guaranteed to get enough movers to process its
  // particles with this allocation.

  max_nm = args->max_nm - ( args->np&15 );

  if ( max_nm < 0 ) max_nm = 0;

  DISTRIBUTE( max_nm, 8, pipeline_rank, n_pipeline, itmp, max_nm );

  if ( pipeline_rank==n_pipeline ) max_nm = args->max_nm - itmp;

  pm   = args->pm + itmp;
  nm   = 0;
  itmp = 0;

  // Determine which accumulator array to use.
  // The host gets the first accumulator array.

  a0 += ( 1 + pipeline_rank ) *
        POW2_CEIL( (args->nx+2)*(args->ny+2)*(args->nz+2), 2 );

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=32 )
  {
    load_16x16_tr( &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
		   &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		   &p[16].dx, &p[18].dx, &p[20].dx, &p[22].dx,
		   &p[24].dx, &p[26].dx, &p[28].dx, &p[30].dx,
		   dx_a, dy_a, dz_a, ii_a, ux_a, uy_a, uz_a, q_a,
		   dx_b, dy_b, dz_b, ii_b, ux_b, uy_b, uz_b, q_b );

    //--------------------------------------------------------------------------
    // Interpolate fields for particles of block a.
    //--------------------------------------------------------------------------
    vp00_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 0) );
    vp01_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 1) );
    vp02_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 2) );
    vp03_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 3) );
    vp04_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 4) );
    vp05_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 5) );
    vp06_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 6) );
    vp07_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 7) );
    vp08_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 8) );
    vp09_a = ( float * ALIGNED(16) ) ( f0 + ii_a( 9) );
    vp10_a = ( float * ALIGNED(16) ) ( f0 + ii_a(10) );
    vp11_a = ( float * ALIGNED(16) ) ( f0 + ii_a(11) );
    vp12_a = ( float * ALIGNED(16) ) ( f0 + ii_a(12) );
    vp13_a = ( float * ALIGNED(16) ) ( f0 + ii_a(13) );
    vp14_a = ( float * ALIGNED(16) ) ( f0 + ii_a(14) );
    vp15_a = ( float * ALIGNED(16) ) ( f0 + ii_a(15) );

    //--------------------------------------------------------------------------
    // Load data for 16 particles for block a.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00_a, vp01_a, vp02_a, vp03_a,
		   vp04_a, vp05_a, vp06_a, vp07_a,
		   vp08_a, vp09_a, vp10_a, vp11_a,
		   vp12_a, vp13_a, vp14_a, vp15_a,
		   hax_a, v00_a, v01_a, v02_a, hay_a, v03_a, v04_a, v05_a,
		   haz_a, v06_a, v07_a, v08_a, cbx_a, v09_a, cby_a, v10_a );

    //--------------------------------------------------------------------------
    // Process particle block a.
    //--------------------------------------------------------------------------
    hax_a = qdt_2mc*
            fma( fma( v02_a, dy_a, v01_a ), dz_a, fma( v00_a, dy_a, hax_a ) );

    hay_a = qdt_2mc*
            fma( fma( v05_a, dz_a, v04_a ), dx_a, fma( v03_a, dz_a, hay_a ) );

    haz_a = qdt_2mc*
            fma( fma( v08_a, dx_a, v07_a ), dy_a, fma( v06_a, dx_a, haz_a ) );

    cbx_a = fma( v09_a, dx_a, cbx_a );

    cby_a = fma( v10_a, dy_a, cby_a );

    load_16x2_tr( vp00_a+16, vp01_a+16, vp02_a+16, vp03_a+16,
		  vp04_a+16, vp05_a+16, vp06_a+16, vp07_a+16,
		  vp08_a+16, vp09_a+16, vp10_a+16, vp11_a+16,
		  vp12_a+16, vp13_a+16, vp14_a+16, vp15_a+16,
		  cbz_a, v05_a );

    cbz_a = fma( v05_a, dz_a, cbz_a );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    ux_a += hax_a;
    uy_a += hay_a;
    uz_a += haz_a;
    v00_a  = qdt_2mc*rsqrt( one + fma( ux_a, ux_a,
                                       fma( uy_a, uy_a, uz_a*uz_a ) ) );
    v01_a  = fma( cbx_a, cbx_a, fma( cby_a, cby_a, cbz_a*cbz_a ) );
    v02_a  = (v00_a*v00_a)*v01_a;
    v03_a  = v00_a*fma( fma( two_fifteenths, v02_a, one_third ), v02_a, one );
    v04_a  = v03_a*rcp( fma( v03_a*v03_a, v01_a, one ) );
    v04_a += v04_a;
    v00_a  = fma( fms(  uy_a, cbz_a,  uz_a*cby_a ), v03_a, ux_a );
    v01_a  = fma( fms(  uz_a, cbx_a,  ux_a*cbz_a ), v03_a, uy_a );
    v02_a  = fma( fms(  ux_a, cby_a,  uy_a*cbx_a ), v03_a, uz_a );
    ux_a   = fma( fms( v01_a, cbz_a, v02_a*cby_a ), v04_a, ux_a );
    uy_a   = fma( fms( v02_a, cbx_a, v00_a*cbz_a ), v04_a, uy_a );
    uz_a   = fma( fms( v00_a, cby_a, v01_a*cbx_a ), v04_a, uz_a );
    ux_a  += hax_a;
    uy_a  += hay_a;
    uz_a  += haz_a;

    // Store ux, uy, uz in v6, v7, v8 so store_16x16_tr can be used below.
    v06_a  = ux_a;
    v07_a  = uy_a;
    v08_a  = uz_a;

    // Update the position of inbnd particles
    v00_a  = rsqrt( one + fma( ux_a, ux_a, fma( uy_a, uy_a, uz_a*uz_a ) ) );
    ux_a  *= cdt_dx;
    uy_a  *= cdt_dy;
    uz_a  *= cdt_dz;
    ux_a  *= v00_a;
    uy_a  *= v00_a;
    uz_a  *= v00_a;        // ux,uy,uz are normalized displ (relative to cell size)
    v00_a  =  dx_a + ux_a;
    v01_a  =  dy_a + uy_a;
    v02_a  =  dz_a + uz_a; // New particle midpoint
    v03_a  = v00_a + ux_a;
    v04_a  = v01_a + uy_a;
    v05_a  = v02_a + uz_a; // New particle position

    outbnd_a = ( v03_a > one ) | ( v03_a < neg_one ) |
               ( v04_a > one ) | ( v04_a < neg_one ) |
               ( v05_a > one ) | ( v05_a < neg_one );

    v03_a  = merge( outbnd_a, dx_a, v03_a ); // Do not update outbnd particles
    v04_a  = merge( outbnd_a, dy_a, v04_a );
    v05_a  = merge( outbnd_a, dz_a, v05_a );

    //--------------------------------------------------------------------------
    // Interpolate fields for particles of block b.
    //--------------------------------------------------------------------------
    vp00_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 0) );
    vp01_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 1) );
    vp02_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 2) );
    vp03_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 3) );
    vp04_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 4) );
    vp05_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 5) );
    vp06_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 6) );
    vp07_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 7) );
    vp08_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 8) );
    vp09_b = ( float * ALIGNED(16) ) ( f0 + ii_b( 9) );
    vp10_b = ( float * ALIGNED(16) ) ( f0 + ii_b(10) );
    vp11_b = ( float * ALIGNED(16) ) ( f0 + ii_b(11) );
    vp12_b = ( float * ALIGNED(16) ) ( f0 + ii_b(12) );
    vp13_b = ( float * ALIGNED(16) ) ( f0 + ii_b(13) );
    vp14_b = ( float * ALIGNED(16) ) ( f0 + ii_b(14) );
    vp15_b = ( float * ALIGNED(16) ) ( f0 + ii_b(15) );

    //--------------------------------------------------------------------------
    // Load data for 16 particles for block b.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00_b, vp01_b, vp02_b, vp03_b,
		   vp04_b, vp05_b, vp06_b, vp07_b,
		   vp08_b, vp09_b, vp10_b, vp11_b,
		   vp12_b, vp13_b, vp14_b, vp15_b,
		   hax_b, v00_b, v01_b, v02_b, hay_b, v03_b, v04_b, v05_b,
		   haz_b, v06_b, v07_b, v08_b, cbx_b, v09_b, cby_b, v10_b );

    //--------------------------------------------------------------------------
    // Process particle block b.
    //--------------------------------------------------------------------------
    hax_b = qdt_2mc*
            fma( fma( v02_b, dy_b, v01_b ), dz_b, fma( v00_b, dy_b, hax_b ) );

    hay_b = qdt_2mc*
            fma( fma( v05_b, dz_b, v04_b ), dx_b, fma( v03_b, dz_b, hay_b ) );

    haz_b = qdt_2mc*
            fma( fma( v08_b, dx_b, v07_b ), dy_b, fma( v06_b, dx_b, haz_b ) );

    cbx_b = fma( v09_b, dx_b, cbx_b );

    cby_b = fma( v10_b, dy_b, cby_b );

    load_16x2_tr( vp00_b+16, vp01_b+16, vp02_b+16, vp03_b+16,
		  vp04_b+16, vp05_b+16, vp06_b+16, vp07_b+16,
		  vp08_b+16, vp09_b+16, vp10_b+16, vp11_b+16,
		  vp12_b+16, vp13_b+16, vp14_b+16, vp15_b+16,
		  cbz_b, v05_b );

    cbz_b = fma( v05_b, dz_b, cbz_b );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    ux_b  += hax_b;
    uy_b  += hay_b;
    uz_b  += haz_b;
    v00_b  = qdt_2mc*rsqrt( one + fma( ux_b, ux_b,
                                       fma( uy_b, uy_b, uz_b*uz_b ) ) );
    v01_b  = fma( cbx_b, cbx_b, fma( cby_b, cby_b, cbz_b*cbz_b ) );
    v02_b  = (v00_b*v00_b)*v01_b;
    v03_b  = v00_b*fma( fma( two_fifteenths, v02_b, one_third ), v02_b, one );
    v04_b  = v03_b*rcp( fma( v03_b*v03_b, v01_b, one ) );
    v04_b += v04_b;
    v00_b  = fma( fms(  uy_b, cbz_b,  uz_b*cby_b ), v03_b, ux_b );
    v01_b  = fma( fms(  uz_b, cbx_b,  ux_b*cbz_b ), v03_b, uy_b );
    v02_b  = fma( fms(  ux_b, cby_b,  uy_b*cbx_b ), v03_b, uz_b );
    ux_b   = fma( fms( v01_b, cbz_b, v02_b*cby_b ), v04_b, ux_b );
    uy_b   = fma( fms( v02_b, cbx_b, v00_b*cbz_b ), v04_b, uy_b );
    uz_b   = fma( fms( v00_b, cby_b, v01_b*cbx_b ), v04_b, uz_b );
    ux_b  += hax_b;
    uy_b  += hay_b;
    uz_b  += haz_b;

    // Store ux, uy, uz in v6, v7, v8 so store_16x16_tr can be used below.
    v06_b  = ux_b;
    v07_b  = uy_b;
    v08_b  = uz_b;

    // Update the position of inbnd particles
    v00_b  = rsqrt( one + fma( ux_b, ux_b, fma( uy_b, uy_b, uz_b*uz_b ) ) );
    ux_b  *= cdt_dx;
    uy_b  *= cdt_dy;
    uz_b  *= cdt_dz;
    ux_b  *= v00_b;
    uy_b  *= v00_b;
    uz_b  *= v00_b;        // ux,uy,uz are normalized displ (relative to cell size)
    v00_b  =  dx_b + ux_b;
    v01_b  =  dy_b + uy_b;
    v02_b  =  dz_b + uz_b; // New particle midpoint
    v03_b  = v00_b + ux_b;
    v04_b  = v01_b + uy_b;
    v05_b  = v02_b + uz_b; // New particle position

    outbnd_b = ( v03_b > one ) | ( v03_b < neg_one ) |
               ( v04_b > one ) | ( v04_b < neg_one ) |
               ( v05_b > one ) | ( v05_b < neg_one );

    v03_b  = merge( outbnd_b, dx_b, v03_b ); // Do not update outbnd particles
    v04_b  = merge( outbnd_b, dy_b, v04_b );
    v05_b  = merge( outbnd_b, dz_b, v05_b );

    //--------------------------------------------------------------------------
    // Store results for 2 blocks of particles.
    //--------------------------------------------------------------------------
    store_16x16_tr( v03_a, v04_a, v05_a, ii_a, v06_a, v07_a, v08_a, q_a,
		    v03_b, v04_b, v05_b, ii_b, v06_b, v07_b, v08_b, q_b,
		    &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
		    &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		    &p[16].dx, &p[18].dx, &p[20].dx, &p[22].dx,
		    &p[24].dx, &p[26].dx, &p[28].dx, &p[30].dx );

    //--------------------------------------------------------------------------
    // Process particle block a.
    //--------------------------------------------------------------------------

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q_a  = czero( outbnd_a, q_a*qsp );   // Do not accumulate outbnd particles
    dx_a = v00_a;                        // Streak midpoint (valid for inbnd only)
    dy_a = v01_a;
    dz_a = v02_a;

    v13_a = q_a*ux_a*uy_a*uz_a*one_third;     // Charge conservation correction

    vp00_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 0) ); // Accumulator pointers
    vp01_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 1) );
    vp02_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 2) );
    vp03_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 3) );
    vp04_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 4) );
    vp05_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 5) );
    vp06_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 6) );
    vp07_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 7) );
    vp08_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 8) );
    vp09_a = ( float * ALIGNED(16) ) ( a0 + ii_a( 9) );
    vp10_a = ( float * ALIGNED(16) ) ( a0 + ii_a(10) );
    vp11_a = ( float * ALIGNED(16) ) ( a0 + ii_a(11) );
    vp12_a = ( float * ALIGNED(16) ) ( a0 + ii_a(12) );
    vp13_a = ( float * ALIGNED(16) ) ( a0 + ii_a(13) );
    vp14_a = ( float * ALIGNED(16) ) ( a0 + ii_a(14) );
    vp15_a = ( float * ALIGNED(16) ) ( a0 + ii_a(15) );

    // Accumulate Jx for 16 particles into the v0-v3 vectors.
    v12_a  = q_a*ux_a;     // v12 = q ux
    v01_a  = v12_a*dy_a;   // v01 = q ux dy
    v00_a  = v12_a-v01_a;  // v00 = q ux (1-dy)
    v01_a += v12_a;        // v01 = q ux (1+dy)
    v12_a  = one+dz_a;     // v12 = 1+dz
    v02_a  = v00_a*v12_a;  // v02 = q ux (1-dy)(1+dz)
    v03_a  = v01_a*v12_a;  // v03 = q ux (1+dy)(1+dz)
    v12_a  = one-dz_a;     // v12 = 1-dz
    v00_a *= v12_a;        // v00 = q ux (1-dy)(1-dz)
    v01_a *= v12_a;        // v01 = q ux (1+dy)(1-dz)
    v00_a += v13_a;        // v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ]
    v01_a -= v13_a;        // v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ]
    v02_a -= v13_a;        // v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ]
    v03_a += v13_a;        // v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ]

    // Accumulate Jy for 16 particles into the v4-v7 vectors.
    v12_a  = q_a*uy_a;     // v12 = q uy
    v05_a  = v12_a*dz_a;   // v05 = q uy dz
    v04_a  = v12_a-v05_a;  // v04 = q uy (1-dz)
    v05_a += v12_a;        // v05 = q uy (1+dz)
    v12_a  = one+dx_a;     // v12 = 1+dx
    v06_a  = v04_a*v12_a;  // v06 = q uy (1-dz)(1+dx)
    v07_a  = v05_a*v12_a;  // v07 = q uy (1+dz)(1+dx)
    v12_a  = one-dx_a;     // v12 = 1-dx
    v04_a *= v12_a;        // v04 = q uy (1-dz)(1-dx)
    v05_a *= v12_a;        // v05 = q uy (1+dz)(1-dx)
    v04_a += v13_a;        // v04 = q uy [ (1-dz)(1-dx) + ux*uz/3 ]
    v05_a -= v13_a;        // v05 = q uy [ (1+dz)(1-dx) - ux*uz/3 ]
    v06_a -= v13_a;        // v06 = q uy [ (1-dz)(1+dx) - ux*uz/3 ]
    v07_a += v13_a;        // v07 = q uy [ (1+dz)(1+dx) + ux*uz/3 ]

    // Accumulate Jz for 16 particles into the v8-v11 vectors.
    v12_a  = q_a*uz_a;     // v12 = q uz
    v09_a  = v12_a*dx_a;   // v09 = q uz dx
    v08_a  = v12_a-v09_a;  // v08 = q uz (1-dx)
    v09_a += v12_a;        // v09 = q uz (1+dx)
    v12_a  = one+dy_a;     // v12 = 1+dy
    v10_a  = v08_a*v12_a;  // v10 = q uz (1-dx)(1+dy)
    v11_a  = v09_a*v12_a;  // v11 = q uz (1+dx)(1+dy)
    v12_a  = one-dy_a;     // v12 = 1-dy
    v08_a *= v12_a;        // v08 = q uz (1-dx)(1-dy)
    v09_a *= v12_a;        // v09 = q uz (1+dx)(1-dy)
    v08_a += v13_a;        // v08 = q uz [ (1-dx)(1-dy) + ux*uy/3 ]
    v09_a -= v13_a;        // v09 = q uz [ (1+dx)(1-dy) - ux*uy/3 ]
    v10_a -= v13_a;        // v10 = q uz [ (1-dx)(1+dy) - ux*uy/3 ]
    v11_a += v13_a;        // v11 = q uz [ (1+dx)(1+dy) + ux*uy/3 ]

    // Zero the v12-v15 vectors prior to transposing the data.
    v12_a = 0.0;
    v13_a = 0.0;
    v14_a = 0.0;
    v15_a = 0.0;

    // Transpose the data in vectors v0-v15 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v00_a, v01_a, v02_a, v03_a, v04_a, v05_a, v06_a, v07_a,
	       v08_a, v09_a, v10_a, v11_a, v12_a, v13_a, v14_a, v15_a );

    // Add the contributions to Jx, Jy and Jz from 16 particles into the
    // accumulator arrays for Jx, Jy and Jz.
    increment_16x1( vp00_a, v00_a );
    increment_16x1( vp01_a, v01_a );
    increment_16x1( vp02_a, v02_a );
    increment_16x1( vp03_a, v03_a );
    increment_16x1( vp04_a, v04_a );
    increment_16x1( vp05_a, v05_a );
    increment_16x1( vp06_a, v06_a );
    increment_16x1( vp07_a, v07_a );
    increment_16x1( vp08_a, v08_a );
    increment_16x1( vp09_a, v09_a );
    increment_16x1( vp10_a, v10_a );
    increment_16x1( vp11_a, v11_a );
    increment_16x1( vp12_a, v12_a );
    increment_16x1( vp13_a, v13_a );
    increment_16x1( vp14_a, v14_a );
    increment_16x1( vp15_a, v15_a );

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
    if ( outbnd_a(N) )                              /* Unlikely */      \
    {                                                                   \
      local_pm->dispx = ux_a(N);                                        \
      local_pm->dispy = uy_a(N);                                        \
      local_pm->dispz = uz_a(N);                                        \
      local_pm->i     = ( p - p0 ) + 2*N;                               \
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

    //--------------------------------------------------------------------------
    // Process particle block b.
    //--------------------------------------------------------------------------

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q_b  = czero( outbnd_b, q_b*qsp );   // Do not accumulate outbnd particles
    dx_b = v00_b;                        // Streak midpoint (valid for inbnd only)
    dy_b = v01_b;
    dz_b = v02_b;

    v13_b = q_b*ux_b*uy_b*uz_b*one_third;     // Charge conservation correction

    vp00_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 0) ); // Accumulator pointers
    vp01_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 1) );
    vp02_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 2) );
    vp03_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 3) );
    vp04_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 4) );
    vp05_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 5) );
    vp06_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 6) );
    vp07_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 7) );
    vp08_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 8) );
    vp09_b = ( float * ALIGNED(16) ) ( a0 + ii_b( 9) );
    vp10_b = ( float * ALIGNED(16) ) ( a0 + ii_b(10) );
    vp11_b = ( float * ALIGNED(16) ) ( a0 + ii_b(11) );
    vp12_b = ( float * ALIGNED(16) ) ( a0 + ii_b(12) );
    vp13_b = ( float * ALIGNED(16) ) ( a0 + ii_b(13) );
    vp14_b = ( float * ALIGNED(16) ) ( a0 + ii_b(14) );
    vp15_b = ( float * ALIGNED(16) ) ( a0 + ii_b(15) );

    // Accumulate Jx for 16 particles into the v0-v3 vectors.
    v12_b  = q_b*ux_b;     // v12 = q ux
    v01_b  = v12_b*dy_b;   // v01 = q ux dy
    v00_b  = v12_b-v01_b;  // v00 = q ux (1-dy)
    v01_b += v12_b;        // v01 = q ux (1+dy)
    v12_b  = one+dz_b;     // v12 = 1+dz
    v02_b  = v00_b*v12_b;  // v02 = q ux (1-dy)(1+dz)
    v03_b  = v01_b*v12_b;  // v03 = q ux (1+dy)(1+dz)
    v12_b  = one-dz_b;     // v12 = 1-dz
    v00_b *= v12_b;        // v00 = q ux (1-dy)(1-dz)
    v01_b *= v12_b;        // v01 = q ux (1+dy)(1-dz)
    v00_b += v13_b;        // v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ]
    v01_b -= v13_b;        // v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ]
    v02_b -= v13_b;        // v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ]
    v03_b += v13_b;        // v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ]

    // Accumulate Jy for 16 particles into the v4-v7 vectors.
    v12_b  = q_b*uy_b;     // v12 = q uy
    v05_b  = v12_b*dz_b;   // v05 = q uy dz
    v04_b  = v12_b-v05_b;  // v04 = q uy (1-dz)
    v05_b += v12_b;        // v05 = q uy (1+dz)
    v12_b  = one+dx_b;     // v12 = 1+dx
    v06_b  = v04_b*v12_b;  // v06 = q uy (1-dz)(1+dx)
    v07_b  = v05_b*v12_b;  // v07 = q uy (1+dz)(1+dx)
    v12_b  = one-dx_b;     // v12 = 1-dx
    v04_b *= v12_b;        // v04 = q uy (1-dz)(1-dx)
    v05_b *= v12_b;        // v05 = q uy (1+dz)(1-dx)
    v04_b += v13_b;        // v04 = q uy [ (1-dz)(1-dx) + ux*uz/3 ]
    v05_b -= v13_b;        // v05 = q uy [ (1+dz)(1-dx) - ux*uz/3 ]
    v06_b -= v13_b;        // v06 = q uy [ (1-dz)(1+dx) - ux*uz/3 ]
    v07_b += v13_b;        // v07 = q uy [ (1+dz)(1+dx) + ux*uz/3 ]

    // Accumulate Jz for 16 particles into the v8-v11 vectors.
    v12_b  = q_b*uz_b;     // v12 = q uz
    v09_b  = v12_b*dx_b;   // v09 = q uz dx
    v08_b  = v12_b-v09_b;  // v08 = q uz (1-dx)
    v09_b += v12_b;        // v09 = q uz (1+dx)
    v12_b  = one+dy_b;     // v12 = 1+dy
    v10_b  = v08_b*v12_b;  // v10 = q uz (1-dx)(1+dy)
    v11_b  = v09_b*v12_b;  // v11 = q uz (1+dx)(1+dy)
    v12_b  = one-dy_b;     // v12 = 1-dy
    v08_b *= v12_b;        // v08 = q uz (1-dx)(1-dy)
    v09_b *= v12_b;        // v09 = q uz (1+dx)(1-dy)
    v08_b += v13_b;        // v08 = q uz [ (1-dx)(1-dy) + ux*uy/3 ]
    v09_b -= v13_b;        // v09 = q uz [ (1+dx)(1-dy) - ux*uy/3 ]
    v10_b -= v13_b;        // v10 = q uz [ (1-dx)(1+dy) - ux*uy/3 ]
    v11_b += v13_b;        // v11 = q uz [ (1+dx)(1+dy) + ux*uy/3 ]

    // Zero the v12-v15 vectors prior to transposing the data.
    v12_b = 0.0;
    v13_b = 0.0;
    v14_b = 0.0;
    v15_b = 0.0;

    // Transpose the data in vectors v0-v15 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v00_b, v01_b, v02_b, v03_b, v04_b, v05_b, v06_b, v07_b,
	       v08_b, v09_b, v10_b, v11_b, v12_b, v13_b, v14_b, v15_b );

    // Add the contributions to Jx, Jy and Jz from 16 particles into the
    // accumulator arrays for Jx, Jy and Jz.
    increment_16x1( vp00_b, v00_b );
    increment_16x1( vp01_b, v01_b );
    increment_16x1( vp02_b, v02_b );
    increment_16x1( vp03_b, v03_b );
    increment_16x1( vp04_b, v04_b );
    increment_16x1( vp05_b, v05_b );
    increment_16x1( vp06_b, v06_b );
    increment_16x1( vp07_b, v07_b );
    increment_16x1( vp08_b, v08_b );
    increment_16x1( vp09_b, v09_b );
    increment_16x1( vp10_b, v10_b );
    increment_16x1( vp11_b, v11_b );
    increment_16x1( vp12_b, v12_b );
    increment_16x1( vp13_b, v13_b );
    increment_16x1( vp14_b, v14_b );
    increment_16x1( vp15_b, v15_b );

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
    if ( outbnd_b(N) )                              /* Unlikely */      \
    {                                                                   \
      local_pm->dispx = ux_b(N);                                        \
      local_pm->dispy = uy_b(N);                                        \
      local_pm->dispz = uz_b(N);                                        \
      local_pm->i     = ( p - p0 ) + 2*N + 1;                           \
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

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
#endif // Method 3

// #if 0
//----------------------------------------------------------------------------//
// Method 4
//----------------------------------------------------------------------------//
// This method processes 16 particles at a time instead of 32.
//----------------------------------------------------------------------------//
// This method processes the particles in the same order as the reference
// implementation and gives good reproducibility. This is achieved using
// modified load_16x16_tr_p and store_16x16_tr_p functions which load or
// store the particle data in the correct order in a single step instead
// of using two steps.
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

  if ( pipeline_rank==n_pipeline ) max_nm = args->max_nm - itmp;

  pm   = args->pm + itmp;
  nm   = 0;
  itmp = 0;

  // Determine which accumulator array to use.
  // The host gets the first accumulator array.

  a0 += ( 1 + pipeline_rank ) *
        POW2_CEIL( (args->nx+2)*(args->ny+2)*(args->nz+2), 2 );

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=16 )
  {
    load_16x8_tr_p( &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
                    &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		    dx, dy, dz, ii, ux, uy, uz, q );

    //--------------------------------------------------------------------------
    // Interpolate fields for particles of block a.
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
    // Load data for 16 particles for block a.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00, vp01, vp02, vp03,
		   vp04, vp05, vp06, vp07,
		   vp08, vp09, vp10, vp11,
		   vp12, vp13, vp14, vp15,
		   hax, v00, v01, v02, hay, v03, v04, v05,
		   haz, v06, v07, v08, cbx, v09, cby, v10 );

    //--------------------------------------------------------------------------
    // Process particles for block a.
    //--------------------------------------------------------------------------
    hax = qdt_2mc*
          fma( fma( v02, dy, v01 ), dz, fma( v00, dy, hax ) );

    hay = qdt_2mc*
          fma( fma( v05, dz, v04 ), dx, fma( v03, dz, hay ) );

    haz = qdt_2mc*
          fma( fma( v08, dx, v07 ), dy, fma( v06, dx, haz ) );

    cbx = fma( v09, dx, cbx );

    cby = fma( v10, dy, cby );

    load_16x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
		  vp04+16, vp05+16, vp06+16, vp07+16,
		  vp08+16, vp09+16, vp10+16, vp11+16,
		  vp12+16, vp13+16, vp14+16, vp15+16,
		  cbz, v05 );

    cbz = fma( v05, dz, cbz );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    ux += hax;
    uy += hay;
    uz += haz;
    v00  = qdt_2mc*rsqrt( one + fma( ux, ux,
				     fma( uy, uy, uz*uz ) ) );
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

    // Store ux, uy, uz in v6, v7, v8 so store_16x8_tr_p can be used below.
    v06  = ux;
    v07  = uy;
    v08  = uz;

    // Update the position of inbnd particles
    v00  = rsqrt( one + fma( ux, ux, fma( uy, uy, uz*uz ) ) );
    ux  *= cdt_dx;
    uy  *= cdt_dy;
    uz  *= cdt_dz;
    ux  *= v00;
    uy  *= v00;
    uz  *= v00;      // ux,uy,uz are normalized displ (relative to cell size)
    v00  =  dx + ux;
    v01  =  dy + uy;
    v02  =  dz + uz; // New particle midpoint
    v03  = v00 + ux;
    v04  = v01 + uy;
    v05  = v02 + uz; // New particle position

    outbnd = ( v03 > one ) | ( v03 < neg_one ) |
             ( v04 > one ) | ( v04 < neg_one ) |
             ( v05 > one ) | ( v05 < neg_one );

    v03  = merge( outbnd, dx, v03 ); // Do not update outbnd particles
    v04  = merge( outbnd, dy, v04 );
    v05  = merge( outbnd, dz, v05 );

    //--------------------------------------------------------------------------
    // Store the particle results.
    //--------------------------------------------------------------------------
    store_16x8_tr_p( v03, v04, v05, ii, v06, v07, v08, q,
                     &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
                     &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx );

    //--------------------------------------------------------------------------
    // Process particle block a.
    //--------------------------------------------------------------------------

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q  = czero( outbnd, q*qsp );   // Do not accumulate outbnd particles
    dx = v00;                      // Streak midpoint (valid for inbnd only)
    dy = v01;
    dz = v02;

    v13 = q*ux*uy*uz*one_third;    // Charge conservation correction

    vp00 = ( float * ALIGNED(16) ) ( a0 + ii( 0) ); // Accumulator pointers
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

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
    if ( outbnd(N) )                              /* Unlikely */        \
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

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
// #endif // Method 4
