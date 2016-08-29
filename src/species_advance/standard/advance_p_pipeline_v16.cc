using namespace v16;

// #if 0
// Method 1.
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

    // Interpolate fields.
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

#   define ACCUMULATE_JX(X,Y,Z)                                        \
    v12_a  = q_a*u##X_a;   /* v12 = q ux                            */ \
    v01_a  = v12_a*d##Y_a; /* v01 = q ux dy                         */ \
    v00_a  = v12_a-v01_a;  /* v00 = q ux (1-dy)                     */ \
    v01_a += v12_a;        /* v01 = q ux (1+dy)                     */ \
    v12_a  = one+d##Z_a;   /* v12 = 1+dz                            */ \
    v02_a  = v00_a*v12_a;  /* v02 = q ux (1-dy)(1+dz)               */ \
    v03_a  = v01_a*v12_a;  /* v03 = q ux (1+dy)(1+dz)               */ \
    v12_a  = one-d##Z_a;   /* v12 = 1-dz                            */ \
    v00_a *= v12_a;        /* v00 = q ux (1-dy)(1-dz)               */ \
    v01_a *= v12_a;        /* v01 = q ux (1+dy)(1-dz)               */ \
    v00_a += v13_a;        /* v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */ \
    v01_a -= v13_a;        /* v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */ \
    v02_a -= v13_a;        /* v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */ \
    v03_a += v13_a;        /* v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JY(X,Y,Z)                                        \
    v12_a  = q_a*u##X_a;   /* v12 = q ux                            */ \
    v05_a  = v12_a*d##Y_a; /* v05 = q ux dy                         */ \
    v04_a  = v12_a-v05_a;  /* v04 = q ux (1-dy)                     */ \
    v05_a += v12_a;        /* v05 = q ux (1+dy)                     */ \
    v12_a  = one+d##Z_a;   /* v12 = 1+dz                            */ \
    v06_a  = v04_a*v12_a;  /* v06 = q ux (1-dy)(1+dz)               */ \
    v07_a  = v05_a*v12_a;  /* v07 = q ux (1+dy)(1+dz)               */ \
    v12_a  = one-d##Z_a;   /* v12 = 1-dz                            */ \
    v04_a *= v12_a;        /* v04 = q ux (1-dy)(1-dz)               */ \
    v05_a *= v12_a;        /* v05 = q ux (1+dy)(1-dz)               */ \
    v04_a += v13_a;        /* v04 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */ \
    v05_a -= v13_a;        /* v05 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */ \
    v06_a -= v13_a;        /* v06 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */ \
    v07_a += v13_a;        /* v07 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JZ(X,Y,Z)                                        \
    v12_a  = q_a*u##X_a;   /* v12 = q ux                            */ \
    v09_a  = v12_a*d##Y_a; /* v09 = q ux dy                         */ \
    v08_a  = v12_a-v09_a;  /* v08 = q ux (1-dy)                     */ \
    v09_a += v12_a;        /* v09 = q ux (1+dy)                     */ \
    v12_a  = one+d##Z_a;   /* v12 = 1+dz                            */ \
    v10_a  = v08_a*v12_a;  /* v10 = q ux (1-dy)(1+dz)               */ \
    v11_a  = v09_a*v12_a;  /* v11 = q ux (1+dy)(1+dz)               */ \
    v12_a  = one-d##Z_a;   /* v12 = 1-dz                            */ \
    v08_a *= v12_a;        /* v08 = q ux (1-dy)(1-dz)               */ \
    v09_a *= v12_a;        /* v09 = q ux (1+dy)(1-dz)               */ \
    v08_a += v13_a;        /* v08 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */ \
    v09_a -= v13_a;        /* v09 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */ \
    v10_a -= v13_a;        /* v10 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */ \
    v11_a += v13_a;        /* v11 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

    // Accumulate Jx for 16 particles into the v0-v3 vectors.
    // ACCUMULATE_JX( x, y, z );
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
    // ACCUMULATE_JY( y, z, x );
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
    // ACCUMULATE_JZ( z, x, y );
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

#   undef ACCUMULATE_JX
#   undef ACCUMULATE_JY
#   undef ACCUMULATE_JZ

    // Update position and accumulate outbnd

// #   define MOVE_OUTBND(N)                                               \
//     if ( outbnd_a(N) )                              /* Unlikely */      \
//     {                                                                   \
//       local_pm->dispx = ux_a(N);                                        \
//       local_pm->dispy = uy_a(N);                                        \
//       local_pm->dispz = uz_a(N);                                        \
//       local_pm->i     = ( p - p0 ) + N;                                 \
//       if ( move_p( p0, local_pm, a0, g, _qsp ) )    /* Unlikely */      \
//       {                                                                 \
//         if ( nm < max_nm )                                              \
//         {                                                               \
// 	  v4::copy_4x1( &pm[nm++], local_pm ); 	                        \
// 	}                                                               \
//         else                                        /* Unlikely */      \
// 	{                                                               \
// 	  itmp++;                                                       \
// 	}                                                               \
//       }                                                                 \
//     }

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

#   define ACCUMULATE_JX(X,Y,Z)                                        \
    v12_b  = q_b*u##X_b;   /* v12 = q ux                            */ \
    v01_b  = v12_b*d##Y_b; /* v01 = q ux dy                         */ \
    v00_b  = v12_b-v01_b;  /* v00 = q ux (1-dy)                     */ \
    v01_b += v12_b;        /* v01 = q ux (1+dy)                     */ \
    v12_b  = one+d##Z_b;   /* v12 = 1+dz                            */ \
    v02_b  = v00_b*v12_b;  /* v02 = q ux (1-dy)(1+dz)               */ \
    v03_b  = v01_b*v12_b;  /* v03 = q ux (1+dy)(1+dz)               */ \
    v12_b  = one-d##Z_b;   /* v12 = 1-dz                            */ \
    v00_b *= v12_b;        /* v00 = q ux (1-dy)(1-dz)               */ \
    v01_b *= v12_b;        /* v01 = q ux (1+dy)(1-dz)               */ \
    v00_b += v13_b;        /* v00 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */ \
    v01_b -= v13_b;        /* v01 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */ \
    v02_b -= v13_b;        /* v02 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */ \
    v03_b += v13_b;        /* v03 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JY(X,Y,Z)                                        \
    v12_b  = q_b*u##X_b;   /* v12 = q ux                            */ \
    v05_b  = v12_b*d##Y_b; /* v05 = q ux dy                         */ \
    v04_b  = v12_b-v05_b;  /* v04 = q ux (1-dy)                     */ \
    v05_b += v12_b;        /* v05 = q ux (1+dy)                     */ \
    v12_b  = one+d##Z_b;   /* v12 = 1+dz                            */ \
    v06_b  = v04_b*v12_b;  /* v06 = q ux (1-dy)(1+dz)               */ \
    v07_b  = v05_b*v12_b;  /* v07 = q ux (1+dy)(1+dz)               */ \
    v12_b  = one-d##Z_b;   /* v12 = 1-dz                            */ \
    v04_b *= v12_b;        /* v04 = q ux (1-dy)(1-dz)               */ \
    v05_b *= v12_b;        /* v05 = q ux (1+dy)(1-dz)               */ \
    v04_b += v13_b;        /* v04 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */ \
    v05_b -= v13_b;        /* v05 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */ \
    v06_b -= v13_b;        /* v06 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */ \
    v07_b += v13_b;        /* v07 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JZ(X,Y,Z)                                        \
    v12_b  = q_b*u##X_b;   /* v12 = q ux                            */ \
    v09_b  = v12_b*d##Y_b; /* v09 = q ux dy                         */ \
    v08_b  = v12_b-v09_b;  /* v08 = q ux (1-dy)                     */ \
    v09_b += v12_b;        /* v09 = q ux (1+dy)                     */ \
    v12_b  = one+d##Z_b;   /* v12 = 1+dz                            */ \
    v10_b  = v08_b*v12_b;  /* v10 = q ux (1-dy)(1+dz)               */ \
    v11_b  = v09_b*v12_b;  /* v11 = q ux (1+dy)(1+dz)               */ \
    v12_b  = one-d##Z_b;   /* v12 = 1-dz                            */ \
    v08_b *= v12_b;        /* v08 = q ux (1-dy)(1-dz)               */ \
    v09_b *= v12_b;        /* v09 = q ux (1+dy)(1-dz)               */ \
    v08_b += v13_b;        /* v08 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */ \
    v09_b -= v13_b;        /* v09 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */ \
    v10_b -= v13_b;        /* v10 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */ \
    v11_b += v13_b;        /* v11 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

    // Accumulate Jx for 16 particles into the v0-v3 vectors.
    // ACCUMULATE_JX( x, y, z );
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
    // ACCUMULATE_JY( y, z, x );
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
    // ACCUMULATE_JZ( z, x, y );
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

#   undef ACCUMULATE_JX
#   undef ACCUMULATE_JY
#   undef ACCUMULATE_JZ

    // Update position and accumulate outbnd

// #   define MOVE_OUTBND(N)                                               \
//     if ( outbnd_b(N) )                              /* Unlikely */      \
//     {                                                                   \
//       local_pm->dispx = ux_b(N);                                        \
//       local_pm->dispy = uy_b(N);                                        \
//       local_pm->dispz = uz_b(N);                                        \
//       local_pm->i     = ( p - p0 ) + N;                                 \
//       if ( move_p( p0, local_pm, a0, g, _qsp ) )    /* Unlikely */      \
//       {                                                                 \
//         if ( nm < max_nm )                                              \
//         {                                                               \
// 	  v4::copy_4x1( &pm[nm++], local_pm ); 	                        \
// 	}                                                               \
//         else                                        /* Unlikely */      \
// 	{                                                               \
// 	  itmp++;                                                       \
// 	}                                                               \
//       }                                                                 \
//     }

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
// #endif // Method 1


















#if 0
// Method 1.
void
advance_p_pipeline_v8( advance_p_pipeline_args_t * args,
		       int pipeline_rank,
		       int n_pipeline )
{
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
  float                * ALIGNED(16)  vp4;
  float                * ALIGNED(16)  vp5;
  float                * ALIGNED(16)  vp6;
  float                * ALIGNED(16)  vp7;

  const v8float qdt_2mc(args->qdt_2mc);
  const v8float cdt_dx(args->cdt_dx);
  const v8float cdt_dy(args->cdt_dy);
  const v8float cdt_dz(args->cdt_dz);
  const v8float qsp(args->qsp);
  const v8float one(1.);
  const v8float one_third(1./3.);
  const v8float two_fifteenths(2./15.);
  const v8float neg_one(-1.);

  const float _qsp = args->qsp;

  v8float dx, dy, dz, ux, uy, uz, q;
  v8float hax, hay, haz, cbx, cby, cbz;
  v8float v0, v1, v2, v3, v4, v5, v6, v7, v8;
  v8int   ii, outbnd;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 3;

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

  // Process the particle quads for this pipeline.

  for( ; nq; nq--, p+=8 )
  {
    load_8x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		 &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
		 dx, dy, dz, ii );

    // Interpolate fields.
    vp0 = ( float * ALIGNED(16) ) ( f0 + ii(0) );
    vp1 = ( float * ALIGNED(16) ) ( f0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( f0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( f0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( f0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( f0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( f0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( f0 + ii(7) );

    load_8x4_tr( vp0, vp1, vp2, vp3,
		 vp4, vp5, vp6, vp7,
		 hax, v0, v1, v2 );
    hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );

    load_8x4_tr( vp0+4, vp1+4, vp2+4, vp3+4,
		 vp4+4, vp5+4, vp6+4, vp7+4,
		 hay, v3, v4, v5 );
    hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );

    load_8x4_tr( vp0+8, vp1+8, vp2+8, vp3+8,
		 vp4+8, vp5+8, vp6+8, vp7+8,
		 haz, v0, v1, v2 );
    haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );

    load_8x4_tr( vp0+12, vp1+12, vp2+12, vp3+12,
		 vp4+12, vp5+12, vp6+12, vp7+12,
		 cbx, v3, cby, v4 );
    cbx = fma( v3, dx, cbx );
    cby = fma( v4, dy, cby );

    load_8x2_tr( vp0+16, vp1+16, vp2+16, vp3+16,
		 vp4+16, vp5+16, vp6+16, vp7+16,
		 cbz, v5 );
    cbz = fma( v5, dz, cbz );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    load_8x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		 &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux,
		 ux, uy, uz, q );

    ux += hax;
    uy += hay;
    uz += haz;
    v0  = qdt_2mc*rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
    v1  = fma( cbx,cbx, fma( cby,cby, cbz*cbz ) );
    v2  = (v0*v0)*v1;
    v3  = v0*fma( fma( two_fifteenths, v2, one_third ), v2, one );
    v4  = v3*rcp(fma( v3*v3, v1, one ));
    v4 += v4;
    v0  = fma( fms( uy, cbz, uz*cby ), v3, ux );
    v1  = fma( fms( uz, cbx, ux*cbz ), v3, uy );
    v2  = fma( fms( ux, cby, uy*cbx ), v3, uz );
    ux  = fma( fms( v1, cbz, v2*cby ), v4, ux );
    uy  = fma( fms( v2, cbx, v0*cbz ), v4, uy );
    uz  = fma( fms( v0, cby, v1*cbx ), v4, uz );
    ux += hax;
    uy += hay;
    uz += haz;

    store_8x4_tr( ux, uy, uz, q,
		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		  &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux );

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

    outbnd = ( v3 > one ) | ( v3 < neg_one ) |
             ( v4 > one ) | ( v4 < neg_one ) |
             ( v5 > one ) | ( v5 < neg_one );

    v3  = merge( outbnd, dx, v3 ); // Do not update outbnd particles
    v4  = merge( outbnd, dy, v4 );
    v5  = merge( outbnd, dz, v5 );
 
    store_8x4_tr( v3, v4, v5, ii,
		  &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		  &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx );

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q  = czero( outbnd, q*qsp );      // Do not accumulate outbnd particles
    dx = v0;                          // Streak midpoint (valid for inbnd only)
    dy = v1;
    dz = v2;

    v8 = q*ux*uy*uz*one_third;     // Charge conservation correction

    vp0 = ( float * ALIGNED(16) ) ( a0 + ii(0) ); // Accumulator pointers
    vp1 = ( float * ALIGNED(16) ) ( a0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( a0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( a0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( a0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( a0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( a0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( a0 + ii(7) );

    // In this macro, we need to set the v4-v7 vectors to zero so that when
    // the transpose is performed, the last 4 entries in each vector will be
    // zero because there are only 4 values per particle to load into the
    // accumulator arrays. There are other ways to accomplish the same result.
    // One would be to load the contents of vectors v0-v3 into 8 V4 arrays
    // and use the functionality of the V4 class. Not sure which is optimal.
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
    v0 += v8;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v8;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v8;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v8;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */      \
    v4  = 0.0;      /* Zero pad                             */	    \
    v5  = 0.0;      /* Zero pad                             */	    \
    v6  = 0.0;      /* Zero pad                             */	    \
    v7  = 0.0;      /* Zero pad                             */	    \
    transpose( v0, v1, v2, v3, v4, v5, v6, v7 );                    \
    increment_8x1( vp0 + offset, v0 );                              \
    increment_8x1( vp1 + offset, v1 );                              \
    increment_8x1( vp2 + offset, v2 );                              \
    increment_8x1( vp3 + offset, v3 );                              \
    increment_8x1( vp4 + offset, v4 );                              \
    increment_8x1( vp5 + offset, v5 );                              \
    increment_8x1( vp6 + offset, v6 );                              \
    increment_8x1( vp7 + offset, v7 );

    ACCUMULATE_J( x, y, z, 0 );
    ACCUMULATE_J( y, z, x, 4 );
    ACCUMULATE_J( z, x, y, 8 );

#   undef ACCUMULATE_J

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
    if ( outbnd(N) )                                /* Unlikely */      \
    {                                                                   \
      local_pm->dispx = ux(N);                                          \
      local_pm->dispy = uy(N);                                          \
      local_pm->dispz = uz(N);                                          \
      local_pm->i     = ( p - p0 ) + N;                                 \
      if ( move_p( p0, local_pm, a0, g, _qsp ) )    /* Unlikely */      \
      {                                                                 \
        if ( nm<max_nm )                                                \
        {                                                               \
	  v4::copy_4x1( &pm[nm++], local_pm ); 	                        \
	}                                                               \
        else                                        /* Unlikely */      \
	{                                                               \
	  itmp++;                                                       \
	}                                                               \
      }                                                                 \
    }

    MOVE_OUTBND(0);
    MOVE_OUTBND(1);
    MOVE_OUTBND(2);
    MOVE_OUTBND(3);
    MOVE_OUTBND(4);
    MOVE_OUTBND(5);
    MOVE_OUTBND(6);
    MOVE_OUTBND(7);

#   undef MOVE_OUTBND
  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
#endif // Method 1

#if 0
// Method 2.
void
advance_p_pipeline_v8( advance_p_pipeline_args_t * args,
		       int pipeline_rank,
		       int n_pipeline )
{
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
  float                * ALIGNED(16)  vp4;
  float                * ALIGNED(16)  vp5;
  float                * ALIGNED(16)  vp6;
  float                * ALIGNED(16)  vp7;

  const v8float qdt_2mc(args->qdt_2mc);
  const v8float cdt_dx(args->cdt_dx);
  const v8float cdt_dy(args->cdt_dy);
  const v8float cdt_dz(args->cdt_dz);
  const v8float qsp(args->qsp);
  const v8float one(1.);
  const v8float one_third(1./3.);
  const v8float two_fifteenths(2./15.);
  const v8float neg_one(-1.);

  const float _qsp = args->qsp;

  v8float dx, dy, dz, ux, uy, uz, q;
  v8float hax, hay, haz, cbx, cby, cbz;
  v8float v0, v1, v2, v3, v4, v5, v6, v7, v8, v9;
  v8int   ii, outbnd;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 3;

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

  // Process the particle quads for this pipeline.

  for( ; nq; nq--, p+=8 )
  {
    load_8x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		 &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
		 dx, dy, dz, ii );

    // Interpolate fields.
    vp0 = ( float * ALIGNED(16) ) ( f0 + ii(0) );
    vp1 = ( float * ALIGNED(16) ) ( f0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( f0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( f0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( f0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( f0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( f0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( f0 + ii(7) );

    load_8x4_tr( vp0, vp1, vp2, vp3,
		 vp4, vp5, vp6, vp7,
		 hax, v0, v1, v2 );
    hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );

    load_8x4_tr( vp0+4, vp1+4, vp2+4, vp3+4,
		 vp4+4, vp5+4, vp6+4, vp7+4,
		 hay, v3, v4, v5 );
    hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );

    load_8x4_tr( vp0+8, vp1+8, vp2+8, vp3+8,
		 vp4+8, vp5+8, vp6+8, vp7+8,
		 haz, v0, v1, v2 );
    haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );

    load_8x4_tr( vp0+12, vp1+12, vp2+12, vp3+12,
		 vp4+12, vp5+12, vp6+12, vp7+12,
		 cbx, v3, cby, v4 );
    cbx = fma( v3, dx, cbx );
    cby = fma( v4, dy, cby );

    load_8x2_tr( vp0+16, vp1+16, vp2+16, vp3+16,
		 vp4+16, vp5+16, vp6+16, vp7+16,
		 cbz, v5 );
    cbz = fma( v5, dz, cbz );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    load_8x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		 &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux,
		 ux, uy, uz, q );

    ux += hax;
    uy += hay;
    uz += haz;
    v0  = qdt_2mc*rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
    v1  = fma( cbx,cbx, fma( cby,cby, cbz*cbz ) );
    v2  = (v0*v0)*v1;
    v3  = v0*fma( fma( two_fifteenths, v2, one_third ), v2, one );
    v4  = v3*rcp(fma( v3*v3, v1, one ));
    v4 += v4;
    v0  = fma( fms( uy, cbz, uz*cby ), v3, ux );
    v1  = fma( fms( uz, cbx, ux*cbz ), v3, uy );
    v2  = fma( fms( ux, cby, uy*cbx ), v3, uz );
    ux  = fma( fms( v1, cbz, v2*cby ), v4, ux );
    uy  = fma( fms( v2, cbx, v0*cbz ), v4, uy );
    uz  = fma( fms( v0, cby, v1*cbx ), v4, uz );
    ux += hax;
    uy += hay;
    uz += haz;

    store_8x4_tr( ux, uy, uz, q,
		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		  &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux );

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

    outbnd = ( v3 > one ) | ( v3 < neg_one ) |
             ( v4 > one ) | ( v4 < neg_one ) |
             ( v5 > one ) | ( v5 < neg_one );

    v3  = merge( outbnd, dx, v3 ); // Do not update outbnd particles
    v4  = merge( outbnd, dy, v4 );
    v5  = merge( outbnd, dz, v5 );

    store_8x4_tr( v3, v4, v5, ii,
		  &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		  &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx );

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q  = czero( outbnd, q*qsp );      // Do not accumulate outbnd particles
    dx = v0;                          // Streak midpoint (valid for inbnd only)
    dy = v1;
    dz = v2;

    v9 = q*ux*uy*uz*one_third;     // Charge conservation correction

    vp0 = ( float * ALIGNED(16) ) ( a0 + ii(0) ); // Accumulator pointers
    vp1 = ( float * ALIGNED(16) ) ( a0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( a0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( a0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( a0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( a0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( a0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( a0 + ii(7) );

#   define ACCUMULATE_JX(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v1  = v8*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v8-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v8;       /* v1 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v2  = v0*v8;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v8;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v0 *= v8;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v8;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v9;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v9;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v9;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v9;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JY(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v5  = v8*d##Y;  /* v5 = q ux dy                         */      \
    v4  = v8-v5;    /* v4 = q ux (1-dy)                     */      \
    v5 += v8;       /* v5 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v6  = v4*v8;    /* v6 = q ux (1-dy)(1+dz)               */      \
    v7  = v5*v8;    /* v7 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v4 *= v8;       /* v4 = q ux (1-dy)(1-dz)               */      \
    v5 *= v8;       /* v5 = q ux (1+dy)(1-dz)               */      \
    v4 += v9;       /* v4 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v5 -= v9;       /* v5 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v6 -= v9;       /* v6 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v7 += v9;       /* v7 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JZ(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v1  = v8*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v8-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v8;       /* v1 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v2  = v0*v8;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v8;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v0 *= v8;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v8;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v9;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v9;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v9;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v9;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

    // Accumulate Jx for 8 particles into the v0-v3 vectors.
    ACCUMULATE_JX( x, y, z );

    // Accumulate Jy for 8 particles into the v4-v7 vectors.
    ACCUMULATE_JY( y, z, x );

    // Transpose the data in vectors v0-v7 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v0, v1, v2, v3, v4, v5, v6, v7 );

    // Add the contributions to Jx and Jy from 8 particles into the
    // accumulator arrays for Jx and Jy.
    increment_8x1( vp0, v0 );
    increment_8x1( vp1, v1 );
    increment_8x1( vp2, v2 );
    increment_8x1( vp3, v3 );
    increment_8x1( vp4, v4 );
    increment_8x1( vp5, v5 );
    increment_8x1( vp6, v6 );
    increment_8x1( vp7, v7 );

    // Accumulate Jz for 8 particles into the v0-v3 vectors.
    ACCUMULATE_JZ( z, x, y );

    // Zero the v4-v7 vectors prior to transposing the data.
    v4 = 0.0;
    v5 = 0.0;
    v6 = 0.0;
    v7 = 0.0;

    // Transpose the data in vectors v0-v7 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v0, v1, v2, v3, v4, v5, v6, v7 );

    // Add the contributions to Jz from 8 particles into the accumulator
    // arrays for Jz.
    increment_8x1( vp0 + 8, v0 );
    increment_8x1( vp1 + 8, v1 );
    increment_8x1( vp2 + 8, v2 );
    increment_8x1( vp3 + 8, v3 );
    increment_8x1( vp4 + 8, v4 );
    increment_8x1( vp5 + 8, v5 );
    increment_8x1( vp6 + 8, v6 );
    increment_8x1( vp7 + 8, v7 );

#   undef ACCUMULATE_JX
#   undef ACCUMULATE_JY
#   undef ACCUMULATE_JZ

    // Update position and accumulate outbnd

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

    MOVE_OUTBND(0);
    MOVE_OUTBND(1);
    MOVE_OUTBND(2);
    MOVE_OUTBND(3);
    MOVE_OUTBND(4);
    MOVE_OUTBND(5);
    MOVE_OUTBND(6);
    MOVE_OUTBND(7);

#   undef MOVE_OUTBND
  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
#endif // Method 2

#if 0
// Method 3.
void
advance_p_pipeline_v8( advance_p_pipeline_args_t * args,
		       int pipeline_rank,
		       int n_pipeline )
{
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
  float                * ALIGNED(16)  vp4;
  float                * ALIGNED(16)  vp5;
  float                * ALIGNED(16)  vp6;
  float                * ALIGNED(16)  vp7;

  const v8float qdt_2mc(args->qdt_2mc);
  const v8float cdt_dx(args->cdt_dx);
  const v8float cdt_dy(args->cdt_dy);
  const v8float cdt_dz(args->cdt_dz);
  const v8float qsp(args->qsp);
  const v8float one(1.);
  const v8float one_third(1./3.);
  const v8float two_fifteenths(2./15.);
  const v8float neg_one(-1.);

  const float _qsp = args->qsp;

  v8float dx, dy, dz, ux, uy, uz, q;
  v8float hax, hay, haz, cbx, cby, cbz;
  v8float v0, v1, v2, v3, v4, v5, v6, v7, v8, v9;
  v8int   ii, outbnd;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 3;

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

  // Process the particle quads for this pipeline.

  for( ; nq; nq--, p+=8 )
  {
    load_8x8_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		 &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
		 dx, dy, dz, ii, ux, uy, uz, q );

    // Interpolate fields.
    vp0 = ( float * ALIGNED(16) ) ( f0 + ii(0) );
    vp1 = ( float * ALIGNED(16) ) ( f0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( f0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( f0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( f0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( f0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( f0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( f0 + ii(7) );

    load_8x8_tr( vp0, vp1, vp2, vp3,
		 vp4, vp5, vp6, vp7,
		 hax, v0, v1, v2, hay, v3, v4, v5 );

    hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );
    hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );

    load_8x8_tr( vp0+8, vp1+8, vp2+8, vp3+8,
		 vp4+8, vp5+8, vp6+8, vp7+8,
		 haz, v0, v1, v2, cbx, v3, cby, v4 );

    haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );
    cbx = fma( v3, dx, cbx );
    cby = fma( v4, dy, cby );

    load_8x2_tr( vp0+16, vp1+16, vp2+16, vp3+16,
                 vp4+16, vp5+16, vp6+16, vp7+16,
                 cbz, v5 );

    // Use this until I can implement load_8x2_tr. This is a kludge
    // and I assume will perform slower than load_8x2_tr when
    // implemented.
    // load_8x8_tr( vp0+16, vp1+16, vp2+16, vp3+16,
    // 		 vp4+16, vp5+16, vp6+16, vp7+16,
    // 		 cbz, v5, v0, v1, v2, v3, v4, v6 );

    cbz = fma( v5, dz, cbz );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    ux += hax;
    uy += hay;
    uz += haz;
    v0  = qdt_2mc*rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
    v1  = fma( cbx,cbx, fma( cby,cby, cbz*cbz ) );
    v2  = (v0*v0)*v1;
    v3  = v0*fma( fma( two_fifteenths, v2, one_third ), v2, one );
    v4  = v3*rcp(fma( v3*v3, v1, one ));
    v4 += v4;
    v0  = fma( fms( uy, cbz, uz*cby ), v3, ux );
    v1  = fma( fms( uz, cbx, ux*cbz ), v3, uy );
    v2  = fma( fms( ux, cby, uy*cbx ), v3, uz );
    ux  = fma( fms( v1, cbz, v2*cby ), v4, ux );
    uy  = fma( fms( v2, cbx, v0*cbz ), v4, uy );
    uz  = fma( fms( v0, cby, v1*cbx ), v4, uz );
    ux += hax;
    uy += hay;
    uz += haz;

    // Store ux, uy, uz in v6, v7, v8 so store_8x8_tr can be used below.
    v6  = ux;
    v7  = uy;
    v8  = uz;

    // store_8x4_tr( ux, uy, uz, q,
    // 		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
    // 		  &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux );

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

    outbnd = ( v3 > one ) | ( v3 < neg_one ) |
             ( v4 > one ) | ( v4 < neg_one ) |
             ( v5 > one ) | ( v5 < neg_one );

    v3  = merge( outbnd, dx, v3 ); // Do not update outbnd particles
    v4  = merge( outbnd, dy, v4 );
    v5  = merge( outbnd, dz, v5 );

    store_8x8_tr( v3, v4, v5, ii, v6, v7, v8, q,
		  &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		  &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx );

    // store_8x4_tr( v3, v4, v5, ii,
    // 		  &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
    // 		  &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx );

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q  = czero( outbnd, q*qsp );   // Do not accumulate outbnd particles
    dx = v0;                       // Streak midpoint (valid for inbnd only)
    dy = v1;
    dz = v2;

    v9 = q*ux*uy*uz*one_third;     // Charge conservation correction

    vp0 = ( float * ALIGNED(16) ) ( a0 + ii(0) ); // Accumulator pointers
    vp1 = ( float * ALIGNED(16) ) ( a0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( a0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( a0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( a0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( a0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( a0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( a0 + ii(7) );

#   define ACCUMULATE_JX(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v1  = v8*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v8-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v8;       /* v1 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v2  = v0*v8;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v8;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v0 *= v8;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v8;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v9;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v9;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v9;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v9;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JY(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v5  = v8*d##Y;  /* v5 = q ux dy                         */      \
    v4  = v8-v5;    /* v4 = q ux (1-dy)                     */      \
    v5 += v8;       /* v5 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v6  = v4*v8;    /* v6 = q ux (1-dy)(1+dz)               */      \
    v7  = v5*v8;    /* v7 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v4 *= v8;       /* v4 = q ux (1-dy)(1-dz)               */      \
    v5 *= v8;       /* v5 = q ux (1+dy)(1-dz)               */      \
    v4 += v9;       /* v4 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v5 -= v9;       /* v5 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v6 -= v9;       /* v6 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v7 += v9;       /* v7 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JZ(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v1  = v8*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v8-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v8;       /* v1 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v2  = v0*v8;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v8;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v0 *= v8;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v8;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v9;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v9;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v9;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v9;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

    // Accumulate Jx for 8 particles into the v0-v3 vectors.
    ACCUMULATE_JX( x, y, z );

    // Accumulate Jy for 8 particles into the v4-v7 vectors.
    ACCUMULATE_JY( y, z, x );

    // Transpose the data in vectors v0-v7 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v0, v1, v2, v3, v4, v5, v6, v7 );

    // Add the contributions to Jx and Jy from 8 particles into the
    // accumulator arrays for Jx and Jy.
    increment_8x1( vp0, v0 );
    increment_8x1( vp1, v1 );
    increment_8x1( vp2, v2 );
    increment_8x1( vp3, v3 );
    increment_8x1( vp4, v4 );
    increment_8x1( vp5, v5 );
    increment_8x1( vp6, v6 );
    increment_8x1( vp7, v7 );

    // Accumulate Jz for 8 particles into the v0-v3 vectors.
    ACCUMULATE_JZ( z, x, y );

    // Zero the v4-v7 vectors prior to transposing the data.
    v4 = 0.0;
    v5 = 0.0;
    v6 = 0.0;
    v7 = 0.0;

    // Transpose the data in vectors v0-v7 so it can be added into the
    // accumulator arrays using vector operations.
    transpose( v0, v1, v2, v3, v4, v5, v6, v7 );

    // Add the contributions to Jz from 8 particles into the accumulator
    // arrays for Jz.
    increment_8x1( vp0 + 8, v0 );
    increment_8x1( vp1 + 8, v1 );
    increment_8x1( vp2 + 8, v2 );
    increment_8x1( vp3 + 8, v3 );
    increment_8x1( vp4 + 8, v4 );
    increment_8x1( vp5 + 8, v5 );
    increment_8x1( vp6 + 8, v6 );
    increment_8x1( vp7 + 8, v7 );

#   undef ACCUMULATE_JX
#   undef ACCUMULATE_JY
#   undef ACCUMULATE_JZ

    // Update position and accumulate outbnd

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

    MOVE_OUTBND(0);
    MOVE_OUTBND(1);
    MOVE_OUTBND(2);
    MOVE_OUTBND(3);
    MOVE_OUTBND(4);
    MOVE_OUTBND(5);
    MOVE_OUTBND(6);
    MOVE_OUTBND(7);

#   undef MOVE_OUTBND
  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
#endif // Method 3

#if 0
// Method 4.
void
advance_p_pipeline_v8( advance_p_pipeline_args_t * args,
		       int pipeline_rank,
		       int n_pipeline )
{
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
  float                * ALIGNED(16)  vp4;
  float                * ALIGNED(16)  vp5;
  float                * ALIGNED(16)  vp6;
  float                * ALIGNED(16)  vp7;

  const v8float qdt_2mc(args->qdt_2mc);
  const v8float cdt_dx(args->cdt_dx);
  const v8float cdt_dy(args->cdt_dy);
  const v8float cdt_dz(args->cdt_dz);
  const v8float qsp(args->qsp);
  const v8float one(1.);
  const v8float one_third(1./3.);
  const v8float two_fifteenths(2./15.);
  const v8float neg_one(-1.);

  const float _qsp = args->qsp;

  v8float dx, dy, dz, ux, uy, uz, q;
  v8float hax, hay, haz, cbx, cby, cbz;
  v8float v0, v1, v2, v3, v4, v5, v6, v7, v8, v9;
  v8int   ii, outbnd;

  __m256 t00, t01, t02, t03, t04, t05;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 3;

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

  // Process the particle quads for this pipeline.

  for( ; nq; nq--, p+=8 )
  {
    //--------------------------------------------------------------------------//
    load_8x8_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		 &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
		 dx, dy, dz, ii, ux, uy, uz, q );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // Interpolate fields.
    //--------------------------------------------------------------------------//

    vp0 = ( float * ALIGNED(16) ) ( f0 + ii(0) );
    vp1 = ( float * ALIGNED(16) ) ( f0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( f0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( f0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( f0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( f0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( f0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( f0 + ii(7) );

    //--------------------------------------------------------------------------//
    load_8x8_tr( vp0, vp1, vp2, vp3,
		 vp4, vp5, vp6, vp7,
		 hax, v0, v1, v2, hay, v3, v4, v5 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );
    //--------------------------------------------------------------------------//

    hax.v = _mm256_mul_ps( qdt_2mc.v,
			   _mm256_fmadd_ps( _mm256_fmadd_ps( v2.v,
							     dy.v,
							     v1.v ),
					    dz.v,
					    _mm256_fmadd_ps( v0.v,
							     dy.v,
							     hax.v ) ) );

    //--------------------------------------------------------------------------//
    // hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );
    //--------------------------------------------------------------------------//

    hay.v = _mm256_mul_ps( qdt_2mc.v,
			   _mm256_fmadd_ps( _mm256_fmadd_ps( v5.v,
							     dz.v,
							     v4.v ),
					    dx.v,
					    _mm256_fmadd_ps( v3.v,
							     dz.v,
							     hay.v ) ) );

    //--------------------------------------------------------------------------//
    load_8x8_tr( vp0+8, vp1+8, vp2+8, vp3+8,
		 vp4+8, vp5+8, vp6+8, vp7+8,
		 haz, v0, v1, v2, cbx, v3, cby, v4 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );
    //--------------------------------------------------------------------------//

    haz.v = _mm256_mul_ps( qdt_2mc.v,
			   _mm256_fmadd_ps( _mm256_fmadd_ps( v2.v,
							     dx.v,
							     v1.v ),
					    dy.v,
					    _mm256_fmadd_ps( v0.v,
							     dx.v,
							     haz.v ) ) );

    //--------------------------------------------------------------------------//
    // cbx = fma( v3, dx, cbx );
    //--------------------------------------------------------------------------//

    cbx.v = _mm256_fmadd_ps( v3.v, dx.v, cbx.v );

    //--------------------------------------------------------------------------//
    // cby = fma( v4, dy, cby );
    //--------------------------------------------------------------------------//

    cby.v = _mm256_fmadd_ps( v4.v, dy.v, cby.v );

    //--------------------------------------------------------------------------//
    // load_8x2_tr( vp0+16, vp1+16, vp2+16, vp3+16,
    //              vp4+16, vp5+16, vp6+16, vp7+16,
    //              cbz, v5 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // Use this until I can implement load_8x2_tr. This is a kludge
    // and I assume will perform slower than load_8x2_tr when
    // implemented.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    load_8x8_tr( vp0+16, vp1+16, vp2+16, vp3+16,
		 vp4+16, vp5+16, vp6+16, vp7+16,
		 cbz, v5, v0, v1, v2, v3, v4, v6 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // cbz = fma( v5, dz, cbz );
    //--------------------------------------------------------------------------//

    cbz.v = _mm256_fmadd_ps( v5.v, dz.v, cbz.v );

    //--------------------------------------------------------------------------//
    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // ux += hax;
    //--------------------------------------------------------------------------//

    ux.v = _mm256_add_ps( ux.v, hax.v );

    //--------------------------------------------------------------------------//
    // uy += hay;
    //--------------------------------------------------------------------------//

    uy.v = _mm256_add_ps( uy.v, hay.v );

    //--------------------------------------------------------------------------//
    // uz += haz;
    //--------------------------------------------------------------------------//

    uz.v = _mm256_add_ps( uz.v, haz.v );

    //--------------------------------------------------------------------------//
    // v0  = qdt_2mc*rsqrt( one + fma( ux, ux, fma( uy, uy, uz*uz ) ) );
    //--------------------------------------------------------------------------//

    t00 = _mm256_add_ps( one.v,
			 _mm256_fmadd_ps( ux.v,
					  ux.v,
					  _mm256_fmadd_ps( uy.v,
							   uy.v,
							   _mm256_mul_ps( uz.v,
									  uz.v )
							   )
					  )
			 );

    t01 = _mm256_rsqrt_ps( t00 );

    t02 = _mm256_fmadd_ps( _mm256_set1_ps( 0.5f ),
			   _mm256_fnmadd_ps( t00,
					     _mm256_mul_ps( t01,
							    _mm256_mul_ps( t01,
									   t01 ) ),
					     t01 ),
			   t01 );

    v0.v = _mm256_mul_ps( qdt_2mc.v, t02 );

    //--------------------------------------------------------------------------//
    // v1  = fma( cbx, cbx, fma( cby, cby, cbz*cbz ) );
    //--------------------------------------------------------------------------//

    v1.v = _mm256_fmadd_ps( cbx.v,
			    cbx.v,
			    _mm256_fmadd_ps( cby.v,
					     cby.v,
					     _mm256_mul_ps( cbz.v,
							    cbz.v ) ) );

    //--------------------------------------------------------------------------//
    // v2  = (v0*v0)*v1;
    //--------------------------------------------------------------------------//

    v2.v = _mm256_mul_ps( _mm256_mul_ps( v0.v, v0.v ), v1.v );

    //--------------------------------------------------------------------------//
    // v3  = v0*fma( fma( two_fifteenths, v2, one_third ), v2, one );
    //--------------------------------------------------------------------------//

    v3.v = _mm256_mul_ps( v0.v,
			  _mm256_fmadd_ps( _mm256_fmadd_ps( two_fifteenths.v,
							    v2.v,
							    one_third.v ),
					   v2.v,
					   one.v ) );

    //--------------------------------------------------------------------------//
    // v4  = v3*rcp( fma( v3*v3, v1, one ) );
    //--------------------------------------------------------------------------//

    t00 = _mm256_fmadd_ps( _mm256_mul_ps( v3.v, v3.v ), v1.v, one.v );

    t01 = _mm256_rcp_ps( t00 );

    v4.v = _mm256_mul_ps( v3.v,
			  _mm256_fnmadd_ps( t00,
					    _mm256_mul_ps( t01, t01 ),
					    _mm256_add_ps( t01, t01 ) ) );

    //--------------------------------------------------------------------------//
    // v4 += v4;
    //--------------------------------------------------------------------------//

    v4.v = _mm256_add_ps( v4.v, v4.v );

    //--------------------------------------------------------------------------//
    // v0  = fma( fms( uy, cbz, uz*cby ), v3, ux );
    //--------------------------------------------------------------------------//

    v0.v = _mm256_fmadd_ps( _mm256_fmsub_ps( uy.v,
					     cbz.v,
					     _mm256_mul_ps( uz.v, cby.v ) ),
			    v3.v,
			    ux.v );

    //--------------------------------------------------------------------------//
    // v1  = fma( fms( uz, cbx, ux*cbz ), v3, uy );
    //--------------------------------------------------------------------------//

    v1.v = _mm256_fmadd_ps( _mm256_fmsub_ps( uz.v,
					     cbx.v,
					     _mm256_mul_ps( ux.v, cbz.v ) ),
			    v3.v,
			    uy.v );

    //--------------------------------------------------------------------------//
    // v2  = fma( fms( ux, cby, uy*cbx ), v3, uz );
    //--------------------------------------------------------------------------//

    v2.v = _mm256_fmadd_ps( _mm256_fmsub_ps( ux.v,
					     cby.v,
					     _mm256_mul_ps( uy.v, cbx.v ) ),
			    v3.v,
			    uz.v );

    //--------------------------------------------------------------------------//
    // ux  = fma( fms( v1, cbz, v2*cby ), v4, ux );
    //--------------------------------------------------------------------------//

    ux.v = _mm256_fmadd_ps( _mm256_fmsub_ps( v1.v,
					     cbz.v,
					     _mm256_mul_ps( v2.v, cby.v ) ),
			    v4.v,
			    ux.v );

    //--------------------------------------------------------------------------//
    // uy  = fma( fms( v2, cbx, v0*cbz ), v4, uy );
    //--------------------------------------------------------------------------//

    uy.v = _mm256_fmadd_ps( _mm256_fmsub_ps( v2.v,
					     cbx.v,
					     _mm256_mul_ps( v0.v, cbz.v ) ),
			    v4.v,
			    uy.v );

    //--------------------------------------------------------------------------//
    // uz  = fma( fms( v0, cby, v1*cbx ), v4, uz );
    //--------------------------------------------------------------------------//

    uz.v = _mm256_fmadd_ps( _mm256_fmsub_ps( v0.v,
					     cby.v,
					     _mm256_mul_ps( v1.v, cbx.v ) ),
			    v4.v,
			    uz.v );

    //--------------------------------------------------------------------------//
    // ux += hax;
    //--------------------------------------------------------------------------//

    ux.v = _mm256_add_ps( ux.v, hax.v );

    //--------------------------------------------------------------------------//
    // uy += hay;
    //--------------------------------------------------------------------------//

    uy.v = _mm256_add_ps( uy.v, hay.v );

    //--------------------------------------------------------------------------//
    // uz += haz;
    //--------------------------------------------------------------------------//

    uz.v = _mm256_add_ps( uz.v, haz.v );

    //--------------------------------------------------------------------------//
    // Store ux, uy, uz in v6, v7, v8 so store_8x8_tr can be used below.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v6  = ux;
    //--------------------------------------------------------------------------//

    v6.v = ux.v;

    //--------------------------------------------------------------------------//
    // v7  = uy;
    //--------------------------------------------------------------------------//

    v7.v = uy.v;

    //--------------------------------------------------------------------------//
    // v8  = uz;
    //--------------------------------------------------------------------------//

    v8.v = uz.v;

    //--------------------------------------------------------------------------//
    // store_8x4_tr( ux, uy, uz, q,
    //  	     &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
    // 		     &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // Update the position of inbnd particles
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v0  = rsqrt( one + fma( ux, ux, fma( uy, uy, uz*uz ) ) );
    //--------------------------------------------------------------------------//

    t00 = _mm256_add_ps( one.v,
			 _mm256_fmadd_ps( ux.v,
					  ux.v,
					  _mm256_fmadd_ps( uy.v,
							   uy.v,
							   _mm256_mul_ps( uz.v,
									  uz.v )
							   )
					  )
			 );

    t01 = _mm256_rsqrt_ps( t00 );

    v0.v = _mm256_fmadd_ps( _mm256_set1_ps( 0.5f ),
			    _mm256_fnmadd_ps( t00,
					      _mm256_mul_ps( t01,
							     _mm256_mul_ps( t01,
									    t01 ) ),
					      t01 ),
			    t01 );

    //--------------------------------------------------------------------------//
    // ux *= cdt_dx;
    //--------------------------------------------------------------------------//

    ux.v = _mm256_mul_ps( ux.v, cdt_dx.v );

    //--------------------------------------------------------------------------//
    // uy *= cdt_dy;
    //--------------------------------------------------------------------------//

    uy.v = _mm256_mul_ps( uy.v, cdt_dy.v );

    //--------------------------------------------------------------------------//
    // uz *= cdt_dz;
    //--------------------------------------------------------------------------//

    uz.v = _mm256_mul_ps( uz.v, cdt_dz.v );

    //--------------------------------------------------------------------------//
    // ux *= v0;
    //--------------------------------------------------------------------------//

    ux.v = _mm256_mul_ps( ux.v, v0.v );

    //--------------------------------------------------------------------------//
    // uy *= v0;
    //--------------------------------------------------------------------------//

    uy.v = _mm256_mul_ps( uy.v, v0.v );

    //--------------------------------------------------------------------------//
    // uz *= v0;      // ux,uy,uz are normalized displ (relative to cell size)
    //--------------------------------------------------------------------------//

    uz.v = _mm256_mul_ps( uz.v, v0.v );

    //--------------------------------------------------------------------------//
    // v0  = dx + ux;
    //--------------------------------------------------------------------------//

    v0.v = _mm256_add_ps( dx.v, ux.v );

    //--------------------------------------------------------------------------//
    // v1  = dy + uy;
    //--------------------------------------------------------------------------//

    v1.v = _mm256_add_ps( dy.v, uy.v );

    //--------------------------------------------------------------------------//
    // v2  = dz + uz; // New particle midpoint
    //--------------------------------------------------------------------------//

    v2.v = _mm256_add_ps( dz.v, uz.v );

    //--------------------------------------------------------------------------//
    // v3  = v0 + ux;
    //--------------------------------------------------------------------------//

    v3.v = _mm256_add_ps( v0.v, ux.v );

    //--------------------------------------------------------------------------//
    // v4  = v1 + uy;
    //--------------------------------------------------------------------------//

    v4.v = _mm256_add_ps( v1.v, uy.v );

    //--------------------------------------------------------------------------//
    // v5  = v2 + uz; // New particle position
    //--------------------------------------------------------------------------//

    v5.v = _mm256_add_ps( v2.v, uz.v );

    //--------------------------------------------------------------------------//
    outbnd = ( v3 > one ) | ( v3 < neg_one ) |
             ( v4 > one ) | ( v4 < neg_one ) |
             ( v5 > one ) | ( v5 < neg_one );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v3  = merge( outbnd, dx, v3 ); // Do not update outbnd particles
    //--------------------------------------------------------------------------//

    v3.v = _mm256_or_ps( _mm256_andnot_ps( outbnd.v,
					   v3.v ),
			 _mm256_and_ps( outbnd.v,
					dx.v ) );

    //--------------------------------------------------------------------------//
    // v4  = merge( outbnd, dy, v4 );
    //--------------------------------------------------------------------------//

    v4.v = _mm256_or_ps( _mm256_andnot_ps( outbnd.v,
					   v4.v ),
			 _mm256_and_ps( outbnd.v,
					dy.v ) );

    //--------------------------------------------------------------------------//
    // v5  = merge( outbnd, dz, v5 );
    //--------------------------------------------------------------------------//

    v5.v = _mm256_or_ps( _mm256_andnot_ps( outbnd.v,
					   v5.v ),
			 _mm256_and_ps( outbnd.v,
					dz.v ) );

    //--------------------------------------------------------------------------//
    store_8x8_tr( v3, v4, v5, ii, v6, v7, v8, q,
		  &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		  &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx );
    //--------------------------------------------------------------------------//

    // store_8x4_tr( v3, v4, v5, ii,
    // 		  &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
    // 		  &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx );

    //--------------------------------------------------------------------------//
    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // q  = czero( outbnd, q*qsp );   // Do not accumulate outbnd particles
    //--------------------------------------------------------------------------//

    q.v = _mm256_andnot_ps( outbnd.v,
			    _mm256_mul_ps( q.v, qsp.v ) );

    //--------------------------------------------------------------------------//
    // dx = v0;                       // Streak midpoint (valid for inbnd only)
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
    // v9 = q*ux*uy*uz*one_third;     // Charge conservation correction
    //--------------------------------------------------------------------------//

    v9.v = _mm256_mul_ps( q.v,
			  _mm256_mul_ps( ux.v,
					 _mm256_mul_ps( uy.v,
							_mm256_mul_ps( uz.v,
								       one_third.v )
							)
					 )
			  );

    //--------------------------------------------------------------------------//
    // Accumulator pointers.
    //--------------------------------------------------------------------------//

    vp0 = ( float * ALIGNED(16) ) ( a0 + ii(0) ); // Accumulator pointers
    vp1 = ( float * ALIGNED(16) ) ( a0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( a0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( a0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( a0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( a0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( a0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( a0 + ii(7) );

#   define ACCUMULATE_JX(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v1  = v8*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v8-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v8;       /* v1 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v2  = v0*v8;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v8;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v0 *= v8;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v8;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v9;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v9;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v9;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v9;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JY(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v5  = v8*d##Y;  /* v5 = q ux dy                         */      \
    v4  = v8-v5;    /* v4 = q ux (1-dy)                     */      \
    v5 += v8;       /* v5 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v6  = v4*v8;    /* v6 = q ux (1-dy)(1+dz)               */      \
    v7  = v5*v8;    /* v7 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v4 *= v8;       /* v4 = q ux (1-dy)(1-dz)               */      \
    v5 *= v8;       /* v5 = q ux (1+dy)(1-dz)               */      \
    v4 += v9;       /* v4 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v5 -= v9;       /* v5 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v6 -= v9;       /* v6 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v7 += v9;       /* v7 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JZ(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v1  = v8*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v8-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v8;       /* v1 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v2  = v0*v8;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v8;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v0 *= v8;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v8;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v9;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v9;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v9;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v9;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

    //--------------------------------------------------------------------------//
    // Accumulate Jx for 8 particles into the v0-v3 vectors.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // ACCUMULATE_JX( x, y, z );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v8  = q*ux;
    //--------------------------------------------------------------------------//

    v8.v = _mm256_mul_ps( q.v, ux.v );

    //--------------------------------------------------------------------------//
    // v1  = v8*dy;
    //--------------------------------------------------------------------------//

    v1.v = _mm256_mul_ps( v8.v, dy.v );

    //--------------------------------------------------------------------------//
    // v0  = v8-v1;
    //--------------------------------------------------------------------------//

    v0.v = _mm256_sub_ps( v8.v, v1.v );

    //--------------------------------------------------------------------------//
    // v1 += v8;
    //--------------------------------------------------------------------------//

    v1.v = _mm256_add_ps( v1.v, v8.v );

    //--------------------------------------------------------------------------//
    // v8  = one+dz;
    //--------------------------------------------------------------------------//

    v8.v = _mm256_add_ps( one.v, dz.v );

    //--------------------------------------------------------------------------//
    // v2  = v0*v8;
    //--------------------------------------------------------------------------//

    v2.v = _mm256_mul_ps( v0.v, v8.v );

    //--------------------------------------------------------------------------//
    // v3  = v1*v8;
    //--------------------------------------------------------------------------//

    v3.v = _mm256_mul_ps( v1.v, v8.v );

    //--------------------------------------------------------------------------//
    // v8  = one-dz;
    //--------------------------------------------------------------------------//

    v8.v = _mm256_sub_ps( one.v, dz.v );

    //--------------------------------------------------------------------------//
    // v0 *= v8;
    //--------------------------------------------------------------------------//

    v0.v = _mm256_mul_ps( v0.v, v8.v );

    //--------------------------------------------------------------------------//
    // v1 *= v8;
    //--------------------------------------------------------------------------//

    v1.v = _mm256_mul_ps( v1.v, v8.v );

    //--------------------------------------------------------------------------//
    // v0 += v9;
    //--------------------------------------------------------------------------//

    v0.v = _mm256_add_ps( v0.v, v9.v );

    //--------------------------------------------------------------------------//
    // v1 -= v9;
    //--------------------------------------------------------------------------//

    v1.v = _mm256_sub_ps( v1.v, v9.v );

    //--------------------------------------------------------------------------//
    // v2 -= v9;
    //--------------------------------------------------------------------------//

    v2.v = _mm256_sub_ps( v2.v, v9.v );

    //--------------------------------------------------------------------------//
    // v3 += v9;
    //--------------------------------------------------------------------------//

    v3.v = _mm256_add_ps( v3.v, v9.v );

    //--------------------------------------------------------------------------//
    // Accumulate Jy for 8 particles into the v4-v7 vectors.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // ACCUMULATE_JY( y, z, x );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v8  = q*uy;
    //--------------------------------------------------------------------------//

    v8.v = _mm256_mul_ps( q.v, uy.v );

    //--------------------------------------------------------------------------//
    // v5  = v8*dz;
    //--------------------------------------------------------------------------//

    v5.v = _mm256_mul_ps( v8.v, dz.v );

    //--------------------------------------------------------------------------//
    // v4  = v8-v5;
    //--------------------------------------------------------------------------//

    v4.v = _mm256_sub_ps( v8.v, v5.v );

    //--------------------------------------------------------------------------//
    // v5 += v8;
    //--------------------------------------------------------------------------//

    v5.v = _mm256_add_ps( v5.v, v8.v );

    //--------------------------------------------------------------------------//
    // v8  = one+dx;
    //--------------------------------------------------------------------------//

    v8.v = _mm256_add_ps( one.v, dx.v );

    //--------------------------------------------------------------------------//
    // v6  = v4*v8;
    //--------------------------------------------------------------------------//

    v6.v = _mm256_mul_ps( v04v, v8.v );

    //--------------------------------------------------------------------------//
    // v7  = v5*v8;
    //--------------------------------------------------------------------------//

    v7.v = _mm256_mul_ps( v5.v, v8.v );

    //--------------------------------------------------------------------------//
    // v8  = one-dx;
    //--------------------------------------------------------------------------//

    v8.v = _mm256_sub_ps( one.v, dx.v );

    //--------------------------------------------------------------------------//
    // v4 *= v8;
    //--------------------------------------------------------------------------//

    v4.v = _mm256_mul_ps( v4.v, v8.v );

    //--------------------------------------------------------------------------//
    // v5 *= v8;
    //--------------------------------------------------------------------------//

    v5.v = _mm256_mul_ps( v5.v, v8.v );

    //--------------------------------------------------------------------------//
    // v4 += v9;
    //--------------------------------------------------------------------------//

    v4.v = _mm256_add_ps( v4.v, v9.v );

    //--------------------------------------------------------------------------//
    // v5 -= v9;
    //--------------------------------------------------------------------------//

    v5.v = _mm256_sub_ps( v5.v, v9.v );

    //--------------------------------------------------------------------------//
    // v6 -= v9;
    //--------------------------------------------------------------------------//

    v6.v = _mm256_sub_ps( v6.v, v9.v );

    //--------------------------------------------------------------------------//
    // v7 += v9;
    //--------------------------------------------------------------------------//

    v7.v = _mm256_add_ps( v7.v, v9.v );

    //--------------------------------------------------------------------------//
    // Transpose the data in vectors v0-v7 so it can be added into the
    // accumulator arrays using vector operations.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    transpose( v0, v1, v2, v3, v4, v5, v6, v7 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // Add the contributions to Jx and Jy from 8 particles into the
    // accumulator arrays for Jx and Jy.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp0, v0 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp1, v1 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp2, v2 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp3, v3 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp4, v4 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp5, v5 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp6, v6 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp7, v7 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // Accumulate Jz for 8 particles into the v0-v3 vectors.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // ACCUMULATE_JZ( z, x, y );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // v8  = q*uz;
    //--------------------------------------------------------------------------//

    v8.v = _mm256_mul_ps( q.v, uz.v );

    //--------------------------------------------------------------------------//
    // v1  = v8*dx;
    //--------------------------------------------------------------------------//

    v1.v = _mm256_mul_ps( v8.v, dx.v );

    //--------------------------------------------------------------------------//
    // v0  = v8-v1;
    //--------------------------------------------------------------------------//

    v0.v = _mm256_sub_ps( v8.v, v1.v );

    //--------------------------------------------------------------------------//
    // v1 += v8;
    //--------------------------------------------------------------------------//

    v1.v = _mm256_add_ps( v1.v, v8.v );

    //--------------------------------------------------------------------------//
    // v8  = one+dy;
    //--------------------------------------------------------------------------//

    v8.v = _mm256_add_ps( one.v, dy.v );

    //--------------------------------------------------------------------------//
    // v2  = v0*v8;
    //--------------------------------------------------------------------------//

    v2.v = _mm256_mul_ps( v0.v, v8.v );

    //--------------------------------------------------------------------------//
    // v3  = v1*v8;
    //--------------------------------------------------------------------------//

    v3.v = _mm256_mul_ps( v1.v, v8.v );

    //--------------------------------------------------------------------------//
    // v8  = one-dy;
    //--------------------------------------------------------------------------//

    v8.v = _mm256_sub_ps( one.v, dy.v );

    //--------------------------------------------------------------------------//
    // v0 *= v8;
    //--------------------------------------------------------------------------//

    v0.v = _mm256_mul_ps( v0.v, v8.v );

    //--------------------------------------------------------------------------//
    // v1 *= v8;
    //--------------------------------------------------------------------------//

    v1.v = _mm256_mul_ps( v1.v, v8.v );

    //--------------------------------------------------------------------------//
    // v0 += v9;
    //--------------------------------------------------------------------------//

    v0.v = _mm256_add_ps( v0.v, v9.v );

    //--------------------------------------------------------------------------//
    // v1 -= v9;
    //--------------------------------------------------------------------------//

    v1.v = _mm256_sub_ps( v1.v, v9.v );

    //--------------------------------------------------------------------------//
    // v2 -= v9;
    //--------------------------------------------------------------------------//

    v2.v = _mm256_sub_ps( v2.v, v9.v );

    //--------------------------------------------------------------------------//
    // v3 += v9;
    //--------------------------------------------------------------------------//

    v3.v = _mm256_add_ps( v3.v, v9.v );

    //--------------------------------------------------------------------------//
    // Zero the v4-v7 vectors prior to transposing the data.
    //--------------------------------------------------------------------------//

    v4 = 0.0;
    v5 = 0.0;
    v6 = 0.0;
    v7 = 0.0;

    //--------------------------------------------------------------------------//
    // Transpose the data in vectors v0-v7 so it can be added into the
    // accumulator arrays using vector operations.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    transpose( v0, v1, v2, v3, v4, v5, v6, v7 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    // Add the contributions to Jz from 8 particles into the accumulator
    // arrays for Jz.
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp0 + 8, v0 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp1 + 8, v1 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp2 + 8, v2 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp3 + 8, v3 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp4 + 8, v4 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp5 + 8, v5 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp6 + 8, v6 );
    //--------------------------------------------------------------------------//

    //--------------------------------------------------------------------------//
    increment_8x1( vp7 + 8, v7 );
    //--------------------------------------------------------------------------//

#   undef ACCUMULATE_JX
#   undef ACCUMULATE_JY
#   undef ACCUMULATE_JZ

    // Update position and accumulate outbnd

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

    MOVE_OUTBND(0);
    MOVE_OUTBND(1);
    MOVE_OUTBND(2);
    MOVE_OUTBND(3);
    MOVE_OUTBND(4);
    MOVE_OUTBND(5);
    MOVE_OUTBND(6);
    MOVE_OUTBND(7);

#   undef MOVE_OUTBND
  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
#endif // Method 4

#if 0
// Method 5.
void
advance_p_pipeline_v8( advance_p_pipeline_args_t * args,
		       int pipeline_rank,
		       int n_pipeline )
{
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
  float                * ALIGNED(16)  vp4;
  float                * ALIGNED(16)  vp5;
  float                * ALIGNED(16)  vp6;
  float                * ALIGNED(16)  vp7;

  const v8float qdt_2mc(args->qdt_2mc);
  const v8float cdt_dx(args->cdt_dx);
  const v8float cdt_dy(args->cdt_dy);
  const v8float cdt_dz(args->cdt_dz);
  const v8float qsp(args->qsp);
  const v8float one(1.);
  const v8float one_third(1./3.);
  const v8float two_fifteenths(2./15.);
  const v8float neg_one(-1.);

  const float _qsp = args->qsp;

  v8float dx, dy, dz, ux, uy, uz, q;
  v8float hax, hay, haz, cbx, cby, cbz;
  v8float v0, v1, v2, v3, v4, v5, v6, v7, v8, v9;
  v8int   ii, outbnd;

  int itmp, nq, nm, max_nm;

  DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

  // Determine which quads of particle quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 3;

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

  // Process the particle quads for this pipeline.

  for( ; nq; nq--, p+=8 )
  {
    load_8x8_tr_v0( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		    &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
		    dx, dy, dz, ii, ux, uy, uz, q );
    // transpose_v0( dx, dy, dz, ii, ux, uy, uz, q );

    // Interpolate fields.
    vp0 = ( float * ALIGNED(16) ) ( f0 + ii(0) );
    vp1 = ( float * ALIGNED(16) ) ( f0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( f0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( f0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( f0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( f0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( f0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( f0 + ii(7) );

    load_8x8_tr_v1( vp0, vp1, vp2, vp3,
		    vp4, vp5, vp6, vp7,
		    hax, v0, v1, v2, hay, v3, v4, v5 );
    transpose_v1( hax, v0, v1, v2, hay, v3, v4, v5 );

    hax = qdt_2mc*fma_v00( fma_v01( v2, dy, v1 ), dz, fma_v02( v0, dy, hax ) );
    hay = qdt_2mc*fma_v03( fma_v04( v5, dz, v4 ), dx, fma_v05( v3, dz, hay ) );

    load_8x8_tr_v2( vp0+8, vp1+8, vp2+8, vp3+8,
		    vp4+8, vp5+8, vp6+8, vp7+8,
		    haz, v0, v1, v2, cbx, v3, cby, v4 );
    transpose_v2( haz, v0, v1, v2, cbx, v3, cby, v4 );

    haz = qdt_2mc*fma_v06( fma_v07( v2, dx, v1 ), dy, fma_v08( v0, dx, haz ) );
    cbx = fma_v09( v3, dx, cbx );
    cby = fma_v10( v4, dy, cby );

    load_8x2_tr_v0( vp0+16, vp1+16, vp2+16, vp3+16,
		    vp4+16, vp5+16, vp6+16, vp7+16,
		    cbz, v5 );

    // Use this until I can implement load_8x2_tr. This is a kludge
    // and I assume will perform slower than load_8x2_tr when
    // implemented.
    // load_8x8_tr( vp0+16, vp1+16, vp2+16, vp3+16,
    // 		 vp4+16, vp5+16, vp6+16, vp7+16,
    // 		 cbz, v5, v0, v1, v2, v3, v4, v6 );

    cbz = fma_v11( v5, dz, cbz );

    // Update momentum
    // If you are willing to eat a 5-10% performance hit,
    // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
    // quite in the noise numerically) for cyclotron frequencies
    // approaching the nyquist frequency.

    ux += hax;
    uy += hay;
    uz += haz;
    v0  = qdt_2mc*rsqrt_v0( one + fma_v12( ux,ux, fma_v13( uy,uy, uz*uz ) ) );
    v1  = fma_v14( cbx,cbx, fma_v15( cby,cby, cbz*cbz ) );
    v2  = (v0*v0)*v1;
    v3  = v0*fma_v16( fma_v17( two_fifteenths, v2, one_third ), v2, one );
    v4  = v3*rcp(fma_v18( v3*v3, v1, one ));
    v4 += v4;
    v0  = fma_v19( fms( uy, cbz, uz*cby ), v3, ux );
    v1  = fma_v20( fms( uz, cbx, ux*cbz ), v3, uy );
    v2  = fma_v21( fms( ux, cby, uy*cbx ), v3, uz );
    ux  = fma_v22( fms( v1, cbz, v2*cby ), v4, ux );
    uy  = fma_v23( fms( v2, cbx, v0*cbz ), v4, uy );
    uz  = fma_v24( fms( v0, cby, v1*cbx ), v4, uz );
    ux += hax;
    uy += hay;
    uz += haz;

    // Store ux, uy, uz in v6, v7, v8 so store_8x8_tr can be used below.
    v6  = ux;
    v7  = uy;
    v8  = uz;

    // store_8x4_tr( ux, uy, uz, q,
    // 		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
    // 		  &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux );

    // Update the position of inbnd particles
    v0  = rsqrt_v1( one + fma_v25( ux,ux, fma_v26( uy,uy, uz*uz ) ) );
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

    outbnd = ( v3 > one ) | ( v3 < neg_one ) |
             ( v4 > one ) | ( v4 < neg_one ) |
             ( v5 > one ) | ( v5 < neg_one );

    v3  = merge( outbnd, dx, v3 ); // Do not update outbnd particles
    v4  = merge( outbnd, dy, v4 );
    v5  = merge( outbnd, dz, v5 );

    // transpose_v3( v3, v4, v5, ii, v6, v7, v8, q );
    store_8x8_tr_v0( v3, v4, v5, ii, v6, v7, v8, q,
		     &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		     &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx );

    // store_8x4_tr( v3, v4, v5, ii,
    // 		  &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
    // 		  &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx );

    // Accumulate current of inbnd particles.
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step.
    q  = czero( outbnd, q*qsp );   // Do not accumulate outbnd particles
    dx = v0;                       // Streak midpoint (valid for inbnd only)
    dy = v1;
    dz = v2;

    v9 = q*ux*uy*uz*one_third;     // Charge conservation correction

    vp0 = ( float * ALIGNED(16) ) ( a0 + ii(0) ); // Accumulator pointers
    vp1 = ( float * ALIGNED(16) ) ( a0 + ii(1) );
    vp2 = ( float * ALIGNED(16) ) ( a0 + ii(2) );
    vp3 = ( float * ALIGNED(16) ) ( a0 + ii(3) );
    vp4 = ( float * ALIGNED(16) ) ( a0 + ii(4) );
    vp5 = ( float * ALIGNED(16) ) ( a0 + ii(5) );
    vp6 = ( float * ALIGNED(16) ) ( a0 + ii(6) );
    vp7 = ( float * ALIGNED(16) ) ( a0 + ii(7) );

#   define ACCUMULATE_JX(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v1  = v8*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v8-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v8;       /* v1 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v2  = v0*v8;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v8;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v0 *= v8;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v8;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v9;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v9;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v9;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v9;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JY(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v5  = v8*d##Y;  /* v5 = q ux dy                         */      \
    v4  = v8-v5;    /* v4 = q ux (1-dy)                     */      \
    v5 += v8;       /* v5 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v6  = v4*v8;    /* v6 = q ux (1-dy)(1+dz)               */      \
    v7  = v5*v8;    /* v7 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v4 *= v8;       /* v4 = q ux (1-dy)(1-dz)               */      \
    v5 *= v8;       /* v5 = q ux (1+dy)(1-dz)               */      \
    v4 += v9;       /* v4 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v5 -= v9;       /* v5 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v6 -= v9;       /* v6 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v7 += v9;       /* v7 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

#   define ACCUMULATE_JZ(X,Y,Z)                                     \
    v8  = q*u##X;   /* v8 = q ux                            */      \
    v1  = v8*d##Y;  /* v1 = q ux dy                         */      \
    v0  = v8-v1;    /* v0 = q ux (1-dy)                     */      \
    v1 += v8;       /* v1 = q ux (1+dy)                     */      \
    v8  = one+d##Z; /* v8 = 1+dz                            */      \
    v2  = v0*v8;    /* v2 = q ux (1-dy)(1+dz)               */      \
    v3  = v1*v8;    /* v3 = q ux (1+dy)(1+dz)               */      \
    v8  = one-d##Z; /* v8 = 1-dz                            */      \
    v0 *= v8;       /* v0 = q ux (1-dy)(1-dz)               */      \
    v1 *= v8;       /* v1 = q ux (1+dy)(1-dz)               */      \
    v0 += v9;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */      \
    v1 -= v9;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */      \
    v2 -= v9;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */      \
    v3 += v9;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

    // Accumulate Jx for 8 particles into the v0-v3 vectors.
    ACCUMULATE_JX( x, y, z );

    // Accumulate Jy for 8 particles into the v4-v7 vectors.
    ACCUMULATE_JY( y, z, x );

    // Transpose the data in vectors v0-v7 so it can be added into the
    // accumulator arrays using vector operations.
    transpose_v4( v0, v1, v2, v3, v4, v5, v6, v7 );

    // Add the contributions to Jx and Jy from 8 particles into the
    // accumulator arrays for Jx and Jy.
    increment_8x1_v00( vp0, v0 );
    increment_8x1_v01( vp1, v1 );
    increment_8x1_v02( vp2, v2 );
    increment_8x1_v03( vp3, v3 );
    increment_8x1_v04( vp4, v4 );
    increment_8x1_v05( vp5, v5 );
    increment_8x1_v06( vp6, v6 );
    increment_8x1_v07( vp7, v7 );

    // Accumulate Jz for 8 particles into the v0-v3 vectors.
    ACCUMULATE_JZ( z, x, y );

    // Zero the v4-v7 vectors prior to transposing the data.
    v4 = 0.0;
    v5 = 0.0;
    v6 = 0.0;
    v7 = 0.0;

    // Transpose the data in vectors v0-v7 so it can be added into the
    // accumulator arrays using vector operations.
    transpose_v5( v0, v1, v2, v3, v4, v5, v6, v7 );

    // Add the contributions to Jz from 8 particles into the accumulator
    // arrays for Jz.
    increment_8x1_v08( vp0 + 8, v0 );
    increment_8x1_v09( vp1 + 8, v1 );
    increment_8x1_v10( vp2 + 8, v2 );
    increment_8x1_v11( vp3 + 8, v3 );
    increment_8x1_v12( vp4 + 8, v4 );
    increment_8x1_v13( vp5 + 8, v5 );
    increment_8x1_v14( vp6 + 8, v6 );
    increment_8x1_v15( vp7 + 8, v7 );

#   undef ACCUMULATE_JX
#   undef ACCUMULATE_JY
#   undef ACCUMULATE_JZ

    // Update position and accumulate outbnd

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

    MOVE_OUTBND(0);
    MOVE_OUTBND(1);
    MOVE_OUTBND(2);
    MOVE_OUTBND(3);
    MOVE_OUTBND(4);
    MOVE_OUTBND(5);
    MOVE_OUTBND(6);
    MOVE_OUTBND(7);

#   undef MOVE_OUTBND
  }

  args->seg[pipeline_rank].pm        = pm;
  args->seg[pipeline_rank].max_nm    = max_nm;
  args->seg[pipeline_rank].nm        = nm;
  args->seg[pipeline_rank].n_ignored = itmp;
}
#endif // Method 5
