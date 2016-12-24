using namespace v16;

void
center_p_pipeline_v16( center_p_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline )
{
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  particle_t           * ALIGNED(128) p;

  const float          * ALIGNED(16)  vp00;
  const float          * ALIGNED(16)  vp01;
  const float          * ALIGNED(16)  vp02;
  const float          * ALIGNED(16)  vp03;
  const float          * ALIGNED(16)  vp04;
  const float          * ALIGNED(16)  vp05;
  const float          * ALIGNED(16)  vp06;
  const float          * ALIGNED(16)  vp07;
  const float          * ALIGNED(16)  vp08;
  const float          * ALIGNED(16)  vp09;
  const float          * ALIGNED(16)  vp10;
  const float          * ALIGNED(16)  vp11;
  const float          * ALIGNED(16)  vp12;
  const float          * ALIGNED(16)  vp13;
  const float          * ALIGNED(16)  vp14;
  const float          * ALIGNED(16)  vp15;

  const v16float qdt_2mc(    args->qdt_2mc);
  const v16float qdt_4mc(0.5*args->qdt_2mc); // For half Boris rotate
  const v16float one(1.0);
  const v16float one_third(1.0/3.0);
  const v16float two_fifteenths(2.0/15.0);

  v16float dx, dy, dz, ux, uy, uz, q;
  v16float hax, hay, haz, cbx, cby, cbz;
  v16float v00, v01, v02, v03, v04, v05, v06, v07, v08, v09, v10;
  v16int   ii;

  int itmp, nq;

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 4;

  // Process the particle quads for this pipeline.

  for( ; nq; nq--, p+=16 )
  {
    load_16x8_tr_p( &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
                    &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		    dx, dy, dz, ii, ux, uy, uz, q );

    // load_8x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
    // 		 &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
    // 		 dx, dy, dz, ii );

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

    // // Interpolate fields
    // vp0 = (const float * ALIGNED(16))(f0 + ii(0));
    // vp1 = (const float * ALIGNED(16))(f0 + ii(1));
    // vp2 = (const float * ALIGNED(16))(f0 + ii(2));
    // vp3 = (const float * ALIGNED(16))(f0 + ii(3));
    // vp4 = (const float * ALIGNED(16))(f0 + ii(4));
    // vp5 = (const float * ALIGNED(16))(f0 + ii(5));
    // vp6 = (const float * ALIGNED(16))(f0 + ii(6));
    // vp7 = (const float * ALIGNED(16))(f0 + ii(7));

    //--------------------------------------------------------------------------
    // Load data for 16 particles for block a.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00, vp01, vp02, vp03,
		   vp04, vp05, vp06, vp07,
		   vp08, vp09, vp10, vp11,
		   vp12, vp13, vp14, vp15,
		   hax, v00, v01, v02, hay, v03, v04, v05,
		   haz, v06, v07, v08, cbx, v09, cby, v10 );

    // load_8x4_tr( vp0, vp1, vp2, vp3,
    // 		 vp4, vp5, vp6, vp7,
    // 		 hax, v0, v1, v2 );

    //--------------------------------------------------------------------------
    // Process particles for block a.
    //--------------------------------------------------------------------------
    hax = qdt_2mc*fma( fma( dy, v02, v01 ), dz, fma( dy, v00, hax ) );

    // load_8x4_tr( vp0+4, vp1+4, vp2+4, vp3+4,
    // 		 vp4+4, vp5+4, vp6+4, vp7+4,
    // 		 hay, v3, v4, v5 );

    hay = qdt_2mc*fma( fma( dz, v05, v04 ), dx, fma( dz, v03, hay ) );

    // load_8x4_tr( vp0+8, vp1+8, vp2+8, vp3+8,
    // 		 vp4+8, vp5+8, vp6+8, vp7+8,
    // 		 haz, v0, v1, v2 );

    haz = qdt_2mc*fma( fma( dx, v08, v07 ), dy, fma( dx, v06, haz ) );

    // load_8x4_tr( vp0+12, vp1+12, vp2+12, vp3+12,
    // 		 vp4+12, vp5+12, vp6+12, vp7+12,
    // 		 cbx, v3, cby, v4 );

    cbx = fma( v09, dx, cbx );
    cby = fma( v10, dy, cby );

    load_16x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
		  vp04+16, vp05+16, vp06+16, vp07+16,
		  vp08+16, vp09+16, vp10+16, vp11+16,
		  vp12+16, vp13+16, vp14+16, vp15+16,
		  cbz, v05 );

    // load_8x2_tr( vp0+16, vp1+16, vp2+16, vp3+16,
    // 		 vp4+16, vp5+16, vp6+16, vp7+16,
    // 		 cbz, v5 );

    cbz = fma( v05, dz, cbz );

    // Update momentum.  Could use load_8x3_tr.
    // load_8x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
    // 		 &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux,
    // 		 ux, uy, uz, q );

    ux  += hax;
    uy  += hay;
    uz  += haz;
    v00  = qdt_4mc*rsqrt( one + fma( ux, ux, fma( uy, uy, uz*uz ) ) );
    v01  = fma( cbx, cbx, fma( cby, cby, cbz*cbz ) );
    v02  = (v00*v00)*v01;
    v03  = v00*fma( v02, fma( v02, two_fifteenths, one_third ), one );
    v04  = v03*rcp( fma( v03*v03, v01, one ) );
    v04 += v04;
    v00  = fma( fms(  uy, cbz,  uz*cby ), v03, ux );
    v01  = fma( fms(  uz, cbx,  ux*cbz ), v03, uy );
    v02  = fma( fms(  ux, cby,  uy*cbx ), v03, uz );
    ux   = fma( fms( v01, cbz, v02*cby ), v04, ux );
    uy   = fma( fms( v02, cbx, v00*cbz ), v04, uy );
    uz   = fma( fms( v00, cby, v01*cbx ), v04, uz );

    // // Could use store_8x3_tr.
    // store_8x4_tr( ux, uy, uz, q,
    // 		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
    // 		  &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux );

    //--------------------------------------------------------------------------
    // Store the particle momenta.
    //--------------------------------------------------------------------------
    store_16x4_tr( ux, uy, uz, q,
                   &p[ 0].ux, &p[ 1].ux, &p[ 2].ux, &p[ 3].ux,
                   &p[ 4].ux, &p[ 5].ux, &p[ 6].ux, &p[ 7].ux,
                   &p[ 8].ux, &p[ 9].ux, &p[10].ux, &p[11].ux,
                   &p[12].ux, &p[13].ux, &p[14].ux, &p[15].ux );

    // store_16x4_tr_p( ux, uy, uz, q,
    //                  &p[ 0].ux, &p[ 2].ux, &p[ 4].ux, &p[ 6].ux,
    //                  &p[ 8].ux, &p[10].ux, &p[12].ux, &p[14].ux );
  }
}
