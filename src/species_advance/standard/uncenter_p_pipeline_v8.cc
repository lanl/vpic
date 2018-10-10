using namespace v8;

void
uncenter_p_pipeline_v8( center_p_pipeline_args_t * args,
                        int pipeline_rank,
                        int n_pipeline )
{
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  particle_t           * ALIGNED(128) p;

  const float          * ALIGNED(32)  vp00;
  const float          * ALIGNED(32)  vp01;
  const float          * ALIGNED(32)  vp02;
  const float          * ALIGNED(32)  vp03;
  const float          * ALIGNED(32)  vp04;
  const float          * ALIGNED(32)  vp05;
  const float          * ALIGNED(32)  vp06;
  const float          * ALIGNED(32)  vp07;

  const v8float qdt_2mc(    -args->qdt_2mc); // For backward half advance.
  const v8float qdt_4mc(-0.5*args->qdt_2mc); // For backward half Boris rotate.
  const v8float one(1.0);
  const v8float one_third(1.0/3.0);
  const v8float two_fifteenths(2.0/15.0);

  v8float dx, dy, dz, ux, uy, uz, q;
  v8float hax, hay, haz, cbx, cby, cbz;
  v8float v00, v01, v02, v03, v04, v05;
  v8int   ii;

  int first, nq;

  // Determine which particle quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, first, nq );

  p = args->p0 + first;

  nq >>= 3;

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=8 )
  {
    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------
    load_8x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		 &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
		 dx, dy, dz, ii );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( const float * ALIGNED(32) ) ( f0 + ii(0) );
    vp01 = ( const float * ALIGNED(32) ) ( f0 + ii(1) );
    vp02 = ( const float * ALIGNED(32) ) ( f0 + ii(2) );
    vp03 = ( const float * ALIGNED(32) ) ( f0 + ii(3) );
    vp04 = ( const float * ALIGNED(32) ) ( f0 + ii(4) );
    vp05 = ( const float * ALIGNED(32) ) ( f0 + ii(5) );
    vp06 = ( const float * ALIGNED(32) ) ( f0 + ii(6) );
    vp07 = ( const float * ALIGNED(32) ) ( f0 + ii(7) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_8x4_tr( vp00, vp01, vp02, vp03,
		 vp04, vp05, vp06, vp07,
		 hax, v00, v01, v02 );

    hax = qdt_2mc*fma( fma( dy, v02, v01 ), dz, fma( dy, v00, hax ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_8x4_tr( vp00+4, vp01+4, vp02+4, vp03+4,
		 vp04+4, vp05+4, vp06+4, vp07+4,
		 hay, v03, v04, v05 );

    hay = qdt_2mc*fma( fma( dz, v05, v04 ), dx, fma( dz, v03, hay ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_8x4_tr( vp00+8, vp01+8, vp02+8, vp03+8,
		 vp04+8, vp05+8, vp06+8, vp07+8,
		 haz, v00, v01, v02 );

    haz = qdt_2mc*fma( fma( dx, v02, v01 ), dy, fma( dx, v00, haz ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_8x4_tr( vp00+12, vp01+12, vp02+12, vp03+12,
		 vp04+12, vp05+12, vp06+12, vp07+12,
		 cbx, v03, cby, v04 );

    cbx = fma( v03, dx, cbx );
    cby = fma( v04, dy, cby );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles, final.
    //--------------------------------------------------------------------------
    load_8x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
		 vp04+16, vp05+16, vp06+16, vp07+16,
		 cbz, v05 );

    cbz = fma( v05, dz, cbz );

    //--------------------------------------------------------------------------
    // Load particle data.  Could use load_8x3_tr.
    //--------------------------------------------------------------------------
    load_8x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		 &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux,
		 ux, uy, uz, q );

    //--------------------------------------------------------------------------
    // Update momentum.
    //--------------------------------------------------------------------------
    v00  = qdt_4mc * rsqrt( one + fma( ux, ux, fma( uy, uy, uz * uz ) ) );
    v01  = fma( cbx, cbx, fma( cby, cby, cbz * cbz ) );
    v02  = ( v00 * v00 ) * v01;
    v03  = v00 * fma( v02, fma( v02, two_fifteenths, one_third ), one );
    v04  = v03 * rcp( fma( v03 * v03, v01, one ) );
    v04 += v04;
    v00  = fma( fms(  uy, cbz,  uz * cby ), v03, ux );
    v01  = fma( fms(  uz, cbx,  ux * cbz ), v03, uy );
    v02  = fma( fms(  ux, cby,  uy * cbx ), v03, uz );
    ux   = fma( fms( v01, cbz, v02 * cby ), v04, ux );
    uy   = fma( fms( v02, cbx, v00 * cbz ), v04, uy );
    uz   = fma( fms( v00, cby, v01 * cbx ), v04, uz );
    ux  += hax;
    uy  += hay;
    uz  += haz;

    //--------------------------------------------------------------------------
    // Store particle data.  Could use store_8x3_tr.
    //--------------------------------------------------------------------------
    store_8x4_tr( ux, uy, uz, q,
		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		  &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux );
  }
}

#if 0
void
uncenter_p_pipeline_v8( center_p_pipeline_args_t * args,
                        int pipeline_rank,
                        int n_pipeline )
{
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  particle_t           * ALIGNED(128) p;

  const float          * ALIGNED(32)  vp00;
  const float          * ALIGNED(32)  vp01;
  const float          * ALIGNED(32)  vp02;
  const float          * ALIGNED(32)  vp03;
  const float          * ALIGNED(32)  vp04;
  const float          * ALIGNED(32)  vp05;
  const float          * ALIGNED(32)  vp06;
  const float          * ALIGNED(32)  vp07;

  // Basic constants.
  const v8float qdt_2mc(    -args->qdt_2mc); // For backward half advance.
  const v8float qdt_4mc(-0.5*args->qdt_2mc); // For backward half Boris rotate.
  const v8float one(1.0);
  const v8float one_third(1.0/3.0);
  const v8float two_fifteenths(2.0/15.0);

  v8float dx, dy, dz, ux, uy, uz, q;
  v8float hax, hay, haz, cbx, cby, cbz;
  v8float v00, v01, v02, v03, v04, v05;
  v8int   ii;

  int first, nq;

  // Determine which particle quads this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, first, nq );

  p = args->p0 + first;

  nq >>= 3;

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=8 )
  {
    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------
    load_8x8_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		 &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
		 dx, dy, dz, ii, ux, uy, uz, q );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( const float * ALIGNED(32) ) ( f0 + ii(0) );
    vp01 = ( const float * ALIGNED(32) ) ( f0 + ii(1) );
    vp02 = ( const float * ALIGNED(32) ) ( f0 + ii(2) );
    vp03 = ( const float * ALIGNED(32) ) ( f0 + ii(3) );
    vp04 = ( const float * ALIGNED(32) ) ( f0 + ii(4) );
    vp05 = ( const float * ALIGNED(32) ) ( f0 + ii(5) );
    vp06 = ( const float * ALIGNED(32) ) ( f0 + ii(6) );
    vp07 = ( const float * ALIGNED(32) ) ( f0 + ii(7) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_8x8_tr( vp00, vp01, vp02, vp03,
		 vp04, vp05, vp06, vp07,
		 hax, v00, v01, v02, hay, v03, v04, v05 );

    hax = qdt_2mc*fma( fma( dy, v02, v01 ), dz, fma( dy, v00, hax ) );

    hay = qdt_2mc*fma( fma( dz, v05, v04 ), dx, fma( dz, v03, hay ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_8x8_tr( vp00+8, vp01+8, vp02+8, vp03+8,
		 vp04+8, vp05+8, vp06+8, vp07+8,
		 haz, v00, v01, v02, cbx, v03, cby, v04 );

    haz = qdt_2mc*fma( fma( dx, v02, v01 ), dy, fma( dx, v00, haz ) );

    cbx = fma( v03, dx, cbx );

    cby = fma( v04, dy, cby );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles, final.
    //--------------------------------------------------------------------------
    load_8x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
		 vp04+16, vp05+16, vp06+16, vp07+16,
		 cbz, v05 );

    cbz = fma( v05, dz, cbz );

    //--------------------------------------------------------------------------
    // Update momentum.
    //--------------------------------------------------------------------------
    v00  = qdt_4mc * rsqrt( one + fma( ux, ux, fma( uy, uy, uz * uz ) ) );
    v01  = fma( cbx, cbx, fma( cby, cby, cbz * cbz ) );
    v02  = ( v00 * v00 ) * v01;
    v03  = v00 * fma( v02, fma( v02, two_fifteenths, one_third ), one );
    v04  = v03 * rcp( fma( v03 * v03, v01, one ) );
    v04 += v04;

    v00  = fma( fms( uy, cbz, uz * cby ), v03, ux );
    v01  = fma( fms( uz, cbx, ux * cbz ), v03, uy );
    v02  = fma( fms( ux, cby, uy * cbx ), v03, uz );

    ux   = fma( fms( v01, cbz, v02 * cby ), v04, ux );
    uy   = fma( fms( v02, cbx, v00 * cbz ), v04, uy );
    uz   = fma( fms( v00, cby, v01 * cbx ), v04, uz );

    ux  += hax;
    uy  += hay;
    uz  += haz;

    //--------------------------------------------------------------------------
    // Store particle data.  Could use store_8x3_tr.
    //--------------------------------------------------------------------------
    store_8x4_tr( ux, uy, uz, q,
		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		  &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux );
  }
}
#endif
