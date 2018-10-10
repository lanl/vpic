using namespace v16;

void
energy_p_pipeline_v16( energy_p_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline )
{
  const interpolator_t * RESTRICT ALIGNED(128) f = args->f;
  const particle_t     * RESTRICT ALIGNED(128) p = args->p;

  const float          * RESTRICT ALIGNED(64)  vp00;
  const float          * RESTRICT ALIGNED(64)  vp01;
  const float          * RESTRICT ALIGNED(64)  vp02;
  const float          * RESTRICT ALIGNED(64)  vp03;
  const float          * RESTRICT ALIGNED(64)  vp04;
  const float          * RESTRICT ALIGNED(64)  vp05;
  const float          * RESTRICT ALIGNED(64)  vp06;
  const float          * RESTRICT ALIGNED(64)  vp07;
  const float          * RESTRICT ALIGNED(64)  vp08;
  const float          * RESTRICT ALIGNED(64)  vp09;
  const float          * RESTRICT ALIGNED(64)  vp10;
  const float          * RESTRICT ALIGNED(64)  vp11;
  const float          * RESTRICT ALIGNED(64)  vp12;
  const float          * RESTRICT ALIGNED(64)  vp13;
  const float          * RESTRICT ALIGNED(64)  vp14;
  const float          * RESTRICT ALIGNED(64)  vp15;

  const v16float qdt_2mc(args->qdt_2mc);
  const v16float msp(args->msp);
  const v16float one(1.0);

  v16float dx, dy, dz;
  v16float ex, ey, ez;
  v16float v00, v01, v02, w;
  v16int i;

  double en00 = 0.0, en01 = 0.0, en02 = 0.0, en03 = 0.0;
  double en04 = 0.0, en05 = 0.0, en06 = 0.0, en07 = 0.0;
  double en08 = 0.0, en09 = 0.0, en10 = 0.0, en11 = 0.0;
  double en12 = 0.0, en13 = 0.0, en14 = 0.0, en15 = 0.0;

  int n0, nq;

  // Determine which particle blocks this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, nq );

  p += n0;

  nq >>= 4;

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=16 )
  {
    //--------------------------------------------------------------------------
    // Load particle position data.
    //--------------------------------------------------------------------------
    load_16x4_tr( &p[ 0].dx, &p[ 1].dx, &p[ 2].dx, &p[ 3].dx,
                  &p[ 4].dx, &p[ 5].dx, &p[ 6].dx, &p[ 7].dx,
                  &p[ 8].dx, &p[ 9].dx, &p[10].dx, &p[11].dx,
                  &p[12].dx, &p[13].dx, &p[14].dx, &p[15].dx,
                  dx, dy, dz, i );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( float * ) ( f + i( 0) );
    vp01 = ( float * ) ( f + i( 1) );
    vp02 = ( float * ) ( f + i( 2) );
    vp03 = ( float * ) ( f + i( 3) );
    vp04 = ( float * ) ( f + i( 4) );
    vp05 = ( float * ) ( f + i( 5) );
    vp06 = ( float * ) ( f + i( 6) );
    vp07 = ( float * ) ( f + i( 7) );
    vp08 = ( float * ) ( f + i( 8) );
    vp09 = ( float * ) ( f + i( 9) );
    vp10 = ( float * ) ( f + i(10) );
    vp11 = ( float * ) ( f + i(11) );
    vp12 = ( float * ) ( f + i(12) );
    vp13 = ( float * ) ( f + i(13) );
    vp14 = ( float * ) ( f + i(14) );
    vp15 = ( float * ) ( f + i(15) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_16x4_tr( vp00, vp01, vp02, vp03,
                  vp04, vp05, vp06, vp07,
                  vp08, vp09, vp10, vp11,
                  vp12, vp13, vp14, vp15,
                  ex, v00, v01, v02 );

    ex = fma( fma( dy, v02, v01 ), dz, fma( dy, v00, ex ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_16x4_tr( vp00+4, vp01+4, vp02+4, vp03+4,
                  vp04+4, vp05+4, vp06+4, vp07+4,
                  vp08+4, vp09+4, vp10+4, vp11+4,
                  vp12+4, vp13+4, vp14+4, vp15+4,
                  ey, v00, v01, v02 );

    ey = fma( fma( dz, v02, v01 ), dx, fma( dz, v00, ey ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_16x4_tr( vp00+8, vp01+8, vp02+8, vp03+8,
                  vp04+8, vp05+8, vp06+8, vp07+8,
                  vp08+8, vp09+8, vp10+8, vp11+8,
                  vp12+8, vp13+8, vp14+8, vp15+8,
                  ez, v00, v01, v02 );

    ez = fma( fma( dx, v02, v01 ), dy, fma( dx, v00, ez ) );

    //--------------------------------------------------------------------------
    // Load particle momentum data.
    //--------------------------------------------------------------------------
    load_16x4_tr( &p[ 0].ux, &p[ 1].ux, &p[ 2].ux, &p[ 3].ux,
                  &p[ 4].ux, &p[ 5].ux, &p[ 6].ux, &p[ 7].ux,
                  &p[ 8].ux, &p[ 9].ux, &p[10].ux, &p[11].ux,
                  &p[12].ux, &p[13].ux, &p[14].ux, &p[15].ux,
                  v00, v01, v02, w );

    //--------------------------------------------------------------------------
    // Update momentum to half step. Note that Boris rotation does not change
    // energy and thus is not necessary.
    //--------------------------------------------------------------------------
    v00 = fma( ex, qdt_2mc, v00 );
    v01 = fma( ey, qdt_2mc, v01 );
    v02 = fma( ez, qdt_2mc, v02 );

    //--------------------------------------------------------------------------
    // Calculate kinetic energy of particles.
    //--------------------------------------------------------------------------
    v00 = fma( v00, v00, fma( v01, v01, v02 * v02 ) );

    v00 = ( msp * w ) * ( v00 / ( one + sqrt( one + v00 ) ) ); 

    //--------------------------------------------------------------------------
    // Accumulate energy for each vector element.
    //--------------------------------------------------------------------------
    en00 += ( double ) v00( 0);
    en01 += ( double ) v00( 1);
    en02 += ( double ) v00( 2);
    en03 += ( double ) v00( 3);
    en04 += ( double ) v00( 4);
    en05 += ( double ) v00( 5);
    en06 += ( double ) v00( 6);
    en07 += ( double ) v00( 7);
    en08 += ( double ) v00( 8);
    en09 += ( double ) v00( 9);
    en10 += ( double ) v00(10);
    en11 += ( double ) v00(11);
    en12 += ( double ) v00(12);
    en13 += ( double ) v00(13);
    en14 += ( double ) v00(14);
    en15 += ( double ) v00(15);
  }

  //--------------------------------------------------------------------------
  // Accumulate energy for each rank or thread.
  //--------------------------------------------------------------------------
  args->en[pipeline_rank] = en00 + en01 + en02 + en03 +
                            en04 + en05 + en06 + en07 +
                            en08 + en09 + en10 + en11 +
                            en12 + en13 + en14 + en15;
}
