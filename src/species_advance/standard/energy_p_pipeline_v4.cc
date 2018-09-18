using namespace v4;

void
energy_p_pipeline_v4( energy_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline )
{
  const interpolator_t * RESTRICT ALIGNED(128) f = args->f;
  const particle_t     * RESTRICT ALIGNED(128) p = args->p;

  const float          * RESTRICT ALIGNED(16)  vp00;
  const float          * RESTRICT ALIGNED(16)  vp01;
  const float          * RESTRICT ALIGNED(16)  vp02;
  const float          * RESTRICT ALIGNED(16)  vp03;

  const v4float qdt_2mc(args->qdt_2mc);
  const v4float msp(args->msp);
  const v4float one(1.0);

  v4float dx, dy, dz;
  v4float ex, ey, ez;
  v4float v00, v01, v02, w;
  v4int i;

  double en00 = 0.0, en01 = 0.0, en02 = 0.0, en03 = 0.0;

  int n0, nq;

  // Determine which particle blocks this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, nq );

  p += n0;

  nq >>= 2;

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=4 )
  {
    //--------------------------------------------------------------------------
    // Load particle position data.
    //--------------------------------------------------------------------------
    load_4x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
                 dx, dy, dz, i );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( float * ) ( f + i(0) );
    vp01 = ( float * ) ( f + i(1) );
    vp02 = ( float * ) ( f + i(2) );
    vp03 = ( float * ) ( f + i(3) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00, vp01, vp02, vp03,
                 ex, v00, v01, v02 );

    ex = fma( fma( dy, v02, v01 ), dz, fma( dy, v00, ex ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00+4, vp01+4, vp02+4, vp03+4,
                 ey, v00, v01, v02);

    ey = fma( fma( dz, v02, v01 ), dx, fma( dz, v00, ey ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00+8, vp01+8, vp02+8, vp03+8,
                 ez, v00, v01, v02);

    ez = fma( fma( dx, v02, v01 ), dy, fma( dx, v00, ez ) );

    //--------------------------------------------------------------------------
    // Load particle momentum data.
    //--------------------------------------------------------------------------
    load_4x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
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
    en00 += ( double ) v00(0);
    en01 += ( double ) v00(1);
    en02 += ( double ) v00(2);
    en03 += ( double ) v00(3);
  }

  //--------------------------------------------------------------------------
  // Accumulate energy for each rank or thread.
  //--------------------------------------------------------------------------
  args->en[pipeline_rank] = en00 + en01 + en02 + en03;
}
