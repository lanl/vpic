using namespace v8;

void
energy_p_pipeline_v8( energy_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline )
{
  const interpolator_t * RESTRICT ALIGNED(128) f = args->f;
  const particle_t     * RESTRICT ALIGNED(128) p = args->p;

  const float          * RESTRICT ALIGNED(16)  vp0;
  const float          * RESTRICT ALIGNED(16)  vp1;
  const float          * RESTRICT ALIGNED(16)  vp2;
  const float          * RESTRICT ALIGNED(16)  vp3;
  const float          * RESTRICT ALIGNED(16)  vp4;
  const float          * RESTRICT ALIGNED(16)  vp5;
  const float          * RESTRICT ALIGNED(16)  vp6;
  const float          * RESTRICT ALIGNED(16)  vp7;

  const v8float qdt_2mc(args->qdt_2mc);
  const v8float msp(args->msp);
  const v8float one(1.);

  v8float dx, dy, dz;
  v8float ex, ey, ez;
  v8float v0, v1, v2, w;
  v8int i;

  double en0 = 0, en1 = 0, en2 = 0, en3 = 0;
  double en4 = 0, en5 = 0, en6 = 0, en7 = 0;

  int n0, nq;

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, nq );
  p += n0;
  nq >>= 3;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=8 )
  {
    load_8x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		 &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
		 dx, dy, dz, i );

    // Interpolate fields

    vp0 = (float *)(f + i(0));
    vp1 = (float *)(f + i(1));
    vp2 = (float *)(f + i(2));
    vp3 = (float *)(f + i(3));
    vp4 = (float *)(f + i(4));
    vp5 = (float *)(f + i(5));
    vp6 = (float *)(f + i(6));
    vp7 = (float *)(f + i(7));

    load_8x4_tr( vp0, vp1, vp2, vp3,
		 vp4, vp5, vp6, vp7,
		 ex, v0, v1, v2 );

    ex = fma( fma( dy, v2, v1 ), dz, fma( dy, v0, ex ) );

    load_8x4_tr( vp0+4, vp1+4, vp2+4, vp3+4,
		 vp4+4, vp5+4, vp6+4, vp7+4,
		 ey, v0, v1, v2 );

    ey = fma( fma( dz, v2, v1 ), dx, fma( dz, v0, ey ) );

    load_8x4_tr( vp0+8, vp1+8, vp2+8, vp3+8,
		 vp4+8, vp5+8, vp6+8, vp7+8,
		 ez, v0, v1, v2 );

    ez = fma( fma( dx, v2, v1 ), dy, fma( dx, v0, ez ) );

    // Update momentum to half step
    // (note Boris rotation does not change energy so it is unnecessary)

    load_8x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		 &p[4].ux, &p[5].ux, &p[6].ux, &p[7].ux,
		 v0, v1, v2, w );

    v0  = fma( ex, qdt_2mc, v0 );
    v1  = fma( ey, qdt_2mc, v1 );
    v2  = fma( ez, qdt_2mc, v2 );

    // Accumulate energy

    v0 = fma( v0,v0, fma( v1,v1, v2*v2 ) );
    v0 = (msp * w) * (v0 / (one + sqrt(one + v0))); 

    en0 += (double)v0(0);
    en1 += (double)v0(1);
    en2 += (double)v0(2);
    en3 += (double)v0(3);
    en4 += (double)v0(4);
    en5 += (double)v0(5);
    en6 += (double)v0(6);
    en7 += (double)v0(7);
  }

  args->en[pipeline_rank] = en0 + en1 + en2 + en3 +
                            en4 + en5 + en6 + en7;
}
