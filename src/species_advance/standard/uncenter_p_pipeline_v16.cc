#define IN_spa

#include "spa_private.h"

#if defined(V16_ACCELERATION)

using namespace v16;

void
uncenter_p_pipeline_v16( center_p_pipeline_args_t * args,
                         int pipeline_rank,
                         int n_pipeline )
{
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  particle_t           * ALIGNED(128) p;

  const float          * ALIGNED(64)  vp00;
  const float          * ALIGNED(64)  vp01;
  const float          * ALIGNED(64)  vp02;
  const float          * ALIGNED(64)  vp03;
  const float          * ALIGNED(64)  vp04;
  const float          * ALIGNED(64)  vp05;
  const float          * ALIGNED(64)  vp06;
  const float          * ALIGNED(64)  vp07;
  const float          * ALIGNED(64)  vp08;
  const float          * ALIGNED(64)  vp09;
  const float          * ALIGNED(64)  vp10;
  const float          * ALIGNED(64)  vp11;
  const float          * ALIGNED(64)  vp12;
  const float          * ALIGNED(64)  vp13;
  const float          * ALIGNED(64)  vp14;
  const float          * ALIGNED(64)  vp15;

  const v16float qdt_2mc(    -args->qdt_2mc); // For backward half advance.
  const v16float qdt_4mc(-0.5*args->qdt_2mc); // For backward half Boris rotate.
  const v16float one(1.0);
  const v16float one_third(1.0/3.0);
  const v16float two_fifteenths(2.0/15.0);

  v16float dx, dy, dz, ux, uy, uz, q;
  v16float hax, hay, haz, cbx, cby, cbz;
  v16float v00, v01, v02, v03, v04, v05, v06, v07, v08, v09, v10;
  v16int   ii;

  int first, nq;

  // Determine which particle blocks this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, first, nq );

  p = args->p0 + first;

  nq >>= 4;

  // Process the particle blocks for this pipeline.

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
    vp00 = ( const float * ALIGNED(64) ) ( f0 + ii( 0) );
    vp01 = ( const float * ALIGNED(64) ) ( f0 + ii( 1) );
    vp02 = ( const float * ALIGNED(64) ) ( f0 + ii( 2) );
    vp03 = ( const float * ALIGNED(64) ) ( f0 + ii( 3) );
    vp04 = ( const float * ALIGNED(64) ) ( f0 + ii( 4) );
    vp05 = ( const float * ALIGNED(64) ) ( f0 + ii( 5) );
    vp06 = ( const float * ALIGNED(64) ) ( f0 + ii( 6) );
    vp07 = ( const float * ALIGNED(64) ) ( f0 + ii( 7) );
    vp08 = ( const float * ALIGNED(64) ) ( f0 + ii( 8) );
    vp09 = ( const float * ALIGNED(64) ) ( f0 + ii( 9) );
    vp10 = ( const float * ALIGNED(64) ) ( f0 + ii(10) );
    vp11 = ( const float * ALIGNED(64) ) ( f0 + ii(11) );
    vp12 = ( const float * ALIGNED(64) ) ( f0 + ii(12) );
    vp13 = ( const float * ALIGNED(64) ) ( f0 + ii(13) );
    vp14 = ( const float * ALIGNED(64) ) ( f0 + ii(14) );
    vp15 = ( const float * ALIGNED(64) ) ( f0 + ii(15) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00, vp01, vp02, vp03,
                   vp04, vp05, vp06, vp07,
                   vp08, vp09, vp10, vp11,
                   vp12, vp13, vp14, vp15,
                   hax, v00, v01, v02, hay, v03, v04, v05,
                   haz, v06, v07, v08, cbx, v09, cby, v10 );

    hax = qdt_2mc*fma( fma( dy, v02, v01 ), dz, fma( dy, v00, hax ) );

    hay = qdt_2mc*fma( fma( dz, v05, v04 ), dx, fma( dz, v03, hay ) );

    haz = qdt_2mc*fma( fma( dx, v08, v07 ), dy, fma( dx, v06, haz ) );

    cbx = fma( v09, dx, cbx );

    cby = fma( v10, dy, cby );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles, final.
    //--------------------------------------------------------------------------
    load_16x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
                  vp04+16, vp05+16, vp06+16, vp07+16,
                  vp08+16, vp09+16, vp10+16, vp11+16,
                  vp12+16, vp13+16, vp14+16, vp15+16,
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
    // Store particle momentum data.  Could use store_16x4_tr_p or
    // store_16x3_tr_p.
    //--------------------------------------------------------------------------
    store_16x8_tr_p( dx, dy, dz, ii, ux, uy, uz, q,
                     &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
                     &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx );
  }
}

#else

void
uncenter_p_pipeline_v16( center_p_pipeline_args_t * args,
                         int pipeline_rank,
                         int n_pipeline )
{
  // No v16 implementation.
  ERROR( ( "No uncenter_p_pipeline_v16 implementation." ) );
}

#endif
