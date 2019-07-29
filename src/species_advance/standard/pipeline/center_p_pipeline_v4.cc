#define IN_spa

#include "spa_private.h"

#if defined(V4_ACCELERATION)

using namespace v4;

#ifdef V4_NEON_ACCELERATION_SNOUT

void
center_p_pipeline_v4( center_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline )
{
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  particle_t           * ALIGNED(128) p;

  const float          * ALIGNED(16)  vp00;
  const float          * ALIGNED(16)  vp01;
  const float          * ALIGNED(16)  vp02;
  const float          * ALIGNED(16)  vp03;

  const v4float qdt_2mc(    args->qdt_2mc);
  const v4float qdt_4mc(0.5*args->qdt_2mc); // For half Boris rotate.
  const v4float one(1.0);
  const v4float one_third(1.0/3.0);
  const v4float two_fifteenths(2.0/15.0);

  v4float dx, dy, dz, ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v00, v01, v02, v03, v04, v05;
  v4float v06, v07, v08, v09, v10;
  v4int   ii;

  int itmp, nq;

  // Determine which particle blocks this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 2;

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=4 )
  {
    //--------------------------------------------------------------------------
    // Load particle position data.
    //--------------------------------------------------------------------------
    load_4x8_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
                 dx, dy, dz, ii, ux, uy, uz, q );

    // load_4x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
    // 		 dx, dy, dz, ii );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( const float * ALIGNED(16) ) ( f0 + ii(0) );
    vp01 = ( const float * ALIGNED(16) ) ( f0 + ii(1) );
    vp02 = ( const float * ALIGNED(16) ) ( f0 + ii(2) );
    vp03 = ( const float * ALIGNED(16) ) ( f0 + ii(3) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x16_tr( vp00, vp01, vp02, vp03,
                  hax, v00, v01, v02,
                  hay, v03, v04, v05,
                  haz, v06, v07, v08,
                  cbx, v09, cby, v10 );

    // load_4x4_tr( vp00, vp01, vp02, vp03,
    // 		 hax, v00, v01, v02 );

    hax = qdt_2mc * fma( fma( dy, v02, v01 ), dz, fma( dy, v00, hax ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    // load_4x4_tr( vp00+4, vp01+4, vp02+4, vp03+4,
    // 		 hay, v03, v04, v05 );

    hay = qdt_2mc * fma( fma( dz, v05, v04 ), dx, fma( dz, v03, hay ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    // load_4x4_tr( vp00+8, vp01+8, vp02+8, vp03+8,
    // 		 haz, v00, v01, v02 );

    haz = qdt_2mc * fma( fma( dx, v08, v07 ), dy, fma( dx, v06, haz ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    // load_4x4_tr( vp00+12, vp01+12, vp02+12, vp03+12,
    // 		 cbx, v03, cby, v04 );

    cbx = fma( v09, dx, cbx );
    cby = fma( v10, dy, cby );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles, final.
    //--------------------------------------------------------------------------
    load_4x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
		 cbz, v05 );

    cbz = fma( v05, dz, cbz );

    //--------------------------------------------------------------------------
    // Load particle momentum data.  Could use load_4x3_tr.
    //--------------------------------------------------------------------------
    // load_4x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
    // 		 ux, uy, uz, q );

    //--------------------------------------------------------------------------
    // Update momentum.
    //--------------------------------------------------------------------------
    ux  += hax;
    uy  += hay;
    uz  += haz;

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

    //--------------------------------------------------------------------------
    // Store particle momentum data.  Could use store_4x3_tr.
    //--------------------------------------------------------------------------
    store_4x4_tr( ux, uy, uz, q,
		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux );
  }
}

#else

void
center_p_pipeline_v4( center_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline )
{
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  particle_t           * ALIGNED(128) p;

  const float          * ALIGNED(16)  vp00;
  const float          * ALIGNED(16)  vp01;
  const float          * ALIGNED(16)  vp02;
  const float          * ALIGNED(16)  vp03;

  const v4float qdt_2mc(    args->qdt_2mc);
  const v4float qdt_4mc(0.5*args->qdt_2mc); // For half Boris rotate.
  const v4float one(1.0);
  const v4float one_third(1.0/3.0);
  const v4float two_fifteenths(2.0/15.0);

  v4float dx, dy, dz, ux, uy, uz, q;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v00, v01, v02, v03, v04, v05;
  v4int   ii;

  int itmp, nq;

  // Determine which particle blocks this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, nq );

  p = args->p0 + itmp;

  nq >>= 2;

  // Process the particle blocks for this pipeline.

  for( ; nq; nq--, p+=4 )
  {
    //--------------------------------------------------------------------------
    // Load particle position data.
    //--------------------------------------------------------------------------
    load_4x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
		 dx, dy, dz, ii );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( const float * ALIGNED(16) ) ( f0 + ii(0) );
    vp01 = ( const float * ALIGNED(16) ) ( f0 + ii(1) );
    vp02 = ( const float * ALIGNED(16) ) ( f0 + ii(2) );
    vp03 = ( const float * ALIGNED(16) ) ( f0 + ii(3) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00, vp01, vp02, vp03,
		 hax, v00, v01, v02 );

    hax = qdt_2mc*fma( fma( dy, v02, v01 ), dz, fma( dy, v00, hax ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00+4, vp01+4, vp02+4, vp03+4,
		 hay, v03, v04, v05 );

    hay = qdt_2mc*fma( fma( dz, v05, v04 ), dx, fma( dz, v03, hay ) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00+8, vp01+8, vp02+8, vp03+8,
		 haz, v00, v01, v02 );

    haz = qdt_2mc*fma( fma( dx, v02, v01 ), dy, fma( dx, v00, haz ) );

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
    // Load particle momentum data.  Could use load_4x3_tr.
    //--------------------------------------------------------------------------
    load_4x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux,
		 ux, uy, uz, q );

    //--------------------------------------------------------------------------
    // Update momentum.
    //--------------------------------------------------------------------------
    ux  += hax;
    uy  += hay;
    uz  += haz;

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

    //--------------------------------------------------------------------------
    // Store particle momentum data.  Could use store_4x3_tr.
    //--------------------------------------------------------------------------
    store_4x4_tr( ux, uy, uz, q,
		  &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux );
  }
}

#endif

#else

void
center_p_pipeline_v4( center_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline )
{
  // No v4 implementation.
  ERROR( ( "No center_p_pipeline_v4 implementation." ) );
}

#endif
