#define IN_spa

#include "spa_private.h"

#if defined(V16_ACCELERATION)

using namespace v16;

void
hydro_p_pipeline_v16( hydro_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline)
{
  const species_t      *              sp = args->sp;
  const grid_t         *              g  = sp->g;
  /**/  hydro_t        * ALIGNED(128) h  = args->h + pipeline_rank*args->h_size;
  const particle_t     * ALIGNED(128) p  = sp->p;
  const interpolator_t * ALIGNED(128) f  = args->f;
  const bool              charge_weight  = args->charge_weight;

  /**/  float          * ALIGNED(64)  vp00;
  /**/  float          * ALIGNED(64)  vp01;
  /**/  float          * ALIGNED(64)  vp02;
  /**/  float          * ALIGNED(64)  vp03;
  /**/  float          * ALIGNED(64)  vp04;
  /**/  float          * ALIGNED(64)  vp05;
  /**/  float          * ALIGNED(64)  vp06;
  /**/  float          * ALIGNED(64)  vp07;
  /**/  float          * ALIGNED(64)  vp08;
  /**/  float          * ALIGNED(64)  vp09;
  /**/  float          * ALIGNED(64)  vp10;
  /**/  float          * ALIGNED(64)  vp11;
  /**/  float          * ALIGNED(64)  vp12;
  /**/  float          * ALIGNED(64)  vp13;
  /**/  float          * ALIGNED(64)  vp14;
  /**/  float          * ALIGNED(64)  vp15;

  const v16float qsp(charge_weight ? sp->q : sp->m);
  const v16float qdt_2mc(args->qdt_2mc);
  const v16float qdt_4mc2(args->qdt_2mc / (2*g->cvac));
  const v16float mspc(args->msp*g->cvac);
  const v16float one(1.);
  const v16float zero(0.);
  const v16float c(sp->g->cvac);
  const v16float r8V(g->r8V);
  const v16float two_fifteenths(2./15.);
  const v16float one_third(1./3.);

  const v16int stride_10(g->sx);                   // (i,j,k) -> (i+1, j,   k  )
  const v16int stride_21(g->sy - g->sx);           // (i,j,k) -> (i-1, j+1, k  )
  const v16int stride_43(g->sz - g->sy - g->sx);   // (i,j,k) -> (i-1, j-1, k+1)

  v16float dx, dy, dz, ux, uy, uz, w, vx, vy, vz, ke_mc;
  v16float hax, hay, haz, cbx, cby, cbz;
  v16float v00, v01, v02, v03, v04, v05, v06, v07;
  v16float v08, v09, v10, v11, v12, v13, v14, v15;
  v16float t, w0, w1, w2, w3, w4, w5, w6, w7;
  v16int ii;

  int n0, nq;

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, nq );
  p += n0;
  nq >>= 4;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=16 )
  {
    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------
    load_16x8_tr_p( &p[ 0].dx, &p[ 2].dx, &p[ 4].dx, &p[ 6].dx,
                    &p[ 8].dx, &p[10].dx, &p[12].dx, &p[14].dx,
		                dx, dy, dz, ii, ux, uy, uz, w );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( float * ALIGNED(64) ) ( f + ii( 0) );
    vp01 = ( float * ALIGNED(64) ) ( f + ii( 1) );
    vp02 = ( float * ALIGNED(64) ) ( f + ii( 2) );
    vp03 = ( float * ALIGNED(64) ) ( f + ii( 3) );
    vp04 = ( float * ALIGNED(64) ) ( f + ii( 4) );
    vp05 = ( float * ALIGNED(64) ) ( f + ii( 5) );
    vp06 = ( float * ALIGNED(64) ) ( f + ii( 6) );
    vp07 = ( float * ALIGNED(64) ) ( f + ii( 7) );
    vp08 = ( float * ALIGNED(64) ) ( f + ii( 8) );
    vp09 = ( float * ALIGNED(64) ) ( f + ii( 9) );
    vp10 = ( float * ALIGNED(64) ) ( f + ii(10) );
    vp11 = ( float * ALIGNED(64) ) ( f + ii(11) );
    vp12 = ( float * ALIGNED(64) ) ( f + ii(12) );
    vp13 = ( float * ALIGNED(64) ) ( f + ii(13) );
    vp14 = ( float * ALIGNED(64) ) ( f + ii(14) );
    vp15 = ( float * ALIGNED(64) ) ( f + ii(15) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles and half advance with E
    //--------------------------------------------------------------------------
    load_16x16_tr( vp00, vp01, vp02, vp03,
                   vp04, vp05, vp06, vp07,
                   vp08, vp09, vp10, vp11,
                   vp12, vp13, vp14, vp15,
                   hax, v00, v01, v02, hay, v03, v04, v05,
                   haz, v06, v07, v08, cbx, v09, cby, v10 );

    hax = qdt_2mc*fma( fma( v02, dy, v01 ), dz, fma( v00, dy, hax ) );
    hay = qdt_2mc*fma( fma( v05, dz, v04 ), dx, fma( v03, dz, hay ) );
    haz = qdt_2mc*fma( fma( v08, dx, v07 ), dy, fma( v06, dx, haz ) );
    cbx = fma( v09, dx, cbx );
    cby = fma( v10, dy, cby );

    load_16x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
                  vp04+16, vp05+16, vp06+16, vp07+16,
                  vp08+16, vp09+16, vp10+16, vp11+16,
                  vp12+16, vp13+16, vp14+16, vp15+16,
                  cbz, v05 );

    cbz = fma( v05, dz, cbz );

    //--------------------------------------------------------------------------
    // Update momentum.
    //--------------------------------------------------------------------------

    ux  += hax;
    uy  += hay;
    uz  += haz;

    //--------------------------------------------------------------------------
    // Compute kinetic energy.
    //--------------------------------------------------------------------------

    ke_mc = fma( ux, ux, fma( uy, uy, uz*uz ) );
    vz    = sqrt( one + ke_mc );
    ke_mc = (c*ke_mc)*rcp( one + vz );
    vz    = c*rcp( vz );

    //--------------------------------------------------------------------------
    // Half Boris Advance.
    //--------------------------------------------------------------------------

    v00  = qdt_4mc2*vz;
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

    //--------------------------------------------------------------------------
    // Compute velocity
    //--------------------------------------------------------------------------

    vx = ux*vz;
    vy = uy*vz;
    vz = uz*vz;

    //--------------------------------------------------------------------------
    // Trilinear coefficients
    //--------------------------------------------------------------------------

    w0  = r8V*w;    // w0 = (1/8)(w/V)
    dx *= w0;       // dx = (1/8)(w/V) x
    w1  = w0+dx;    // w1 = (1/8)(w/V) + (1/8)(w/V)x = (1/8)(w/V)(1+x)
    w0 -= dx;       // w0 = (1/8)(w/V) - (1/8)(w/V)x = (1/8)(w/V)(1-x)
    w3  = 1+dy;     // w3 = 1+y
    w2  = w0*w3;    // w2 = (1/8)(w/V)(1-x)(1+y)
    w3 *= w1;       // w3 = (1/8)(w/V)(1+x)(1+y)
    dy  = 1-dy;     // dy = 1-y
    w0 *= dy;       // w0 = (1/8)(w/V)(1-x)(1-y)
    w1 *= dy;       // w1 = (1/8)(w/V)(1+x)(1-y)
    w7  = 1+dz;     // w7 = 1+z
    w4  = w0*w7;    // w4 = (1/8)(w/V)(1-x)(1-y)(1+z) = (w/V) trilin_0 *Done
    w5  = w1*w7;    // w5 = (1/8)(w/V)(1+x)(1-y)(1+z) = (w/V) trilin_1 *Done
    w6  = w2*w7;    // w6 = (1/8)(w/V)(1-x)(1+y)(1+z) = (w/V) trilin_2 *Done
    w7 *= w3;       // w7 = (1/8)(w/V)(1+x)(1+y)(1+z) = (w/V) trilin_3 *Done
    dz  = 1-dz;     // dz = 1-z
    w0 *= dz;       // w0 = (1/8)(w/V)(1-x)(1-y)(1-z) = (w/V) trilin_4 *Done
    w1 *= dz;       // w1 = (1/8)(w/V)(1+x)(1-y)(1-z) = (w/V) trilin_5 *Done
    w2 *= dz;       // w2 = (1/8)(w/V)(1-x)(1+y)(1-z) = (w/V) trilin_6 *Done
    w3 *= dz;       // w3 = (1/8)(w/V)(1+x)(1+y)(1-z) = (w/V) trilin_7 *Done

    //--------------------------------------------------------------------------
    // Accumulation.
    //--------------------------------------------------------------------------

#   define ACCUM_HYDRO(w)                                   \
    vp00 = ( float * ALIGNED(64) ) ( h + ii( 0) );          \
    vp01 = ( float * ALIGNED(64) ) ( h + ii( 1) );          \
    vp02 = ( float * ALIGNED(64) ) ( h + ii( 2) );          \
    vp03 = ( float * ALIGNED(64) ) ( h + ii( 3) );          \
    vp04 = ( float * ALIGNED(64) ) ( h + ii( 4) );          \
    vp05 = ( float * ALIGNED(64) ) ( h + ii( 5) );          \
    vp06 = ( float * ALIGNED(64) ) ( h + ii( 6) );          \
    vp07 = ( float * ALIGNED(64) ) ( h + ii( 7) );          \
    vp08 = ( float * ALIGNED(64) ) ( h + ii( 8) );          \
    vp09 = ( float * ALIGNED(64) ) ( h + ii( 9) );          \
    vp10 = ( float * ALIGNED(64) ) ( h + ii(10) );          \
    vp11 = ( float * ALIGNED(64) ) ( h + ii(11) );          \
    vp12 = ( float * ALIGNED(64) ) ( h + ii(12) );          \
    vp13 = ( float * ALIGNED(64) ) ( h + ii(13) );          \
    vp14 = ( float * ALIGNED(64) ) ( h + ii(14) );          \
    vp15 = ( float * ALIGNED(64) ) ( h + ii(15) );          \
    t   = qsp*w;                                            \
    v00 = t*vx;       /* w vx */                            \
    v01 = t*vy;       /* w vy */                            \
    v02 = t*vz;       /* w vz */                            \
    v03 = t;          /* w */                               \
    t   = mspc*w;                                           \
    dx  = t*ux;                                             \
    dy  = t*uy;                                             \
    dz  = t*uz;                                             \
    v04 = dx;         /* m c w ux */                        \
    v05 = dy;         /* m c w uy */                        \
    v06 = dz;         /* m c w uz */                        \
    v07 = t*ke_mc;    /* m c w (gamma-1) */                 \
    v08 = dx*vx;      /* m c w ux vx */                     \
    v09 = dy*vy;      /* m c w uy vy */                     \
    v10 = dz*vz;      /* m c w uz vz */                     \
    v11 = dy*vz;      /* m c w uy vz */                     \
    v12 = dz*vx;      /* m c w uz vx */                     \
    v13 = dx*vy;      /* m c w ux vy */                     \
    v14 = zero;       /* pad[0] */                          \
    v15 = zero;       /* pad[1] */                          \
    transpose( v00, v01, v02, v03, v04, v05, v06, v07,	    \
               v08, v09, v10, v11, v12, v13, v14, v15);     \
    increment_16x1( vp00, v00 );                            \
    increment_16x1( vp01, v01 );                            \
    increment_16x1( vp02, v02 );                            \
    increment_16x1( vp03, v03 );                            \
    increment_16x1( vp04, v04 );                            \
    increment_16x1( vp05, v05 );                            \
    increment_16x1( vp06, v06 );                            \
    increment_16x1( vp07, v07 );                            \
    increment_16x1( vp08, v08 );                            \
    increment_16x1( vp09, v09 );                            \
    increment_16x1( vp10, v10 );                            \
    increment_16x1( vp11, v11 );                            \
    increment_16x1( vp12, v12 );                            \
    increment_16x1( vp13, v13 );                            \
    increment_16x1( vp14, v14 );                            \
    increment_16x1( vp15, v15 )

    /**/             ACCUM_HYDRO(w0); // Cell i,j,k
    ii += stride_10; ACCUM_HYDRO(w1); // Cell i+1,j,k
    ii += stride_21; ACCUM_HYDRO(w2); // Cell i,j+1,k
    ii += stride_10; ACCUM_HYDRO(w3); // Cell i+1,j+1,k
    ii += stride_43; ACCUM_HYDRO(w4); // Cell i,j,k+1
    ii += stride_10; ACCUM_HYDRO(w5); // Cell i+1,j,k+1
    ii += stride_21; ACCUM_HYDRO(w6); // Cell i,j+1,k+1
    ii += stride_10; ACCUM_HYDRO(w7); // Cell i+1,j+1,k+1

#   undef ACCUM_HYDRO

  }

}

#else

void
hydro_p_pipeline_v16( hydro_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline) {
  ERROR(("No hydro_p_pipeline_v16 implementation."));
}

#endif
