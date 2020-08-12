#define IN_spa

#include "spa_private.h"

#if defined(V8_ACCELERATION)

using namespace v8;

void
hydro_p_pipeline_v8( hydro_p_pipeline_args_t * args,
                     int pipeline_rank,
                     int n_pipeline)
{
  const species_t      *              sp = args->sp;
  const grid_t         *              g  = sp->g;
  /**/  hydro_t        * ALIGNED(128) h  = args->h + pipeline_rank*args->h_size;
  const particle_t     * ALIGNED(128) p  = sp->p;
  const interpolator_t * ALIGNED(128) f  = args->f;
  const bool              charge_weight  = args->charge_weight;

  /**/  float          * ALIGNED(32)  vp00;
  /**/  float          * ALIGNED(32)  vp01;
  /**/  float          * ALIGNED(32)  vp02;
  /**/  float          * ALIGNED(32)  vp03;
  /**/  float          * ALIGNED(32)  vp04;
  /**/  float          * ALIGNED(32)  vp05;
  /**/  float          * ALIGNED(32)  vp06;
  /**/  float          * ALIGNED(32)  vp07;

  const v8float qsp(charge_weight ? sp->q : sp->m);
  const v8float qdt_2mc(args->qdt_2mc);
  const v8float qdt_4mc2(args->qdt_2mc / (2*g->cvac));
  const v8float mspc(args->msp*g->cvac);
  const v8float one(1.);
  const v8float zero(0.);
  const v8float c(sp->g->cvac);
  const v8float r8V(g->r8V);
  const v8float two_fifteenths(2./15.);
  const v8float one_third(1./3.);

  const v8int stride_10(g->sx);                   // (i,j,k) -> (i+1, j,   k  )
  const v8int stride_21(g->sy - g->sx);           // (i,j,k) -> (i-1, j+1, k  )
  const v8int stride_43(g->sz - g->sy - g->sx);   // (i,j,k) -> (i-1, j-1, k+1)

  v8float dx, dy, dz, ux, uy, uz, w, vx, vy, vz, ke_mc;
  v8float hax, hay, haz, cbx, cby, cbz;
  v8float v00, v01, v02, v03, v04, v05, v06, v07;
  v8float t, w0, w1, w2, w3, w4, w5, w6, w7;
  v8int ii;

  int n0, nq;

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, nq );
  p += n0;
  nq >>= 3;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=8 )
  {
    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------
    load_8x8_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx,
                 &p[4].dx, &p[5].dx, &p[6].dx, &p[7].dx,
                 dx, dy, dz, ii, ux, uy, uz, w );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( float * ALIGNED(32) ) ( f + ii( 0) );
    vp01 = ( float * ALIGNED(32) ) ( f + ii( 1) );
    vp02 = ( float * ALIGNED(32) ) ( f + ii( 2) );
    vp03 = ( float * ALIGNED(32) ) ( f + ii( 3) );
    vp04 = ( float * ALIGNED(32) ) ( f + ii( 4) );
    vp05 = ( float * ALIGNED(32) ) ( f + ii( 5) );
    vp06 = ( float * ALIGNED(32) ) ( f + ii( 6) );
    vp07 = ( float * ALIGNED(32) ) ( f + ii( 7) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles and half advance with E
    //--------------------------------------------------------------------------
    load_8x8_tr( vp00, vp01, vp02, vp03,
		             vp04, vp05, vp06, vp07,
		             hax, v00, v01, v02, hay, v03, v04, v05 );

    hax = qdt_2mc*fma( fma( v02, dy, v01 ), dz, fma( v00, dy, hax ) );
    hay = qdt_2mc*fma( fma( v05, dz, v04 ), dx, fma( v03, dz, hay ) );


    load_8x8_tr( vp00+8, vp01+8, vp02+8, vp03+8,
		             vp04+8, vp05+8, vp06+8, vp07+8,
		             haz, v00, v01, v02, cbx, v03, cby, v04 );

    haz = qdt_2mc*fma( fma( v02, dx, v01 ), dy, fma( v00, dx, haz ) );

    //--------------------------------------------------------------------------
    // Load B interpolation data for particles.
    //--------------------------------------------------------------------------
    cbx = fma( v03, dx, cbx );
    cby = fma( v04, dy, cby );

    load_8x2_tr( vp00+16, vp01+16, vp02+16, vp03+16,
                 vp04+16, vp05+16, vp06+16, vp07+16,
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

#   define INCREMENT(offset)                                \
    transpose( v00, v01, v02, v03, v04, v05, v06, v07 );    \
    increment_8x1( vp00+offset, v00 );                      \
    increment_8x1( vp01+offset, v01 );                      \
    increment_8x1( vp02+offset, v02 );                      \
    increment_8x1( vp03+offset, v03 );                      \
    increment_8x1( vp04+offset, v04 );                      \
    increment_8x1( vp05+offset, v05 );                      \
    increment_8x1( vp06+offset, v06 );                      \
    increment_8x1( vp07+offset, v07 )

#   define ACCUM_HYDRO(w)                                   \
    vp00 = ( float * ALIGNED(32) ) ( h + ii( 0) );          \
    vp01 = ( float * ALIGNED(32) ) ( h + ii( 1) );          \
    vp02 = ( float * ALIGNED(32) ) ( h + ii( 2) );          \
    vp03 = ( float * ALIGNED(32) ) ( h + ii( 3) );          \
    vp04 = ( float * ALIGNED(32) ) ( h + ii( 4) );          \
    vp05 = ( float * ALIGNED(32) ) ( h + ii( 5) );          \
    vp06 = ( float * ALIGNED(32) ) ( h + ii( 6) );          \
    vp07 = ( float * ALIGNED(32) ) ( h + ii( 7) );          \
    t   = qsp*w;                                            \
    v00 = t*vx;       /* w vx */                            \
    v01 = t*vy;       /* w vy */                            \
    v02 = t*vz;       /* w vz */                            \
    v03 = t;          /* w */                              \
    if(charge_weight) {                                    \
    t   = mspc*w;                                           \
    dx  = t*ux;                                             \
    dy  = t*uy;                                             \
    dz  = t*uz;                                             \
    v04 = dx;         /* m c w ux */                        \
    v05 = dy;         /* m c w uy */                        \
    v06 = dz;         /* m c w uz */                        \
    v07 = t*ke_mc;    /* m c w (gamma-1) */                 \
    INCREMENT(0);                                           \
    v00 = dx*vx;      /* m c w ux vx */                     \
    v01 = dy*vy;      /* m c w uy vy */                     \
    v02 = dz*vz;      /* m c w uz vz */                     \
    v03 = dy*vz;      /* m c w uy vz */                     \
    v04 = dz*vx;      /* m c w uz vx */                     \
    v05 = dx*vy;      /* m c w ux vy */                     \
    v06 = zero;       /* pad[0] */                          \
    v07 = zero;       /* pad[1] */                          \
    INCREMENT(8);                                           \
    } else {INCREMENT(0); }

    // when charge_weight is false and we are doing mass weighting we could use
    // v04..v07 for properties of the next particles, but that would require a
    // lot more surgery

    /**/             ACCUM_HYDRO(w0); // Cell i,j,k
    ii += stride_10; ACCUM_HYDRO(w1); // Cell i+1,j,k
    ii += stride_21; ACCUM_HYDRO(w2); // Cell i,j+1,k
    ii += stride_10; ACCUM_HYDRO(w3); // Cell i+1,j+1,k
    ii += stride_43; ACCUM_HYDRO(w4); // Cell i,j,k+1
    ii += stride_10; ACCUM_HYDRO(w5); // Cell i+1,j,k+1
    ii += stride_21; ACCUM_HYDRO(w6); // Cell i,j+1,k+1
    ii += stride_10; ACCUM_HYDRO(w7); // Cell i+1,j+1,k+1

#   undef ACCUM_HYDRO
#   undef INCREMENT

  }

}

#else

void
hydro_p_pipeline_v8( hydro_p_pipeline_args_t * args,
                     int pipeline_rank,
                     int n_pipeline) {
  ERROR(("No hydro_p_pipeline_v8 implementation."));
}

#endif
