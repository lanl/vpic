#define IN_spa

#include "spa_private.h"

#if defined(V4_ACCELERATION)

using namespace v4;

void
accumulate_hydro_p_pipeline_v4( accumulate_hydro_p_pipeline_args_t * args,
                     int pipeline_rank,
                     int n_pipeline)
{
  const species_t      *              sp = args->sp;
  const grid_t         *              g  = sp->g;
  /**/  hydro_t        * ALIGNED(128) h  = args->h + pipeline_rank*args->h_size;
  const particle_t     * ALIGNED(128) p  = sp->p;
  const interpolator_t * ALIGNED(128) f  = args->f;

  /**/  float          * ALIGNED(16)  vp00;
  /**/  float          * ALIGNED(16)  vp01;
  /**/  float          * ALIGNED(16)  vp02;
  /**/  float          * ALIGNED(16)  vp03;

  const v4float qsp(sp->q);
  const v4float qdt_2mc(args->qdt_2mc);
  const v4float qdt_4mc2(args->qdt_2mc / (2*g->cvac));
  const v4float mspc(args->msp*g->cvac);
  const v4float one(1.);
  const v4float zero(0.);
  const v4float c(sp->g->cvac);
  const v4float r8V(g->r8V);
  const v4float two_fifteenths(2./15.);
  const v4float one_third(1./3.);

  const v4int stride_10(g->sx);                   // (i,j,k) -> (i+1, j,   k  )
  const v4int stride_21(g->sy - g->sx);           // (i,j,k) -> (i-1, j+1, k  )
  const v4int stride_43(g->sz - g->sy - g->sx);   // (i,j,k) -> (i-1, j-1, k+1)

  v4float dx, dy, dz, ux, uy, uz, w, vx, vy, vz, ke_mc;
  v4float hax, hay, haz, cbx, cby, cbz;
  v4float v00, v01, v02, v03, v04, v05, v06;
  v4float t, w0, w1, w2, w3, w4, w5, w6, w7;
  v4int ii;

  int n0, nq;

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, nq );
  p += n0;
  nq >>= 2;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=4 )
  {
    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------
    load_4x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx, dx, dy, dz, ii );

    //--------------------------------------------------------------------------
    // Set field interpolation pointers.
    //--------------------------------------------------------------------------
    vp00 = ( float * ALIGNED(16) ) ( f + ii( 0) );
    vp01 = ( float * ALIGNED(16) ) ( f + ii( 1) );
    vp02 = ( float * ALIGNED(16) ) ( f + ii( 2) );
    vp03 = ( float * ALIGNED(16) ) ( f + ii( 3) );

    //--------------------------------------------------------------------------
    // Load interpolation data for particles and half advance with E
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00, vp01, vp02, vp03, hax, v00, v01, v02 );
    hax = qdt_2mc*fma( fma( v02, dy, v01 ), dz, fma( v00, dy, hax ) );


    load_4x4_tr( vp00+4, vp01+4, vp02+4, vp03+4, hay, v03, v04, v05 );
    hay = qdt_2mc*fma( fma( v05, dz, v04 ), dx, fma( v03, dz, hay ) );


    load_4x4_tr( vp00+8, vp01+8, vp02+8, vp03+8, haz, v00, v01, v02 );
    haz = qdt_2mc*fma( fma( v02, dx, v01 ), dy, fma( v00, dx, haz ) );

    //--------------------------------------------------------------------------
    // Load B interpolation data for particles.
    //--------------------------------------------------------------------------
    load_4x4_tr( vp00+12, vp01+12, vp02+12, vp03+12, cbx, v03, cby, v04 );
    cbx = fma( v03, dx, cbx );
    cby = fma( v04, dy, cby );

    load_4x2_tr( vp00+16, vp01+16, vp02+16, vp03+16, cbz, v05 );
    cbz = fma( v05, dz, cbz );

    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------
    load_4x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux, ux, uy, uz, w );

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
    w3  = one+dy;     // w3 = 1+y
    w2  = w0*w3;    // w2 = (1/8)(w/V)(1-x)(1+y)
    w3 *= w1;       // w3 = (1/8)(w/V)(1+x)(1+y)
    dy  = one-dy;     // dy = 1-y
    w0 *= dy;       // w0 = (1/8)(w/V)(1-x)(1-y)
    w1 *= dy;       // w1 = (1/8)(w/V)(1+x)(1-y)
    w7  = one+dz;     // w7 = 1+z
    w4  = w0*w7;    // w4 = (1/8)(w/V)(1-x)(1-y)(1+z) = (w/V) trilin_0 *Done
    w5  = w1*w7;    // w5 = (1/8)(w/V)(1+x)(1-y)(1+z) = (w/V) trilin_1 *Done
    w6  = w2*w7;    // w6 = (1/8)(w/V)(1-x)(1+y)(1+z) = (w/V) trilin_2 *Done
    w7 *= w3;       // w7 = (1/8)(w/V)(1+x)(1+y)(1+z) = (w/V) trilin_3 *Done
    dz  = one-dz;     // dz = 1-z
    w0 *= dz;       // w0 = (1/8)(w/V)(1-x)(1-y)(1-z) = (w/V) trilin_4 *Done
    w1 *= dz;       // w1 = (1/8)(w/V)(1+x)(1-y)(1-z) = (w/V) trilin_5 *Done
    w2 *= dz;       // w2 = (1/8)(w/V)(1-x)(1+y)(1-z) = (w/V) trilin_6 *Done
    w3 *= dz;       // w3 = (1/8)(w/V)(1+x)(1+y)(1-z) = (w/V) trilin_7 *Done

    //--------------------------------------------------------------------------
    // Accumulation.
    //--------------------------------------------------------------------------

#   define INCREMENT(offset)			\
    transpose( v00, v01, v02, v03 );		\
    increment_4x1( vp00+offset, v00 );		\
    increment_4x1( vp01+offset, v01 );		\
    increment_4x1( vp02+offset, v02 );		\
    increment_4x1( vp03+offset, v03 )

#   define ACCUM_HYDRO(w)                                   \
    vp00 = ( float * ALIGNED(16) ) ( h + ii( 0) );          \
    vp01 = ( float * ALIGNED(16) ) ( h + ii( 1) );          \
    vp02 = ( float * ALIGNED(16) ) ( h + ii( 2) );          \
    vp03 = ( float * ALIGNED(16) ) ( h + ii( 3) );          \
    t   = qsp*w;                                            \
    v00 = t*vx;       /* w vx */                            \
    v01 = t*vy;       /* w vy */                            \
    v02 = t*vz;       /* w vz */                            \
    v03 = t;          /* w */                               \
    INCREMENT(0);                                           \
    t   = mspc*w;                                           \
    dx  = t*ux;                                             \
    dy  = t*uy;                                             \
    dz  = t*uz;                                             \
    v00 = dx;         /* m c w ux */                        \
    v01 = dy;         /* m c w uy */                        \
    v02 = dz;         /* m c w uz */                        \
    v03 = t*ke_mc;    /* m c w (gamma-1) */                 \
    INCREMENT(4);                                           \
    v00 = dx*vx;      /* m c w ux vx */                     \
    v01 = dy*vy;      /* m c w uy vy */                     \
    v02 = dz*vz;      /* m c w uz vz */                     \
    v03 = dy*vz;      /* m c w uy vz */                     \
    INCREMENT(8);                                           \
    v00 = dz*vx;      /* m c w uz vx */                     \
    v01 = dx*vy;      /* m c w ux vy */                     \
    v02 = zero;       /* pad[0] */                          \
    v03 = zero;       /* pad[1] */                          \
    INCREMENT(12)

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
accumulate_hydro_p_pipeline_v4( accumulate_hydro_p_pipeline_args_t * args,
                     int pipeline_rank,
                     int n_pipeline) {
  ERROR(("No accumulate_hydro_p_pipeline_v4 implementation."));
}

#endif
