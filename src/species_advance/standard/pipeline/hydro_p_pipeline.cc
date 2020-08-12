#define IN_spa

#define HAS_V4_PIPELINE
#define HAS_V8_PIPELINE
#define HAS_V16_PIPELINE

#include "spa_private.h"

#include "../../../util/pipelines/pipelines_exec.h"

// accumulate_hydro_p adds the hydrodynamic fields associated with the
// supplied particle_list to the hydro array.  Trilinear interpolation
// is used.  hydro is known at the nodes at the same time as particle
// positions. No effort is made to fix up edges of the computational
// domain.  All particles on the list must be inbounds.  Note, the
// hydro jx,jy,jz are for diagnostic purposes only; they are not
// accumulated with a charge conserving algorithm.

void
accumulate_hydro_p_pipeline_scalar( accumulate_hydro_p_pipeline_args_t * args,
                                    int pipeline_rank,
                                    int n_pipeline )
{
  const species_t      *              sp = args->sp;
  const grid_t         *              g  = sp->g;
  /**/  hydro_t        * ALIGNED(128) h  = args->h + pipeline_rank * args->h_size;
  const particle_t     * ALIGNED(128) p  = sp->p;
  const interpolator_t * ALIGNED(128) f  = args->f;
  const bool               charge_weight = args->charge_weight;

  // Constants.

  const float qsp      = sp->q;
  const float qdt_2mc  = args->qdt_2mc;
  const float qdt_4mc2 = qdt_2mc / ( 2 * sp->g->cvac );
  const float mspc     = sp->g->cvac * args->msp;
  const float c        = sp->g->cvac;
  const float r8V      = sp->g->r8V;

  const float zero           = 0.0;
  const float one            = 1.0;
  const float one_third      = 1.0 / 3.0;
  const float two_fifteenths = 2.0 / 15.0;

  float dx, dy, dz, ux, uy, uz, w, vx, vy, vz, ke_mc;
  float t, w0, w1, w2, w3, w4, w5, w6, w7;
  int   i, n, n1, n0, np;

  int   stride_10;
  int   stride_21;
  int   stride_43;

  // Determine which particles this pipeline processes.

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, n0, n1 );

  n1 += n0;

  // c        = sp->g->cvac;
  // qsp      = sp->q;
  // mspc     = args->msp * c;
  // qdt_2mc  = args->qdt_2mc;
  // qdt_4mc2 = qdt_2mc / ( 2 * c );
  // r8V      = sp->g->r8V;

  stride_10 = VOXEL( 1, 0, 0, sp->g->nx, sp->g->ny, sp->g->nz ) -
              VOXEL( 0, 0, 0, sp->g->nx, sp->g->ny, sp->g->nz );

  stride_21 = VOXEL( 0, 1, 0, sp->g->nx, sp->g->ny, sp->g->nz ) -
              VOXEL( 1, 0, 0, sp->g->nx, sp->g->ny, sp->g->nz );

  stride_43 = VOXEL( 0, 0, 1, sp->g->nx, sp->g->ny, sp->g->nz ) -
              VOXEL( 1, 1, 0, sp->g->nx, sp->g->ny, sp->g->nz );

  for( n = n0; n < n1; n++ )
  {
    //--------------------------------------------------------------------------
    // Load particle data.
    //--------------------------------------------------------------------------

    dx = p[n].dx;
    dy = p[n].dy;
    dz = p[n].dz;
    i  = p[n].i;

    ux = p[n].ux;
    uy = p[n].uy;
    uz = p[n].uz;
    w  = p[n].w;

    //--------------------------------------------------------------------------
    // Load interpolation data for particles and half advance with E.
    //--------------------------------------------------------------------------

    ux += qdt_2mc * (      ( f[i].ex    + dy * f[i].dexdy    ) +
		      dz * ( f[i].dexdz + dy * f[i].d2exdydz ) );

    uy += qdt_2mc * (      ( f[i].ey    + dz * f[i].deydz    ) +
		      dx * ( f[i].deydx + dz * f[i].d2eydzdx ) );

    uz += qdt_2mc * (      ( f[i].ez    + dx * f[i].dezdx    ) +
		      dy * ( f[i].dezdy + dx * f[i].d2ezdxdy ) );

    //--------------------------------------------------------------------------
    // Load B interpolation data for particles.
    //--------------------------------------------------------------------------

    w5 = f[i].cbx + dx * f[i].dcbxdx;
    w6 = f[i].cby + dy * f[i].dcbydy;
    w7 = f[i].cbz + dz * f[i].dcbzdz;

    //--------------------------------------------------------------------------
    // Compute kinetic energy.
    //--------------------------------------------------------------------------

    // Boris rotation: curl scalars ( 0.5 in v0 for half rotate ) and
    // kinetic energy computation. Note: gamma - 1 = |u|^2 / ( gamma + 1 )
    // is the numerically accurate way to compute gamma - 1.
    ke_mc  = ux * ux + uy * uy + uz * uz; // ke_mc = |u|^2   (invariant)
    vz     = sqrt( one + ke_mc );         // vz    = gamma   (invariant)
    ke_mc *= c / ( vz + one );            // ke_mc = c|u|^2 / (gamma+1) = c * (gamma-1)
    vz     = c / vz;                      // vz    = c/gamma

    //--------------------------------------------------------------------------
    // Half Boris advance.
    //--------------------------------------------------------------------------

    w0  = qdt_4mc2 * vz;
    w1  = w5 * w5 + w6 * w6 + w7 * w7;  // |cB|^2
    w2  = w0 * w0 * w1;
    w3  = w0 * ( one + ( one_third ) * w2 * ( one + 0.4f * w2 ) );
    w4  = w3 / ( one + w1 * w3 * w3 );
    w4 += w4;

    // Boris rotation: uprime.
    w0 = ux + w3 * ( uy * w7 - uz * w6 );
    w1 = uy + w3 * ( uz * w5 - ux * w7 );
    w2 = uz + w3 * ( ux * w6 - uy * w5 );

    // Boris rotation: u.
    ux += w4 * ( w1 * w7 - w2 * w6 );
    uy += w4 * ( w2 * w5 - w0 * w7 );
    uz += w4 * ( w0 * w6 - w1 * w5 );

    //--------------------------------------------------------------------------
    // Compute velocity.
    //--------------------------------------------------------------------------

    vx = ux * vz;
    vy = uy * vz;
    vz = uz * vz;

    //--------------------------------------------------------------------------
    // Compute the trilinear coefficients.
    //--------------------------------------------------------------------------

    w0  = r8V * w;  // w0 = (1/8) (w/V)
    dx *= w0;       // dx = (1/8) (w/V) x
    w1  = w0 + dx;  // w1 = (1/8) (w/V) + (1/8) (w/V) x = (1/8) (w/V) (1+x)
    w0 -= dx;       // w0 = (1/8) (w/V) - (1/8) (w/V) x = (1/8) (w/V) (1-x)
    w3  = one + dy; // w3 = 1+y
    w2  = w0 * w3;  // w2 = (1/8) (w/V) (1-x) (1+y)
    w3 *= w1;       // w3 = (1/8) (w/V) (1+x) (1+y)
    dy  = one - dy; // dy = 1-y
    w0 *= dy;       // w0 = (1/8) (w/V) (1-x) (1-y)
    w1 *= dy;       // w1 = (1/8) (w/V) (1+x) (1-y)
    w7  = one + dz; // w7 = 1+z
    w4  = w0 * w7;  // w4 = (1/8) (w/V) (1-x) (1-y) (1+z) = (w/V) trilin_0 * Done
    w5  = w1 * w7;  // w5 = (1/8) (w/V) (1+x) (1-y) (1+z) = (w/V) trilin_1 * Done
    w6  = w2 * w7;  // w6 = (1/8) (w/V) (1-x) (1+y) (1+z) = (w/V) trilin_2 * Done
    w7 *= w3;       // w7 = (1/8) (w/V) (1+x) (1+y) (1+z) = (w/V) trilin_3 * Done
    dz  = one - dz; // dz = 1-z
    w0 *= dz;       // w0 = (1/8) (w/V) (1-x) (1-y) (1-z) = (w/V) trilin_4 * Done
    w1 *= dz;       // w1 = (1/8) (w/V) (1+x) (1-y) (1-z) = (w/V) trilin_5 * Done
    w2 *= dz;       // w2 = (1/8) (w/V) (1-x) (1+y) (1-z) = (w/V) trilin_6 * Done
    w3 *= dz;       // w3 = (1/8) (w/V) (1+x) (1+y) (1-z) = (w/V) trilin_7 * Done

    //--------------------------------------------------------------------------
    // Accumulate the hydro fields.
    //--------------------------------------------------------------------------

    #define ACCUM_HYDRO( wn )                                        \
    t         = qsp * wn;        /* t  = ( qsp   w / V ) trilin_n */ \
    h[i].jx  += t * vx;                                              \
    h[i].jy  += t * vy;                                              \
    h[i].jz  += t * vz;                                              \
    h[i].rho += t;                                                   \
    t         = mspc * wn;       /* t  = ( msp c w / V ) trilin_n */ \
    dx        = t * ux;          /* dx = ( px    w / V ) trilin_n */ \
    dy        = t * uy;                                              \
    dz        = t * uz;                                              \
    h[i].px  += dx;                                                  \
    h[i].py  += dy;                                                  \
    h[i].pz  += dz;                                                  \
    h[i].ke  += t * ke_mc;                                           \
    h[i].txx += dx * vx;                                             \
    h[i].tyy += dy * vy;                                             \
    h[i].tzz += dz * vz;                                             \
    h[i].tyz += dy * vz;                                             \
    h[i].tzx += dz * vx;                                             \
    h[i].txy += dx * vy

    /**/            ACCUM_HYDRO( w0 ); // Cell i,   j,   k
    i += stride_10; ACCUM_HYDRO( w1 ); // Cell i+1, j,   k
    i += stride_21; ACCUM_HYDRO( w2 ); // Cell i,   j+1, k
    i += stride_10; ACCUM_HYDRO( w3 ); // Cell i+1, j+1, k
    i += stride_43; ACCUM_HYDRO( w4 ); // Cell i,   j,   k+1
    i += stride_10; ACCUM_HYDRO( w5 ); // Cell i+1, j,   k+1
    i += stride_21; ACCUM_HYDRO( w6 ); // Cell i,   j+1, k+1
    i += stride_10; ACCUM_HYDRO( w7 ); // Cell i+1, j+1, k+1

    #undef ACCUM_HYDRO
  }
}

//----------------------------------------------------------------------------//
// Top level function to select and call the proper accumulate_hydro_p
// pipeline function.
//----------------------------------------------------------------------------//

void
accumulate_hydro_p_pipeline( hydro_array_t * RESTRICT ha,
                             const species_t * RESTRICT sp,
                             const interpolator_array_t * RESTRICT ia )
{
  DECLARE_ALIGNED_ARRAY( accumulate_hydro_p_pipeline_args_t, 128, args, 1 );

  if ( ! ha           ||
       ! sp           ||
       ! ia           ||
       ha->g != sp->g ||
       ha->g != ia->g )
  {
    ERROR( ( "Bad args." ) );
  }

  args->sp      = sp;
  args->f       = ia->i;
  args->h       = ha->h;
  args->h_size  = ha->stride;
  args->qdt_2mc = ( sp->q * sp->g->dt ) / ( 2 * sp->m * sp->g->cvac );
  args->msp     = sp->m;
  args->np      = sp->np;

  EXEC_PIPELINES( accumulate_hydro_p, args, 0 );

  WAIT_PIPELINES();
}
