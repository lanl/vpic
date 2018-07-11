#define IN_spa

#define HAS_V4_PIPELINE
#define HAS_V8_PIPELINE
// #define HAS_V16_PIPELINE

#include "spa_private.h"

//----------------------------------------------------------------------------//
// Reference implementation for a center_p pipeline function which does not
// make use of explicit calls to vector intrinsic functions.
//----------------------------------------------------------------------------//

void
center_p_pipeline( center_p_pipeline_args_t * args,
                   int pipeline_rank,
                   int n_pipeline )
{
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  particle_t           * ALIGNED(32)  p;
  const interpolator_t * ALIGNED(16)  f;

  const float qdt_2mc        =     args->qdt_2mc;
  const float qdt_4mc        = 0.5*args->qdt_2mc; // For half Boris rotate
  const float one            = 1.;
  const float one_third      = 1./3.;
  const float two_fifteenths = 2./15.;

  float dx, dy, dz, ux, uy, uz;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4;
  int   ii;

  int first, n;

  // Determine which particles this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, first, n );
  p = args->p0 + first;

  // Process particles for this pipeline

  for(;n;n--,p++)
  {
    dx   = p->dx;                            // Load position
    dy   = p->dy;
    dz   = p->dz;
    ii   = p->i;
    f    = f0 + ii;                          // Interpolate E
    hax  = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                     dz*( f->dexdz + dy*f->d2exdydz ) );
    hay  = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                     dx*( f->deydx + dz*f->d2eydzdx ) );
    haz  = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                     dy*( f->dezdy + dx*f->d2ezdxdy ) );
    cbx  = f->cbx + dx*f->dcbxdx;            // Interpolate B
    cby  = f->cby + dy*f->dcbydy;
    cbz  = f->cbz + dz*f->dcbzdz;
    ux   = p->ux;                            // Load momentum
    uy   = p->uy;
    uz   = p->uz;
    ux  += hax;                              // Half advance E
    uy  += hay;
    uz  += haz;
    v0   = qdt_4mc/(float)sqrt(one + (ux*ux + (uy*uy + uz*uz)));
    /**/                                     // Boris - scalars
    v1   = cbx*cbx + (cby*cby + cbz*cbz);
    v2   = (v0*v0)*v1;
    v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
    v4   = v3/(one+v1*(v3*v3));
    v4  += v4;
    v0   = ux + v3*( uy*cbz - uz*cby );      // Boris - uprime
    v1   = uy + v3*( uz*cbx - ux*cbz );
    v2   = uz + v3*( ux*cby - uy*cbx );
    ux  += v4*( v1*cbz - v2*cby );           // Boris - rotation
    uy  += v4*( v2*cbx - v0*cbz );
    uz  += v4*( v0*cby - v1*cbx );
    p->ux = ux;                              // Store momentum
    p->uy = uy;
    p->uz = uz;
  }
}

//----------------------------------------------------------------------------//
// If using v4, include an implementation for center_p_pipeline_v4.
//----------------------------------------------------------------------------//

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#include "center_p_pipeline_v4.cc"

#endif

//----------------------------------------------------------------------------//
// If using v8, include an implementation for center_p_pipeline_v8.
//----------------------------------------------------------------------------//

#if defined(V8_ACCELERATION) && defined(HAS_V8_PIPELINE)

#include "center_p_pipeline_v8.cc"

#endif

//----------------------------------------------------------------------------//
// If using v16, include an implementation for center_p_pipeline_v16.
//----------------------------------------------------------------------------//

#if defined(V16_ACCELERATION) && defined(HAS_V16_PIPELINE)

// #include "center_p_pipeline_v16.cc"

#endif

//----------------------------------------------------------------------------//
// Top level function to select and call the proper center_p pipeline
// function.
//----------------------------------------------------------------------------//

void
center_p( /**/  species_t            * RESTRICT sp,
          const interpolator_array_t * RESTRICT ia )
{
  DECLARE_ALIGNED_ARRAY( center_p_pipeline_args_t, 128, args, 1 );

  if( !sp || !ia || sp->g!=ia->g ) ERROR(( "Bad args" ));

  // Have the pipelines do the bulk of particles in blocks and have the
  // host do the final incomplete block.

  args->p0      = sp->p;
  args->f0      = ia->i;
  args->qdt_2mc = (sp->q*sp->g->dt)/(2*sp->m*sp->g->cvac);
  args->np      = sp->np;

  EXEC_PIPELINES( center_p, args, 0 );
  WAIT_PIPELINES();
}
