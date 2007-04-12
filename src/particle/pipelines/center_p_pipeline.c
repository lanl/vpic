#include <particle_pipelines.h>

/* The host passes the negative of n_pipeline to this function for the
   host to process the final incomplete particle quad.  Note that this
   function _cannot_ use n_pipeline directly as it does not know which
   pipeline dispatcher called it! */

void
center_p_pipeline( center_p_pipeline_args_t * args,
                   int pipeline_rank,
                   int n_pipeline ) {
  particle_t           * ALIGNED p   = args->p;
  int                            n   = args->n;
  const float                    q_m = args->q_m;
  const interpolator_t * ALIGNED f0  = args->f;
  const grid_t *                 g   = args->g;

  const float qdt_2mc        = 0.5 *q_m*g->dt/g->cvac;
  const float qdt_4mc        = 0.25*q_m*g->dt/g->cvac; /* For half Boris rotate */
  const float one            = 1.;
  const float one_third      = 1./3.;
  const float two_fifteenths = 2./15.;

  const interpolator_t * ALIGNED f;

  int ii;
  float dx, dy, dz, ux, uy, uz;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4;

  if( pipeline_rank==n_pipeline ) { /* Host does left over cleanup */

    /* Determine which particles the host processes */

    p += n;
    n &= 3;
    p -= n;

  } else { /* Pipelines do any rough equal number of particle quads */

    /* Determine which particles to process in this pipeline */

    double n_target = (double)(n>>2)/(double)n_pipeline;
    n  = (int)( n_target*(double) pipeline_rank    + 0.5 );
    p += 4*n;
    n  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - n;
    n *= 4;

  }

  /* Process particles for this pipeline */

  for(;n;n--,p++) {
    dx = p->dx;                              /* Load position */
    dy = p->dy;
    dz = p->dz;
    ii = p->i;
    f = f0 + ii;                             /* Interpolate E */
    hax = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                    dz*( f->dexdz + dy*f->d2exdydz ) );
    hay = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                    dx*( f->deydx + dz*f->d2eydzdx ) );
    haz = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                    dy*( f->dezdy + dx*f->d2ezdxdy ) );
    cbx = f->cbx + dx*f->dcbxdx;             /* Interpolate B */
    cby = f->cby + dy*f->dcbydy;
    cbz = f->cbz + dz*f->dcbzdz;
    ux = p->ux;                              /* Load momentum */
    uy = p->uy;
    uz = p->uz;
    ux += hax;                               /* Half advance E */
    uy += hay;
    uz += haz;
    v0 = qdt_4mc/(float)sqrt(one + (ux*ux + (uy*uy + uz*uz))); /* Boris - scalars */
    v1 = cbx*cbx + (cby*cby + cbz*cbz);
    v2 = (v0*v0)*v1;
    v3 = v0*(one+v2*(one_third+v2*two_fifteenths));
    v4 = v3/(one+v1*(v3*v3));
    v4 += v4;
    v0 = ux + v3*( uy*cbz - uz*cby );        /* Boris - uprime */
    v1 = uy + v3*( uz*cbx - ux*cbz );
    v2 = uz + v3*( ux*cby - uy*cbx );
    ux += v4*( v1*cbz - v2*cby );            /* Boris - rotation */
    uy += v4*( v2*cbx - v0*cbz );
    uz += v4*( v0*cby - v1*cbx );
    p->ux = ux;                              /* Store momentum */
    p->uy = uy;
    p->uz = uz;
  }
}
