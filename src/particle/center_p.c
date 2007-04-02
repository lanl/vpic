#include <particle_pipelines.h>
#include <math.h> /* For sqrt */

static void
center_p_host( particle_t           * ALIGNED p,    /* Particle array */
               int                            n,    /* Number of particles */
               const float                    q_m,  /* Charge to mass ratio */
               const interpolator_t * ALIGNED f0,   /* Interpolator array */
               const grid_t         *         g ) { /* Local domain grid parameters */
  int ii;
  float dx, dy, dz, ux, uy, uz;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4;
  const float qdt_2mc    = 0.5 *q_m*g->dt/g->cvac;
  const float qdt_4mc    = 0.25*q_m*g->dt/g->cvac; /* For half Boris rotate */
  const float one        = 1;
  const float one_third  = 1./3.;
  const float two_fifths = 0.4;
  const interpolator_t * ALIGNED f;

  p += n; n &= 3; p -= n;

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
    v0 = qdt_4mc/(float)sqrt(one + (ux*ux + uy*uy + uz*uz)); /* Boris - scalars */
    v1 = cbx*cbx + cby*cby + cbz*cbz;
    v2 = v0*v0*v1;
    v3 = v0*(one+one_third*v2*(one+two_fifths*v2));
    v4 = v3/(one + v1*v3*v3);
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

void
center_p( particle_t           * ALIGNED p,
          const int                      n,
          const float                    q_m,
          const interpolator_t * ALIGNED f,
          const grid_t         *         g ) {
  center_p_pipeline_args_t args[1];
  pipeline_request_t request[1];

  if( n<0     ) { ERROR(("Bad number of particles")); return; }
  if( f==NULL ) { ERROR(("Bad interpolator"));        return; }
  if( g==NULL ) { ERROR(("Bad grid"));                return; }

  /* Have the pipelines do the bulk of particles in quads and
     have the host do the final incomplete quad. */

  args->p   = p;
  args->n   = n;
  args->q_m = q_m;
  args->f   = f;
  args->g   = g;

  dispatch_pipelines( center_p_pipeline, args, 0, request );

  center_p_host( p, n, q_m, f, g );

  wait_for_pipelines( request );
}

