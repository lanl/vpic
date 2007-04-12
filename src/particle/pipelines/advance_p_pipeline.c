#include <particle_pipelines.h>

/* The host passes the negative of n_pipeline to this function for the
   host to process the final incomplete particle quad.  Note that this
   function _cannot_ use n_pipeline directly as it does not know which
   pipeline dispatcher called it! */

void
advance_p_pipeline( advance_p_pipeline_args_t * args,
                    int pipeline_rank ) {
  particle_t           * ALIGNED p0  = args->p;
  int                            n   = args->n;
  const float                    q_m = args->q_m;
  particle_mover_t     * ALIGNED pm  = args->pm;
  int                            nm  = args->nm;
  accumulator_t        * ALIGNED a0  = args->a;
  const interpolator_t * ALIGNED f0  = args->f;
  const grid_t *                 g   = args->g;

  const float qdt_2mc        = 0.5*q_m*g->dt/g->cvac;
  const float cdt_dx         = g->cvac*g->dt/g->dx;
  const float cdt_dy         = g->cvac*g->dt/g->dy;
  const float cdt_dz         = g->cvac*g->dt/g->dz;
  const float one            = 1.;
  const float one_third      = 1./3.;
  const float two_fifteenths = 2./15.;

  particle_t           * ALIGNED p;
  const interpolator_t * ALIGNED f;
  float                * ALIGNED a;

  int ii;
  float dx, dy, dz, ux, uy, uz, q;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4, v5;

  if( pipeline_rank<0 ) { /* Host does left over cleanup */
    pipeline_rank = -pipeline_rank; /* pipeline_rank == n_pipeline */

    /* Determine which particles the host processes, which movers are
       reserved for the host.  Note the host uses the first
       accumulator array. */

    p  = p0 + n; pm += nm;
    n &= 3;      nm  = n>nm ? nm : n;
    p -= n;      pm -= nm;

  } else { /* Pipelines do any rough equal number of particle quads */

    double n_target;

    /* Determine which particles to process in this pipeline */

    n_target = (double)(n>>2)/(double)n_pipeline;
    n  = (int)( n_target*(double) pipeline_rank    + 0.5 );
    p  = p0 + 4*n;
    n  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - n;
    n *= 4;

    /* Determine which movers are reserved for this pipeline */

    nm -= (args->n&3)>nm ? nm : (args->n&3); /* Reserve last movers for host */

    n_target = (double)nm / (double)n_pipeline; 
    nm  = (int)( n_target*(double) pipeline_rank    + 0.5 );
    pm += nm;
    nm  = (int)( n_target*(double)(pipeline_rank+1) + 0.5 ) - nm;

    /* Determine which accumulator array to use */

    a0 += (1+pipeline_rank)*(g->nx+2)*(g->ny+2)*(g->nz+2);

  }

  args->seg[pipeline_rank].pm = pm;
  args->seg[pipeline_rank].nm = nm;

  /* Process particles quads for this pipeline */

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
    q  = p->q;
    ux += hax;                               /* Half advance E */
    uy += hay;
    uz += haz;
    v0 = qdt_2mc/(float)sqrt(one + (ux*ux + (uy*uy + uz*uz))); /* Boris - scalars */
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
    ux += hax;                               /* Half advance E */
    uy += hay;
    uz += haz;
    p->ux = ux;                              /* Store momentum */
    p->uy = uy;
    p->uz = uz;
    v0 = one/(float)sqrt(one + (ux*ux+ (uy*uy + uz*uz))); /* Get norm displacement */
    ux *= cdt_dx;
    uy *= cdt_dy;
    uz *= cdt_dz;
    ux *= v0;
    uy *= v0;
    uz *= v0;
    v0 = dx + ux;                            /* Streak midpoint (inbnds) */
    v1 = dy + uy;
    v2 = dz + uz;
    v3 = v0 + ux;                            /* New position */
    v4 = v1 + uy;
    v5 = v2 + uz;
    if(  v3<=one && v4<=one && v5<=one &&          /* Check if inbnds */
        -v3<=one && -v4<=one && -v5<=one ) {

      /* Common case (inbnds)
         Note: accumulator values are 4 times the total physical charge that
         passed through the appropriate current quadrant in a time-step */

      p->dx = v3;                            /* Store new position */
      p->dy = v4;
      p->dz = v5;
      dx = v0;                               /* Streak mid */
      dy = v1;
      dz = v2;
      v5 = q*ux*uy*uz*one_third;             /* Compute correction */
      a = (float *)( a0 + ii );              /* Get accumulator */
#     define accumulate_j(X,Y,Z)                                  \
      v4  = q*u##X;   /* v2 = q ux                            */  \
      v1  = v4*d##Y;  /* v1 = q ux dy                         */  \
      v0  = v4-v1;    /* v0 = q ux (1-dy)                     */  \
      v1 += v4;       /* v1 = q ux (1+dy)                     */  \
      v4  = one+d##Z; /* v4 = 1+dz                            */  \
      v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */  \
      v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */  \
      v4  = one-d##Z; /* v4 = 1-dz                            */  \
      v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */  \
      v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */  \
      v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */  \
      v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */  \
      v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */  \
      v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */  \
      a[0] += v0;                                                 \
      a[1] += v1;                                                 \
      a[2] += v2;                                                 \
      a[3] += v3
      accumulate_j(x,y,z); a += 4;
      accumulate_j(y,z,x); a += 4;
      accumulate_j(z,x,y);
#     undef accumulate_j

    } else if( nm>0 ) {

      pm->dispx = ux;
      pm->dispy = uy;
      pm->dispz = uz;
      pm->i     = p - p0;
      if( move_p( p0, pm, a0, g ) ) pm++, nm--;

    }
  }

  args->seg[pipeline_rank].nm -= nm;
}
