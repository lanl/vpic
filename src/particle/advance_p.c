/* FIXME: PARTICLE MOVERS NEED TO BE OVERALLOCATED IN STRUCTORS
   TO ACCOUNT FOR SPLITTING THE LIST BETWEEN HOST AND PARTICLE PIPELINES */

#include <particle_pipelines.h>
#include <unistd.h>

static void
advance_p_host( particle_t           * ALIGNED p,   /* Particle array */
                int                            n,   /* Number of particles */
                const float                    q_m, /* Charge to mass ratio */
                particle_mover_t     * ALIGNED pm,  /* Mover array */
                int                            nm,  /* Number of movers */
                accumulator_t        * ALIGNED a0,  /* Accumulator arrays */
                const interpolator_t * ALIGNED f0,  /* Interpolator array */
                const grid_t         *         g,   /* Local domain grid parameters */
                particle_mover_t    ** ALIGNED pm_seg,
                int                  *         nm_seg ) {  
  int ii;
  float dx, dy, dz, ux, uy, uz, q;
  float hax, hay, haz, cbx, cby, cbz;
  float v0, v1, v2, v3, v4, v5;
  const float qdt_2mc        = 0.5*q_m*g->dt/g->cvac;
  const float cdt_dx         = g->cvac*g->dt/g->dx;
  const float cdt_dy         = g->cvac*g->dt/g->dy;
  const float cdt_dz         = g->cvac*g->dt/g->dz;
  const float one            = 1.;
  const float one_third      = 1./3.;
  const float two_fifteenths = 2./15.;
  particle_t           * ALIGNED p0 = p;
  const interpolator_t * ALIGNED f;
  float                * ALIGNED a;

  /* Find the particles and movers reserved for the host */

  p  += n; pm += nm;
  n  &= 3; nm  = n>nm ? nm : n; *nm_seg = nm;
  p  -= n; pm -= nm;            *pm_seg = pm;

  /* Process particles assigned to the host */

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
  
  *nm_seg -= nm;
}

int
advance_p( particle_t           * ALIGNED p,
           const int                      n,
           const float                    q_m,
           particle_mover_t     * ALIGNED pm,
           int                            nm,       
           accumulator_t        * ALIGNED a,
           const interpolator_t * ALIGNED f,
           const grid_t         *         g ) {
  advance_p_pipeline_args_t args[1];
  int rank;

  if( p==NULL  ) ERROR(("Bad particle array"));
  if( n<0      ) ERROR(("Bad number of particles"));
  if( pm==NULL ) ERROR(("Bad particle mover"));
  if( nm<0     ) ERROR(("Bad number of movers"));
  if( a==NULL  ) ERROR(("Bad accumulator"));
  if( f==NULL  ) ERROR(("Bad interpolator"));
  if( g==NULL  ) ERROR(("Bad grid"));

  args->p   = p;
  args->n   = n;
  args->q_m = q_m;
  args->pm  = pm;
  args->nm  = nm;
  args->a   = a;
  args->f   = f;
  args->g   = g;

  dispatch_pipelines( advance_p_pipeline, args, 0 );

  /* Have the host processor do the incomplete quad if necessary.
     Note: This is overlapped with the pipelined processing.  As such,
     it uses an entire accumulator.  Reserving an entirely accumulator
     for the host processor to handle at most 3 particles is wasteful
     of memory.  It is anticipated that it may be useful at some point
     in the future have pipelines accumulating currents while the host
     processor is doing other more substantive work (e.g. accumulating
     currents from particles received from neighboring nodes).
     However, it is worth reconsidering this at some point in the
     future. */

  advance_p_host( p, n, q_m, pm, nm, a, f, g,
                  &args->seg[n_pipeline].pm,
                  &args->seg[n_pipeline].nm );

  wait_for_pipelines();

  /* FIXME: HIDEOUS HACK UNTIL BETTER PARTICLE MOVER SEMANTICS
     INSTALLED FOR DEALING WITH PIPELINES.  COMPACT THE PARTICLE
     MOVERS TO ELIMINATE HOLES IN THE ALLOCATION. */

  nm = 0;
  for( rank=0; rank<=n_pipeline; rank++ ) {
    if( pm+nm!=args->seg[rank].pm ) /* FIXME: C99 MEMCPY DOES NOT PERMIT OVERLAP */
      memcpy( pm+nm, args->seg[rank].pm, args->seg[rank].nm*sizeof(*pm) );
    nm += args->seg[rank].nm;
  }

  return nm;
}

