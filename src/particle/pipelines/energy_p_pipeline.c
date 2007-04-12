#include <particle_pipelines.h>

/* The host passes the negative of n_pipeline to this function for the
   host to process the final incomplete particle quad.  Note that this
   function _cannot_ use n_pipeline directly as it does not know which
   pipeline dispatcher called it! */

void
energy_p_pipeline( energy_p_pipeline_args_t * args,
                   int pipeline_rank,
                   int n_pipeline ) {
  const particle_t     * ALIGNED p   = args->p;
  int                            n   = args->n;
  const float                    q_m = args->q_m;
  const interpolator_t * ALIGNED f0  = args->f;
  const grid_t *                 g   = args->g;

  const float qdt_2mc = 0.5*q_m*g->dt/g->cvac;
  const float one     = 1.;

  const interpolator_t * ALIGNED f;

  double en = 0;

  float dx, dy, dz;
  float v0, v1, v2;

  if( pipeline_rank==n_pipeline ) { /* Host does left over cleanup */

    /* Determine which particles the host processes, which movers are
       reserved for the host.  Note the host uses the first
       accumulator array. */

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

  /* Process particles quads for this pipeline */

  for(;n;n--,p++) {
    dx  = p->dx;
    dy  = p->dy;
    dz  = p->dz;
    f   = f0 + p->i;
    v0  = p->ux + qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                            dz*( f->dexdz + dy*f->d2exdydz ) );
    v1  = p->uy + qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                            dx*( f->deydx + dz*f->d2eydzdx ) );
    v2  = p->uz + qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                            dy*( f->dezdy + dx*f->d2ezdxdy ) );
    v0  = v0*v0 + v1*v1 + v2*v2;
    v0 /= (float)sqrt(one+v0)+one;
    en += (double)v0*(double)p->q;
  }

  args->en[pipeline_rank] = en;
}
