#include <particle_pipelines.h>

static double
energy_p_host( const particle_t     * ALIGNED p,    /* First particle to advance */
               int                            n,    /* Number particles to advance */
               const float                    q_m,  /* Charge to mass ratio */
               const interpolator_t * ALIGNED f0,   /* Interpolator array */
               const grid_t         *         g ) { /* Local domain grid parameters */
  float dx, dy, dz;
  float v0, v1, v2;
  const float qdt_2mc = 0.5*q_m*g->dt/g->cvac;
  const float one     = 1;
  const interpolator_t * ALIGNED f;
  double en = 0;

  p += n; n &= 3; p -= n;

  for(;n;n--,p++) {
    dx = p->dx; 
    dy = p->dy;
    dz = p->dz;
    f  = f0 + p->i;
    v0 = p->ux + qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                           dz*( f->dexdz + dy*f->d2exdydz ) );
    v1 = p->uy + qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                           dx*( f->deydx + dz*f->d2eydzdx ) );
    v2 = p->uz + qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                           dy*( f->dezdy + dx*f->d2ezdxdy ) );
    v0 = v0*v0 + v1*v1 + v2*v2;
    v0 /= (float)sqrt(one+v0)+one;
    en += (double)v0*(double)p->q;
  }

  return en;
}

double
energy_p( const particle_t     * ALIGNED p,
          const int                      n,
          const float                    q_m,
          const interpolator_t * ALIGNED f,
          const grid_t         *         g ) {
  energy_p_pipeline_args_t args[1];
  double local, global;
  int rank;

  if( n<0     ) ERROR(("Bad number of particles"));
  if( f==NULL ) ERROR(("Bad interpolator"));
  if( g==NULL ) ERROR(("Bad grid"));

  /* Have the pipelines do the bulk of particles in quads and
     have the host do the final incomplete quad. */

  args->p   = p;
  args->n   = n;
  args->q_m = q_m;
  args->f   = f;
  args->g   = g;

  dispatch_pipelines( energy_p_pipeline, args, 0 );

  local = energy_p_host( p, n, q_m, f, g );

  wait_for_pipelines();

  for( rank=0; rank<n_pipeline; rank++ ) local += args->en[rank];

  mp_allsum_d( &local, &global, 1, g->mp );

  return (double)g->cvac*(double)g->cvac*global/(double)q_m;
}

