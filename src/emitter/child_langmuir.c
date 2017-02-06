#define IN_emitter
#include "emitter_private.h"

/* Private interface *********************************************************/

// FIXME: IN CCUBE SHOULD J_norm BE PROP TO E^(3/2) or (E-THRESH)^(3/2)??

typedef struct child_langmuir {
  /**/  species_t            * sp;
  const interpolator_array_t * ia;
  /**/  field_array_t        * fa;
  /**/  accumulator_array_t  * aa;
  /**/  rng_t                * rng;
  int n_emit_per_face;
  float ut_para;
  float ut_perp;
  float thresh_e_norm;
  float norm;
} child_langmuir_t;

// Notes:
// - ex, ey and ez in the interpolator for a cell are the average
//   values of the fields in that cell
// - based on Child Law. i.e.
//     J = 4 eps0 (2q/m)^(1/2) V^(3/2) / ( 9 L^2 )
//   Total charge inject uses L as the width of a cell in the
//   direction of emission normal. V is the avergae voltage drop
//   across the cell (neglecting gauge issues).
// - Particles are emitted with a half-Maxwellian distribution.
//   See maxwellian_reflux for a derivation how this works.
// - Particles are randomly distributed across the inject surface
//   and have random ages.

void
emit_child_langmuir( child_langmuir_t * RESTRICT              cl,
                     const int        * RESTRICT ALIGNED(128) component,
                     int                                      n_component ) {
  /**/  species_t        * RESTRICT              sp  = cl->sp;
  const interpolator_t   * RESTRICT ALIGNED(128) fi  = cl->ia->i;
  /**/  field_t          * RESTRICT ALIGNED(128) f   = cl->fa->f;
  /**/  accumulator_t    * RESTRICT ALIGNED(128) a   = cl->aa->a;
  /**/  rng_t            * RESTRICT              rng = cl->rng;

  /**/  particle_t       * RESTRICT ALIGNED(128) p   = sp->p;
  /**/  particle_mover_t * RESTRICT ALIGNED(128) pm  = sp->pm;
  /**/  grid_t           * RESTRICT              g   = sp->g;

  const int max_np           = sp->max_np;
  const int max_nm           = sp->max_nm;
  const int np_emit_per_face = cl->n_emit_per_face;

  const float qsp     = sp->q;
  const float rdx     = g->rdx;
  const float rdy     = g->rdy;
  const float rdz     = g->rdz;
  const float cdt     = g->cvac*g->dt;
  const float norm    = ( cl->norm*g->eps0*g->dt ) /
                        ( sqrtf(fabsf(qsp*sp->m))*(float)np_emit_per_face );
  const float norm_x  = norm*sqrtf(rdx)*g->dy*g->dz;
  const float norm_y  = norm*sqrtf(rdy)*g->dz*g->dx;
  const float norm_z  = norm*sqrtf(rdz)*g->dx*g->dy;
  const float ut_para = cl->ut_para;
  const float ut_perp = cl->ut_perp;
  const float thresh  = fabsf(qsp)*cl->thresh_e_norm;

  int np = sp->np, np_skipped = 0;
  int nm = sp->nm, nm_skipped = 0;

  float w, ux, uy, uz;
  int c, cc, i, np_emit;

  // Loop over all components of the region

  for( c=0; c<n_component; c++ ) {
    cc = component[c];
    i  = EXTRACT_LOCAL_CELL( cc );

    // FIXME: COULD PROBABLY ACCELERATE BY GETTING RID OF SWITCH (USE
    // MAXWELLIAN_REFLUX TRICKS?)

#   define EMIT_PARTICLES(X,Y,Z,dir)                                    \
    w = fi[i].e##X;                                                     \
    if( dir qsp*w > thresh ) { /* This face can emit */                 \
      w = norm_##X*sqrtf(fabsf(w*w*w));                                 \
      for( np_emit=np_emit_per_face; np_emit; np_emit-- ) {             \
                                                                        \
        /* Emit the particle */                                         \
                                                                        \
        if( np>=max_np ) { np_skipped++; continue; }                    \
        u##X = dir ut_para*sqrtf(2*frande(rng));                        \
        u##Y = ut_perp*frandn(rng);                                     \
        u##Z = ut_perp*frandn(rng);                                     \
        p[np].d##X = -(dir 1);                                          \
        p[np].d##Y = 2*frand_c0(rng)-1;                                 \
        p[np].d##Z = 2*frand_c0(rng)-1;                                 \
        p[np].i    = i;                                                 \
        p[np].u##X = u##X;                                              \
        p[np].u##Y = u##Y;                                              \
        p[np].u##Z = u##Z;                                              \
        p[np].w    = w;                                                 \
        accumulate_rhob( f, p+np, g, -qsp );                            \
        np++;                                                           \
                                                                        \
        /* Age the particle */                                          \
                                                                        \
        if( nm>=max_nm ) { nm_skipped++; continue; }                    \
        w = ( frand_c0(rng)*cdt ) /                                     \
          sqrtf( ( u##X*u##X + u##Y*u##Y ) + ( u##Z*u##Z + 1 ) );       \
        pm[nm].disp##X = w*u##X*rd##X;                                  \
        pm[nm].disp##Y = w*u##Y*rd##Y;                                  \
        pm[nm].disp##Z = w*u##Z*rd##Z;                                  \
        pm[nm].i       = np-1;                                          \
        nm += move_p( p, pm, a, g, qsp );                               \
      }                                                                 \
    }

    switch( EXTRACT_COMPONENT_TYPE( cc ) ) {
    case BOUNDARY(-1, 0, 0): EMIT_PARTICLES(x,y,z,+) break;
    case BOUNDARY( 0,-1, 0): EMIT_PARTICLES(y,z,x,+) break;
    case BOUNDARY( 0, 0,-1): EMIT_PARTICLES(z,x,y,+) break;
    case BOUNDARY( 1, 0, 0): EMIT_PARTICLES(x,y,z,-) break;
    case BOUNDARY( 0, 1, 0): EMIT_PARTICLES(y,z,x,-) break;
    case BOUNDARY( 0, 0, 1): EMIT_PARTICLES(z,x,y,-) break;
    default: /* Not a cell face ... do not emit */   break;
    }

#   undef EMIT_PARTICLES

  }

  sp->np = np;
  sp->nm = nm;

  if( np_skipped ) WARNING(( "Insufficient local particle storage.  Did not emit %i "
                             "particles in emit_child_langmuir", np_skipped ));
  if( nm_skipped ) WARNING(( "Insufficient local particle mover storage.  Did not age %i "
                             "emitted particles in emit_child_langmuir", nm_skipped ));
}

void
checkpt_child_langmuir( const emitter_t * e ) {
  const child_langmuir_t * cl = (const child_langmuir_t *)e->params;
  CHECKPT( cl, 1 );
  CHECKPT_PTR( cl->sp );
  CHECKPT_PTR( cl->ia );
  CHECKPT_PTR( cl->fa );
  CHECKPT_PTR( cl->aa );
  CHECKPT_PTR( cl->rng );
  checkpt_emitter_internal( e );
}

emitter_t *
restore_child_langmuir( void ) {
  child_langmuir_t * cl;
  RESTORE( cl );
  RESTORE_PTR( cl->sp );
  RESTORE_PTR( cl->ia );
  RESTORE_PTR( cl->fa );
  RESTORE_PTR( cl->aa );
  RESTORE_PTR( cl->rng );
  return restore_emitter_internal( cl );
}

void
delete_child_langmuir( emitter_t * e ) {
  FREE( e->params );
  delete_emitter_internal( e );
}

/* Public interface **********************************************************/

emitter_t *
child_langmuir( /**/  species_t            * RESTRICT sp,
                const interpolator_array_t * RESTRICT ia,
                /**/  field_array_t        * RESTRICT fa,
                /**/  accumulator_array_t  * RESTRICT aa,
                /**/  rng_pool_t           * RESTRICT rp,
                int n_emit_per_face,
                float ut_para,
                float ut_perp,
                float thresh_e_norm,
                float norm ) {
  child_langmuir_t * cl;

  if( !sp || !ia || !fa || !aa || !rp ||
      sp->g!=ia->g || sp->g!=fa->g || sp->g!=aa->g ||
      n_emit_per_face<1 || ut_para<0  || ut_perp<0 || thresh_e_norm<0 )
    ERROR(( "Bad args" ));

  MALLOC( cl, 1 );
  cl->sp              = sp;
  cl->ia              = ia;
  cl->fa              = fa;
  cl->aa              = aa;
  cl->rng             = rp->rng[0];
  cl->n_emit_per_face = n_emit_per_face;
  cl->ut_para         = ut_para;
  cl->ut_perp         = ut_perp;
  cl->thresh_e_norm   = thresh_e_norm;
  cl->norm            = norm;
  return new_emitter_internal( cl,
                               (emit_func_t)emit_child_langmuir,
                               delete_child_langmuir,
                               (checkpt_func_t)checkpt_child_langmuir,
                               (restore_func_t)restore_child_langmuir,
                               NULL );
}

