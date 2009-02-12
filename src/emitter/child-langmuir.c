#include "emitter.h"

// FIXME: BOUNDARY CONDITIONS NEED SERIOUS REVAMP

void
child_langmuir( emitter_t            *              e,      // Actual emitter
                const interpolator_t * ALIGNED(128) fi,     // field interp
                field_t              * ALIGNED(16)  f,      // rhob accum
                accumulator_t        * ALIGNED(128) a,      // inject J accum
                grid_t               *              g,      // matching grid
                mt_rng_t             *              rng ) { // Rand num gen
  child_langmuir_t * args = (child_langmuir_t *)e->model_parameters;
  particle_t       * p;
  particle_mover_t * pm;
  int i, n;

  // Loop over all components of the region

  for( n=0; n<e->n_component; n++ ) {

     i = EXTRACT_LOCAL_CELL( e->component[n] );

     // Notes:
     // - ex, ey and ez in the interpolator for a cell are the average
     //   values of the fields in that cell
     // - based on Child Law. i.e.
     //     J = 4 eps0 (2q/m)^(1/2) V^(3/2) / ( 9 L^2 )
     //   Total charge inject uses L as the width of a cell in the
     //   direction of emission normal. V is the avergae voltage drop
     //   across the cell (neglecting gauge issues).
     // - Particles are emitted with a bi-Maxwellian distribution.
     // - Particles are randomly distributed across the inject surface
     //   and have random ages.

     /* For use in debugging:
        MESSAGE(("\n(x,y,z,i)=%f %f %f %i"
                 "\n(ux,uy,uz,q)=%f %f %f %f"
                 "\n(dx,dy,dz)= %f %f %f",
                 pi.dx, pi.dy, pi.dz, pi.i,
                 pi.ux, pi.uy, pi.uz, pi.q,
                 pi.dispx, pi.dispy, pi.dispz)); */

#    define EMIT_PARTICLES(X,Y,Z,dir) BEGIN_PRIMITIVE {                 \
       if( e->sp->q_m * (dir fi[i].e##X) > 0 ) {                        \
         float qp, age;                                                 \
         int m = args->n_emit_per_face;                                 \
         /* This face can emit ... compute the charge */                \
         /* of the emitted particles */                                 \
         qp = g->eps0*g->d##Y*g->d##Z*g->dt*                            \
            sqrt((32./81.)*fabs(e->sp->q_m*fi[i].e##X*fi[i].e##X*fi[i].e##X)/g->d##X)/ \
           (float)m;                                                    \
         if( e->sp->q_m<0 ) qp = -qp;                                   \
         if( e->sp->np+m>=e->sp->max_np ) {                             \
           WARNING(( "Insufficent room for emission; some particles skipped" )); \
           m = e->sp->max_np - e->sp->np;                               \
           if( m<0 ) m=0;                                               \
         }                                                              \
         for( ; m; m-- ) {                                              \
           p           = e->sp->p  + e->sp->np++;                       \
           p->d##X     = -(dir 1);                                      \
           p->d##Y     = 2*mt_drand_c(rng)-1;                           \
           p->d##Z     = 2*mt_drand_c(rng)-1;                           \
           p->i        = i;                                             \
           p->u##X     = dir fabs( args->ut_para*mt_drandn(rng) );      \
           p->u##Y     = args->ut_perp*mt_drandn(rng);                  \
           p->u##Z     = args->ut_perp*mt_drandn(rng);                  \
           p->q        = -qp; accumulate_rhob( f, p, g ); p->q = qp;    \
           if( e->sp->nm>=e->sp->max_nm ) {                             \
             WARNING(( "Insufficient movers to age emitted particle" )); \
             continue;                                                  \
           }                                                            \
           pm          = e->sp->pm + e->sp->nm;                         \
           age         = mt_drand_c0(rng);                              \
           age        *= g->cvac*g->dt /                                \
             sqrt( p->u##X*p->u##X +p->u##Y*p->u##Y + p->u##Z*p->u##Z + 1 ); \
           pm->disp##X = p->u##X*age/g->d##X;                           \
           pm->disp##Y = p->u##Y*age/g->d##Y;                           \
           pm->disp##Z = p->u##Z*age/g->d##Z;                           \
           pm->i       = e->sp->np-1;                                   \
           e->sp->nm  += move_p( e->sp->p, pm, a, g );                  \
         }                                                              \
       }                                                                \
     } END_PRIMITIVE

     switch( EXTRACT_COMPONENT_TYPE( e->component[n] ) ) {
     case BOUNDARY(-1, 0, 0): EMIT_PARTICLES(x,y,z,+); break;
     case BOUNDARY( 0,-1, 0): EMIT_PARTICLES(y,z,x,+); break;
     case BOUNDARY( 0, 0,-1): EMIT_PARTICLES(z,x,y,+); break;
     case BOUNDARY( 1, 0, 0): EMIT_PARTICLES(x,y,z,-); break;
     case BOUNDARY( 0, 1, 0): EMIT_PARTICLES(y,z,x,-); break;
     case BOUNDARY( 0, 0, 1): EMIT_PARTICLES(z,x,y,-); break;
     default: /* Not a cell face ... do not emit */    break;
     }

#   undef EMIT_PARTICLES

  }
}

