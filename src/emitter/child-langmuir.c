#include <math.h>    /* For sqrt and fabs */
#include <emitter.h>

/* FIXME: ERROR TRAPPING AND ERROR RETURNS */

void
child_langmuir( emitter_t * e,                     /* Actual emitter */
                const interpolator_t * ALIGNED fi, /* field interp */
                field_t * ALIGNED f,               /* rhob accum */
                accumulator_t * ALIGNED a,         /* injected current accum*/
                grid_t * g,                        /* matching grid */
                mt_handle rng ) {                           /* Random number generator */
  child_langmuir_t * args = (child_langmuir_t *)e->model_parameters;
  particle_injector_t pi;
  int i, n;

  /* Loop over all components of the region */

  for( n=0; n<e->n_component; n++ ) {

     i = EXTRACT_LOCAL_CELL( e->component[n] );

     /*  Notes:
        - ex, ey and ez in the interpolator for a cell are the average values of the
          fields in that cell
        - based on Child Law. i.e.
            J = 4 eps0 (2q/m)^(1/2) V^(3/2) / ( 9 L^2 )
          Total charge inject uses L as the width of a cell in the direction of emission
          normal. V is the avergae voltage drop across the cell (neglecting gauge
          issues).
        - Particles are emitted with a bi-Maxwellian distribution.
        - Particles are randomly distributed across the inject surface and have random
          ages. */

        /* For use in debugging: MESSAGE(("\n(x,y,z,i)=%f %f %f %i" \
                    "\n(ux,uy,uz,q)=%f %f %f %f"\
                    "\n(dx,dy,dz)= %f %f %f",\
                    pi.dx, pi.dy, pi.dz, pi.i,\
                    pi.ux, pi.uy, pi.uz, pi.q,\
                    pi.dispx, pi.dispy, pi.dispz));\*/

#    define EMIT_PARTICLES(X,Y,Z,dir) BEGIN_PRIMITIVE {                                \
       if( e->sp->q_m * (dir fi[i].e##X) > 0 ) {                                       \
         float qp, age;                                                                \
         int m = args->n_emit_per_face;                                                \
         /* This face can emit ... compute the charge of the emitted particles */      \
         qp = g->eps0*g->d##Y*g->d##Z*g->dt*                                           \
            sqrt((32./81.)*fabs(e->sp->q_m*fi[i].e##X*fi[i].e##X*fi[i].e##X)/g->d##X)/ \
            (float)m;                                                                  \
         if( e->sp->q_m<0 ) qp = -qp;                                                  \
         for( ; m; m-- ) {                                                             \
           pi.d##X     = -(dir 1);                                                     \
           pi.d##Y     = 2*mt_drand(rng)-1;                                            \
           pi.d##Z     = 2*mt_drand(rng)-1;                                            \
           pi.i        = i;                                                            \
           pi.u##X     = dir fabs( args->ut_para*mt_normal_drand(rng) );               \
           pi.u##Y     = args->ut_perp*mt_normal_drand(rng);                           \
           pi.u##Z     = args->ut_perp*mt_normal_drand(rng);                           \
           pi.q        = qp;                                                           \
           age         = mt_drand(rng);                                                \
           age        *= g->cvac*g->dt /                                               \
             sqrt( pi.u##X*pi.u##X + pi.u##Y*pi.u##Y + pi.u##Z*pi.u##Z + 1 );          \
           pi.disp##X  = pi.u##X*age/g->d##X;                                          \
           pi.disp##Y  = pi.u##Y*age/g->d##Y;                                          \
           pi.disp##Z  = pi.u##Z*age/g->d##Z;                                          \
           e->sp->nm  += inject_p( e->sp->p, e->sp->np, e->sp->pm+e->sp->nm,           \
                                   f, a, &pi, g );                                     \
           e->sp->np++;                                                                \
           /* FIXME: if inject_p fails, sp->np should not be incremented */            \
         }                                                                             \
       }                                                                               \
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

