/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <material.h>

static float minf( float a, float b ) {
  return a<b ? a : b;
}

material_coefficient_t * ALIGNED
new_material_coefficients( const grid_t *g,
                           const material_t *m_list ) {
  float ax, ay, az, cg2;
  material_coefficient_t *mc, *material_coefficient;
  const material_t *m;

  /* Check the input parameters */
  if( g==NULL ) ERROR(("Invalid grid."));
  if( m_list==NULL ) ERROR(("Empty material list."));

  /* Run sanity checks on the material list */
  ax = g->nx>1 ? g->cvac*g->dt/g->dx : 0; ax *= ax;
  ay = g->ny>1 ? g->cvac*g->dt/g->dy : 0; ay *= ay;
  az = g->nz>1 ? g->cvac*g->dt/g->dz : 0; az *= az;
  LIST_FOR_EACH(m,m_list) {
    if( m->sigmax/m->epsx<0 )
      WARNING(("Material \"%s\" is an active medium along x", m->name));
    if( m->epsy*m->muz<0 )
      WARNING(("Material \"%s\" has an imaginary x speed of light (ey)",
               m->name));
    if( m->epsz*m->muy<0 )
      WARNING(("Material \"%s\" has an imaginary x speed of light (ez)",
               m->name));
    if( m->sigmay/m->epsy<0 )
      WARNING(("Material \"%s\" is an active medium along y", m->name));
    if( m->epsz*m->mux<0 )
      WARNING(("Material \"%s\" has an imaginary y speed of light (ez)",
               m->name));
    if( m->epsx*m->muz<0 )
      WARNING(("Material \"%s\" has an imaginary y speed of light (ex)",
               m->name));
    if( m->sigmaz/m->epsz<0 )
      WARNING(("Material \"%s\" is an an active medium along z", m->name));
    if( m->epsx*m->muy<0 )
      WARNING(("Material \"%s\" has an imaginary z speed of light (ex)",
               m->name));
    if( m->epsy*m->mux<0 )
      WARNING(("Material \"%s\" has an imaginary z speed of light (ey)",
               m->name));
    cg2 = ax/minf(m->epsy*m->muz,m->epsz*m->muy) +
          ay/minf(m->epsz*m->mux,m->epsx*m->muz) +
          az/minf(m->epsx*m->muy,m->epsy*m->mux);
    if( cg2>=1 )
      WARNING(( "Material \"%s\" Courant condition estimate = %e",
                m->name, sqrt(cg2) ));
  }

  /* Allocate the material coefficients */
  material_coefficient = (material_coefficient_t * ALIGNED)
    malloc_aligned( (m_list->id+1)*sizeof(material_coefficient_t),
                    preferred_alignment );
  if( material_coefficient==NULL )
    ERROR(("Could not allocate material coefficient array"));

  /* Fill up the material coefficient array */
  LIST_FOR_EACH( m, m_list ) {
    mc = material_coefficient + m->id;

    /* Advance E coefficients
       Note: m ->sigma{x,y,z} = 0 -> Non conductive
             mc->decay{x,y,z} = 0 -> Perfect conductor to numerical precision
             otherwise            -> Conductive */
    ax = ( m->sigmax*g->dt ) / ( m->epsx*g->eps0 );
    ay = ( m->sigmay*g->dt ) / ( m->epsy*g->eps0 );
    az = ( m->sigmaz*g->dt ) / ( m->epsz*g->eps0 );
    mc->decayx = exp(-ax);
    mc->decayy = exp(-ay);
    mc->decayz = exp(-az);
    if( ax==0 )              mc->drivex = 1./m->epsx;
    else if( mc->decayx==0 ) mc->drivex = 0;
    else mc->drivex = 2.*exp(-0.5*ax)*sinh(0.5*ax) / (ax*m->epsx);
    if( ay==0 )              mc->drivey = 1./m->epsy;
    else if( mc->decayy==0 ) mc->drivey = 0;
    else mc->drivey = 2.*exp(-0.5*ay)*sinh(0.5*ay) / (ay*m->epsy);
    if( az==0 )              mc->drivez = 1./m->epsz;
    else if( mc->decayz==0 ) mc->drivez = 0;
    else mc->drivez = 2.*exp(-0.5*az)*sinh(0.5*az) / (az*m->epsz);
    mc->rmux = 1./m->mux;
    mc->rmuy = 1./m->muy;
    mc->rmuz = 1./m->muz;

    /* Clean div E coefficients
       Note: The charge density due to J = sigma E currents is not computed.
       Consequently, the divergence error inside conductors cannot computed.
       The divergence error multiplier is thus set to zero to ignore
       divergence errors inside conducting materials. */
    mc->nonconductive = ( ax==0 && ay==0 && az==0 ) ? 1. : 0.;
    mc->epsx = m->epsx;
    mc->epsy = m->epsy;
    mc->epsz = m->epsz;
  }

  return material_coefficient;  
}

void delete_material_coefficients( material_coefficient_t ** ALIGNED mc ) {
  if( mc==NULL ) return;
  if( *mc!=NULL ) free_aligned(*mc);
  *mc = NULL;
}

