/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <math.h>   /* For exp, sinh */
#include <material.h>

static float minf( float a, float b ) {
  return a<b ? a : b;
}

material_coefficient_t * ALIGNED
new_material_coefficients( const grid_t *g,
                           const material_t *m_list ) {
  float ax, ay, az, cg2;
  material_coefficient_t *mc, *material_coefficient;
  const material_t *mp;

  /* Check the input parameters */
  if( g==NULL ) {
    ERROR(("Invalid grid."));
    return NULL;
  }
  if( m_list==NULL ) {
    ERROR(("Empty material list."));
    return NULL;
  }

  /* Run sanity checks on the material list */
  ax = g->nx>1 ? g->cvac*g->dt/g->dx : 0; ax *= ax;
  ay = g->ny>1 ? g->cvac*g->dt/g->dy : 0; ay *= ay;
  az = g->nz>1 ? g->cvac*g->dt/g->dz : 0; az *= az;
  LIST_FOR_EACH(mp,m_list) {
    if( mp->sigmax/mp->epsx<0 )
      WARNING(("Material \"%s\" is an active medium along x", mp->name));
    if( mp->epsy*mp->muz<0 )
      WARNING(("Material \"%s\" has an imaginary x speed of light (ey)",
               mp->name));
    if( mp->epsz*mp->muy<0 )
      WARNING(("Material \"%s\" has an imaginary x speed of light (ez)",
               mp->name));
    if( mp->sigmay/mp->epsy<0 )
      WARNING(("Material \"%s\" is an active medium along y", mp->name));
    if( mp->epsz*mp->mux<0 )
      WARNING(("Material \"%s\" has an imaginary y speed of light (ez)",
               mp->name));
    if( mp->epsx*mp->muz<0 )
      WARNING(("Material \"%s\" has an imaginary y speed of light (ex)",
               mp->name));
    if( mp->sigmaz/mp->epsz<0 )
      WARNING(("Material \"%s\" is an an active medium along z", mp->name));
    if( mp->epsx*mp->muy<0 )
      WARNING(("Material \"%s\" has an imaginary z speed of light (ex)",
               mp->name));
    if( mp->epsy*mp->mux<0 )
      WARNING(("Material \"%s\" has an imaginary z speed of light (ey)",
               mp->name));
    cg2 = ax/minf(mp->epsy*mp->muz,mp->epsz*mp->muy) +
          ay/minf(mp->epsz*mp->mux,mp->epsx*mp->muz) +
          az/minf(mp->epsx*mp->muy,mp->epsy*mp->mux);
    if( cg2>=1 )
      WARNING(("Material \"%s\" Courant condition estimate = %e",
              mp->name, sqrt(cg2)));
  }

  /* Allocate the material coefficients */
  material_coefficient = (material_coefficient_t *ALIGNED)
    malloc_aligned( (m_list->id+1)*sizeof(material_coefficient_t),
                    preferred_alignment );
  if( material_coefficient==NULL ) {
    ERROR(("Could not allocate material coefficient array"));
    return NULL;
  }

  /* Fill up the material coefficient array */
  LIST_FOR_EACH( mp, m_list ) {
    mc = material_coefficient + mp->id;

    /* Advance E coefficients
       Note: mp->sigma{x,y,z} = 0 -> Non conductive
             mc->decay{x,y,z} = 0 -> Perfect conductor to numerical precision
             otherwise            -> Conductive */
    ax = ( mp->sigmax*g->dt ) / ( mp->epsx*g->eps0 );
    ay = ( mp->sigmay*g->dt ) / ( mp->epsy*g->eps0 );
    az = ( mp->sigmaz*g->dt ) / ( mp->epsz*g->eps0 );
    mc->decayx = exp(-ax);
    mc->decayy = exp(-ay);
    mc->decayz = exp(-az);
    if( ax==0 )              mc->drivex = 1./mp->epsx;
    else if( mc->decayx==0 ) mc->drivex = 0;
    else mc->drivex = 2.*exp(-0.5*ax)*sinh(0.5*ax) / (ax*mp->epsx);
    if( ay==0 )              mc->drivey = 1./mp->epsy;
    else if( mc->decayy==0 ) mc->drivey = 0;
    else mc->drivey = 2.*exp(-0.5*ay)*sinh(0.5*ay) / (ay*mp->epsy);
    if( az==0 )              mc->drivez = 1./mp->epsz;
    else if( mc->decayz==0 ) mc->drivez = 0;
    else mc->drivez = 2.*exp(-0.5*az)*sinh(0.5*az) / (az*mp->epsz);
    mc->rmux = 1./mp->mux;
    mc->rmuy = 1./mp->muy;
    mc->rmuz = 1./mp->muz;

    /* Clean div E coefficients
       Note: The charge density due to J = sigma E currents is not computed.
       Consequently, the divergence error inside conductors cannot computed.
       The divergence error multiplier is thus set to zero to ignore
       divergence errors inside conducting materials. */
    mc->nonconductive = ( ax==0 && ay==0 && az==0 ) ? 1. : 0.;
    mc->epsx = mp->epsx;
    mc->epsy = mp->epsy;
    mc->epsz = mp->epsz;
  }

  return material_coefficient;  
}

void delete_material_coefficients( material_coefficient_t ** ALIGNED mc ) {
  if( mc==NULL ) return;
  if( *mc!=NULL ) free_aligned(*mc);
  *mc = NULL;
}
