#define IN_sfa
#include "sfa_private.h"

field_advance_methods_t
_standard_field_advance[1] = { {

  // Field array construction / destruction

  new_field,
  delete_field,

  // Materal coefficient construction / destruction

  new_material_coefficients,
  delete_material_coefficients,

  // Time stepping interfaces

  advance_b,
  advance_e,

  // Diagnostic interfaces

  energy_f,

  // Accumulator interfaces

  clear_jf,   synchronize_jf,
  clear_rhof, synchronize_rho,

  // Initialize interface

  compute_rhob,
  compute_curl_b,

  // Shared face cleaning interface

  synchronize_tang_e_norm_b,
  
  // Electric field divergence cleaning interface

  compute_div_e_err,
  compute_rms_div_e_err,
  clean_div_e,

  // Magnetic field divergence cleaning interface

  compute_div_b_err,
  compute_rms_div_b_err,
  clean_div_b

} };

/*****************************************************************************/

field_t * ALIGNED(128)
new_field( grid_t * g ) {
  field_t * ALIGNED(128) f;

  if( g==NULL ) ERROR(("Bad grid."));
  if( g->nx<1 || g->ny<1 || g->nz<1 ) ERROR(("Bad resolution."));

  MALLOC_ALIGNED( f, (g->nx+2)*(g->ny+2)*(g->nz+2), 128 );
  CLEAR( f, (g->nx+2)*(g->ny+2)*(g->nz+2) );

  return f;
}

void
delete_field( field_t * ALIGNED(128) f ) {
  FREE_ALIGNED(f);
}

/*****************************************************************************/

static float minf( float a, float b ) {
  return a<b ? a : b;
}

material_coefficient_t * ALIGNED(128)
new_material_coefficients( grid_t * g,
                           material_t * m_list ) {
  float ax, ay, az, cg2;
  material_coefficient_t *mc, *material_coefficient;
  const material_t *m;
  int n_mat;

  // Check the input parameters
  if( g==NULL ) ERROR(("Invalid grid."));
  if( m_list==NULL ) ERROR(("Empty material list."));

  // Run sanity checks on the material list
  ax = g->nx>1 ? g->cvac*g->dt*g->rdx : 0; ax *= ax;
  ay = g->ny>1 ? g->cvac*g->dt*g->rdy : 0; ay *= ay;
  az = g->nz>1 ? g->cvac*g->dt*g->rdz : 0; az *= az;
  n_mat = 0;
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
    if( m->zetax!=0 || m->zetay!=0 || m->zetaz!=0 )
      WARNING(( "Standard field advance does not support magnetic conductivity yet." ));

    n_mat++;
  }

  // Allocate the material coefficients
  MALLOC_ALIGNED( material_coefficient, n_mat, 128 );

  // Fill up the material coefficient array
  LIST_FOR_EACH( m, m_list ) {
    mc = material_coefficient + m->id;

    // Advance E coefficients
    // Note: m ->sigma{x,y,z} = 0 -> Non conductive
    //       mc->decay{x,y,z} = 0 -> Perfect conductor to numerical precision
    //       otherwise            -> Conductive
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

    // Clean div E coefficients.  Note: The charge density due to J =
    // sigma E currents is not computed.  Consequently, the divergence
    // error inside conductors cannot computed.  The divergence error
    // multiplier is thus set to zero to ignore divergence errors
    // inside conducting materials.

    mc->nonconductive = ( ax==0 && ay==0 && az==0 ) ? 1. : 0.;
    mc->epsx = m->epsx;
    mc->epsy = m->epsy;
    mc->epsz = m->epsz;
  }

  return material_coefficient;  
}

void
delete_material_coefficients( material_coefficient_t * ALIGNED(128) mc ) {
  FREE_ALIGNED(mc);
}

/*****************************************************************************/

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
clear_jf( field_t      * ALIGNED(128) f,
          const grid_t *              g ) {
  int nx, ny, nz, x, y, z;
  field_t *f0;

  if( f==NULL ) ERROR(("Bad field")); 
  if( g==NULL ) ERROR(("Bad grid"));

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  for( z=0; z<=nz+1; z++ ) {
    for( y=0; y<=ny+1; y++ ) {
      f0 = &f(0,y,z);
      for( x=0; x<=nx+1; x++ ) {
	f0->jfx = 0;
	f0->jfy = 0;
	f0->jfz = 0;
	f0++;
      }
    }
  }
}

void
clear_rhof( field_t      * ALIGNED(128) f,
            const grid_t *              g ) {
  int nx, ny, nz, x, y, z;
  field_t *f0;

  if( f==NULL ) ERROR(("Bad field"));
  if( g==NULL ) ERROR(("Bad grid"));

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  for( z=0; z<=nz+1; z++ ) {
    for( y=0; y<=ny+1; y++ ) {
      f0 = &f(0,y,z);
      for( x=0; x<=nx+1; x++ ) {
	f0->rhof = 0;
	f0++;
      }
    }
  }
}

// FIXME: clear_jf_and_rhof CALL? OR MERGE ABOVE (MORE EFFICIENT TOO).
