#define IN_sfa
#include "sfa_private.h"

static field_advance_kernels_t sfa_kernels = {

  // Destructor

  delete_standard_field_array,

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

};

static float
minf( float a, 
      float b ) {
  return a<b ? a : b;
}

static sfa_params_t *
create_sfa_params( grid_t           * g,
                   const material_t * m_list,
                   float              damp ) {
  sfa_params_t * p;
  float ax, ay, az, cg2;
  material_coefficient_t *mc;
  const material_t *m;
  int n_mc;

  // Run sanity checks on the material list

  ax = g->nx>1 ? g->cvac*g->dt*g->rdx : 0; ax *= ax;
  ay = g->ny>1 ? g->cvac*g->dt*g->rdy : 0; ay *= ay;
  az = g->nz>1 ? g->cvac*g->dt*g->rdz : 0; az *= az;
  n_mc = 0;
  LIST_FOR_EACH(m,m_list) {
    if( m->sigmax/m->epsx<0 )
      WARNING(("\"%s\" is an active medium along x", m->name));
    if( m->epsy*m->muz<0 )
      WARNING(("\"%s\" has an imaginary x speed of light (ey)", m->name));
    if( m->epsz*m->muy<0 )
      WARNING(("\"%s\" has an imaginary x speed of light (ez)", m->name));
    if( m->sigmay/m->epsy<0 )
      WARNING(("\"%s\" is an active medium along y", m->name));
    if( m->epsz*m->mux<0 )
      WARNING(("\"%s\" has an imaginary y speed of light (ez)", m->name));
    if( m->epsx*m->muz<0 )
      WARNING(("\"%s\" has an imaginary y speed of light (ex)", m->name));
    if( m->sigmaz/m->epsz<0 )
      WARNING(("\"%s\" is an an active medium along z", m->name));
    if( m->epsx*m->muy<0 )
      WARNING(("\"%s\" has an imaginary z speed of light (ex)", m->name));
    if( m->epsy*m->mux<0 )
      WARNING(("\"%s\" has an imaginary z speed of light (ey)", m->name));
    cg2 = ax/minf(m->epsy*m->muz,m->epsz*m->muy) +
          ay/minf(m->epsz*m->mux,m->epsx*m->muz) +
          az/minf(m->epsx*m->muy,m->epsy*m->mux);
    if( cg2>=1 )
      WARNING(( "\"%s\" Courant condition estimate = %e", m->name, sqrt(cg2) ));
    if( m->zetax!=0 || m->zetay!=0 || m->zetaz!=0 )
      WARNING(( "\"%s\" magnetic conductivity is not supported" ));
    n_mc++;
  }

  // Allocate the sfa parameters

  MALLOC( p, 1 );
  MALLOC_ALIGNED( p->mc, n_mc+2, 128 );
  p->n_mc = n_mc;
  p->damp = damp;

  // Fill up the material coefficient array
  // FIXME: THIS IMPLICITLY ASSUMES MATERIALS ARE NUMBERED CONSECUTIVELY FROM
  // O.

  LIST_FOR_EACH( m, m_list ) {
    mc = p->mc + m->id;

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

  return p;
}

void
destroy_sfa_params( sfa_params_t * p ) {
  FREE_ALIGNED( p->mc );
  FREE( p );
}

/*****************************************************************************/

void
checkpt_standard_field_array( const field_array_t * fa ) {
  sfa_params_t * p = (sfa_params_t *)fa->params; 
  CHECKPT( fa, 1 );
  CHECKPT_ALIGNED( fa->f, fa->g->nv, 128 );
  CHECKPT_PTR( fa->g );
  CHECKPT( p, 1 );
  CHECKPT_ALIGNED( p->mc, p->n_mc, 128 );
  checkpt_field_advance_kernels( fa->kernel );
}

// FIXME: Use same new/delete/checkpt/restore structure as found in emitter
// and boundary(e.g. restore_field_advance_kernels =>
// return field_array_internal( params )).

field_array_t *
restore_standard_field_array( void ) {
  field_array_t * fa; 
  sfa_params_t * p;
  RESTORE( fa );
  RESTORE_ALIGNED( fa->f );
  RESTORE_PTR( fa->g );
  RESTORE( p );
  RESTORE_ALIGNED( p->mc );
  fa->params = p;
  restore_field_advance_kernels( fa->kernel );
  return fa;
}

field_array_t *
new_standard_field_array( grid_t           * RESTRICT g,
                          const material_t * RESTRICT m_list,
                          float                       damp ) {
  field_array_t * fa;
  if( !g || !m_list || damp<0 ) ERROR(( "Bad args" ));
  MALLOC( fa, 1 );
  MALLOC_ALIGNED( fa->f, g->nv, 128 );
  CLEAR( fa->f, g->nv );
  fa->g = g;
  fa->params = create_sfa_params( g, m_list, damp );
  fa->kernel[0] = sfa_kernels;
  if( !m_list->next ) {
    /* If there is only one material, then this material permeates all
       space and we can use high performance versions of some kernels. */
    fa->kernel->advance_e         = vacuum_advance_e;
    fa->kernel->energy_f          = vacuum_energy_f;
    fa->kernel->compute_rhob      = vacuum_compute_rhob;
    fa->kernel->compute_curl_b    = vacuum_compute_curl_b;
    fa->kernel->compute_div_e_err = vacuum_compute_div_e_err;
    fa->kernel->clean_div_e       = vacuum_clean_div_e;
  }

  REGISTER_OBJECT( fa, checkpt_standard_field_array,
                       restore_standard_field_array, NULL );
  return fa;
}

void
delete_standard_field_array( field_array_t * fa ) {
  if( !fa ) return;
  UNREGISTER_OBJECT( fa );
  destroy_sfa_params( (sfa_params_t *)fa->params );
  FREE_ALIGNED( fa->f );
  FREE( fa );
}

/*****************************************************************************/

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

void
clear_jf( field_array_t * RESTRICT fa ) {
  if( !fa ) ERROR(( "Bad args" ));
  field_t * RESTRICT ALIGNED(128) f = fa->f;
  const int nv = fa->g->nv;
  for( int v=0; v<nv; v++ ) f[v].jfx = 0, f[v].jfy = 0, f[v].jfz = 0;
}

void
clear_rhof( field_array_t * RESTRICT fa ) {
  if( !fa ) ERROR(( "Bad args" )); 
  field_t * RESTRICT ALIGNED(128) f = fa->f;
  const int nv = fa->g->nv;
  for( int v=0; v<nv; v++ ) f[v].rhof = 0;
}

// FIXME: ADD clear_jf_and_rhof CALL AND/OR ELIMINATE SOME OF THE ABOVE
// (MORE EFFICIENT TOO).
