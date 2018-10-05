#define IN_sf_interface

#include "sf_interface_private.h"

void
checkpt_interpolator_array( const interpolator_array_t * ia )
{
  CHECKPT( ia, 1 );
  CHECKPT_ALIGNED( ia->i, ia->g->nv, 128 );
  CHECKPT_PTR( ia->g );
}

interpolator_array_t *
restore_interpolator_array( void )
{
  interpolator_array_t * ia;
  RESTORE( ia );
  RESTORE_ALIGNED( ia->i );
  RESTORE_PTR( ia->g );
  return ia;
}

interpolator_array_t *
new_interpolator_array( grid_t * g )
{
  interpolator_array_t * ia;
  if( !g ) ERROR(( "NULL grid" ));
  MALLOC( ia, 1 );
  MALLOC_ALIGNED( ia->i, g->nv, 128 );
  CLEAR( ia->i, g->nv );
  ia->g = g;
  REGISTER_OBJECT( ia, checkpt_interpolator_array, restore_interpolator_array,
                   NULL );
  return ia;
}

void
delete_interpolator_array( interpolator_array_t * ia )
{
  if( !ia ) return;
  UNREGISTER_OBJECT( ia );
  FREE_ALIGNED( ia->i );
  FREE( ia );
}

//----------------------------------------------------------------------------//
// Include various programming model implementation files. For now, this is
// just the pipeline model. When there are more models, probably want these
// to be conditionally included.
//----------------------------------------------------------------------------//

#include "interpolator_array_pipeline.cc"

//----------------------------------------------------------------------------//
// Top level function to select and call the proper load_interpolator_array
// function.
//----------------------------------------------------------------------------//

void
load_interpolator_array( interpolator_array_t * RESTRICT ia,
                         const field_array_t * RESTRICT fa )
{
  if ( !ia              ||
       !fa              ||
       ia->g != fa->g )
  {
    ERROR( ( "Bad args" ) );
  }

  // Conditionally execute this when more abstractions are available.
  load_interpolator_array_pipeline( ia, fa );

# if 0 // Original non-pipelined version
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {

      pi = &fi(1,y,z);
      pf0 = &f(1,y,z);
      pfx = &f(2,y,z);
      pfy = &f(1,y+1,z);
      pfz = &f(1,y,z+1);
      pfyz = &f(1,y+1,z+1);
      pfzx = &f(2,y,z+1);
      pfxy = &f(2,y+1,z);

      for( x=1; x<=nx; x++ ) {

        // ex interpolation coefficients
        w0 = pf0->ex;
        w1 = pfy->ex;
        w2 = pfz->ex;
        w3 = pfyz->ex;
        pi->ex       = 0.25*(  w0 + w1 + w2 + w3 );
        pi->dexdy    = 0.25*( -w0 + w1 - w2 + w3 );
        pi->dexdz    = 0.25*( -w0 - w1 + w2 + w3 );
        pi->d2exdydz = 0.25*(  w0 - w1 - w2 + w3 );
        
        // ey interpolation coefficients
        w0 = pf0->ey;
        w1 = pfz->ey;
        w2 = pfx->ey;
        w3 = pfzx->ey;
        pi->ey       = 0.25*(  w0 + w1 + w2 + w3 );
        pi->deydz    = 0.25*( -w0 + w1 - w2 + w3 );
        pi->deydx    = 0.25*( -w0 - w1 + w2 + w3 );
        pi->d2eydzdx = 0.25*(  w0 - w1 - w2 + w3 );
        
        // ez interpolation coefficients
        w0 = pf0->ez;
        w1 = pfx->ez;
        w2 = pfy->ez;
        w3 = pfxy->ez;
        pi->ez       = 0.25*(  w0 + w1 + w2 + w3 );
        pi->dezdx    = 0.25*( -w0 + w1 - w2 + w3 );
        pi->dezdy    = 0.25*( -w0 - w1 + w2 + w3 );
        pi->d2ezdxdy = 0.25*(  w0 - w1 - w2 + w3 );
        
        // bx interpolation coefficients
        w0 = pf0->cbx;
        w1 = pfx->cbx;
        pi->cbx    = 0.5*(  w0 + w1 );
        pi->dcbxdx = 0.5*( -w0 + w1 );
        
        // by interpolation coefficients
        w0 = pf0->cby;
        w1 = pfy->cby;
        pi->cby    = 0.5*(  w0 + w1 );
        pi->dcbydy = 0.5*( -w0 + w1 );
        
        // bz interpolation coefficients
        w0 = pf0->cbz;
        w1 = pfz->cbz;
        pi->cbz    = 0.5*(  w0 + w1 );
        pi->dcbzdz = 0.5*( -w0 + w1 );

        pi++; pf0++; pfx++; pfy++; pfz++; pfyz++; pfzx++; pfxy++;
      }
    }
  }
# endif
}
