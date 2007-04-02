#if 0 // Original non-pipelined non-vectorized version
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

        /* ex interpolation coefficients */
        w0 = pf0->ex;
        w1 = pfy->ex;
        w2 = pfz->ex;
        w3 = pfyz->ex;
        pi->ex       = 0.25*(  w0 + w1 + w2 + w3 );
        pi->dexdy    = 0.25*( -w0 + w1 - w2 + w3 );
        pi->dexdz    = 0.25*( -w0 - w1 + w2 + w3 );
        pi->d2exdydz = 0.25*(  w0 - w1 - w2 + w3 );
        
        /* ey interpolation coefficients */
        w0 = pf0->ey;
        w1 = pfz->ey;
        w2 = pfx->ey;
        w3 = pfzx->ey;
        pi->ey       = 0.25*(  w0 + w1 + w2 + w3 );
        pi->deydz    = 0.25*( -w0 + w1 - w2 + w3 );
        pi->deydx    = 0.25*( -w0 - w1 + w2 + w3 );
        pi->d2eydzdx = 0.25*(  w0 - w1 - w2 + w3 );
        
        /* ez interpolation coefficients */
        w0 = pf0->ez;
        w1 = pfx->ez;
        w2 = pfy->ez;
        w3 = pfxy->ez;
        pi->ez       = 0.25*(  w0 + w1 + w2 + w3 );
        pi->dezdx    = 0.25*( -w0 + w1 - w2 + w3 );
        pi->dezdy    = 0.25*( -w0 - w1 + w2 + w3 );
        pi->d2ezdxdy = 0.25*(  w0 - w1 - w2 + w3 );
        
        /* bx interpolation coefficients */
        w0 = pf0->cbx;
        w1 = pfx->cbx;
        pi->cbx    = 0.5*(  w0 + w1 );
        pi->dcbxdx = 0.5*( -w0 + w1 );
        
        /* by interpolation coefficients */
        w0 = pf0->cby;
        w1 = pfy->cby;
        pi->cby    = 0.5*(  w0 + w1 );
        pi->dcbydy = 0.5*( -w0 + w1 );
        
        /* bz interpolation coefficients */
        w0 = pf0->cbz;
        w1 = pfz->cbz;
        pi->cbz    = 0.5*(  w0 + w1 );
        pi->dcbzdz = 0.5*( -w0 + w1 );

        pi++; pf0++; pfx++; pfy++; pfz++; pfyz++; pfzx++; pfxy++;
      }
    }
  }
#endif

#include <field_pipelines.h>
#include CONCAT3(<,V4VERSION,>)

using namespace v4;

#define fi(x,y,z) fi[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
#define f(x,y,z)  f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
load_interpolator_pipeline( load_interpolator_pipeline_args_t * args,
			    int pipeline_rank ) {
  interpolator_t * ALIGNED fi = args->fi;
  const field_t  * ALIGNED f  = args->f;
  const grid_t   *         g  = args->g;

  int x, y, z, n_voxel;
  interpolator_t *pi;
  const field_t *pf0, *pfx, *pfy, *pfz, *pfyz, *pfzx, *pfxy;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const v4float fourth(0.25);
  const v4float half(  0.5 );

  const v4int   sgn_1_2(  0, 1<<31, 1<<31,     0 );
  const v4int   sgn_2_3(  0,     0, 1<<31, 1<<31 );
  const v4int   sgn_1_3(  0, 1<<31,     0, 1<<31 );
  const v4int   sel_0_1( -1,   -1,     0,     0 );

  v4float w0, w1, w2, w3;

  // Process the voxels assigned to this pipeline
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );
                               
# define LOAD_PTRS()       \
  pi   = &fi(x,  y,  z  ); \
  pf0  =  &f(x,  y,  z  ); \
  pfx  =  &f(x+1,y,  z  ); \
  pfy  =  &f(x,  y+1,z  ); \
  pfz  =  &f(x,  y,  z+1); \
  pfyz =  &f(x,  y+1,z+1); \
  pfzx =  &f(x+1,y,  z+1); \
  pfxy =  &f(x+1,y+1,z  )

  LOAD_PTRS();
  
  for( ; n_voxel; n_voxel-- ) {

    // ex interpolation coefficients 
    w0 = toggle_bits( sgn_1_2, v4float( pf0->ex) ); // [ w0 -w0 -w0 w0 ]
    w1 =                       v4float( pfy->ex);   // [ w1  w1  w1 w1 ]
    w2 = toggle_bits( sgn_1_2, v4float( pfz->ex) ); // [ w2 -w2 -w2 w2 ]
    w3 =                       v4float(pfyz->ex);   // [ w3  w3  w3 w3 ]
    store_4x1( fourth*( ( w3 + w0 ) + toggle_bits( sgn_2_3, w1 + w2 ) ), &pi->ex );

    // ey interpolation coefficients 
    w0 = toggle_bits( sgn_1_2, v4float( pf0->ey) ); // [ w0 -w0 -w0 w0 ]
    w1 =                       v4float( pfz->ey);   // [ w1  w1  w1 w1 ]
    w2 = toggle_bits( sgn_1_2, v4float( pfx->ey) ); // [ w2 -w2 -w2 w2 ]
    w3 =                       v4float(pfzx->ey);   // [ w3  w3  w3 w3 ]
    store_4x1( fourth*( ( w3 + w0 ) + toggle_bits( sgn_2_3, w1 + w2 ) ), &pi->ey );

    // ez interpolation coefficients 
    w0 = toggle_bits( sgn_1_2, v4float( pf0->ez) ); // [ w0 -w0 -w0 w0 ]
    w1 =                       v4float( pfx->ez);   // [ w1  w1  w1 w1 ]
    w2 = toggle_bits( sgn_1_2, v4float( pfy->ez) ); // [ w2 -w2 -w2 w2 ]
    w3 =                       v4float(pfxy->ez);   // [ w3  w3  w3 w3 ]
    store_4x1( fourth*( ( w3 + w0 ) + toggle_bits( sgn_2_3, w1 + w2 ) ), &pi->ez );

    // bx and by interpolation coefficients 
    w0  = toggle_bits( sgn_1_3, merge( sel_0_1, v4float(pf0->cbx), v4float(pf0->cby) ) ); // [ w0x -w0x w0y -w0y ]
    w1  =                       merge( sel_0_1, v4float(pfx->cbx), v4float(pfy->cby) );   // [ w1x  w1x w1y  w1y ]
    store_4x1( half*( w1 + w0 ), &pi->cbx );

    // bz interpolation coefficients 
    w0  = toggle_bits( sgn_1_3, v4float(pf0->cbz) ); // [ w0 -w0 d/c d/c ]
    w1  =                       v4float(pfz->cbz);   // [ w1 -w1 d/c d/c ]
    store_4x1( half*( w1 + w0 ), &pi->cbz );         // Note after bz coefficients is padding!

    pi++; pf0++; pfx++; pfy++; pfz++; pfyz++; pfzx++; pfxy++;
    
    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      LOAD_PTRS();
    }
  }
}
