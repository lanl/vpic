#define IN_sf_interface
#define HAS_V4_PIPELINE
#include "sf_interface_private.h"

#define fi(x,y,z) fi[  INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
#define f(x,y,z)  f [  INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
#define nb(x,y,z) nb[6*INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define HAS_SPU_INTERPOLATOR \
  ( defined(CELL_PPU_BUILD) || defined(CELL_SPU_BUILD) ) && defined(USE_CELL_SPUS)

void
load_interpolator_pipeline( load_interpolator_pipeline_args_t * args,
			    int pipeline_rank,
                            int n_pipeline ) {
  interpolator_t * ALIGNED(128) fi = args->fi;
  const field_t  * ALIGNED(128) f  = args->f;
  const grid_t   *              g  = args->g;

# if HAS_SPU_INTERPOLATOR
  const int64_t * ALIGNED(128) nb = args->g->neighbor;
  const int64_t * ALIGNED(16)  pnb;
# endif

  interpolator_t * ALIGNED(16) pi;

  const field_t  * ALIGNED(16) pf0;
  const field_t  * ALIGNED(16) pfx,  * ALIGNED(16) pfy,  * ALIGNED(16) pfz;
  const field_t  * ALIGNED(16) pfyz, * ALIGNED(16) pfzx, * ALIGNED(16) pfxy;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float fourth = 0.25;
  const float half   = 0.5;

  float w0, w1, w2, w3;

  // Process the voxels assigned to this pipeline
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

# if HAS_SPU_INTERPOLATOR
# define LOAD_STENCIL()    \
  pi   = &fi(x,  y,  z  ); \
  pf0  =  &f(x,  y,  z  ); \
  pfx  =  &f(x+1,y,  z  ); \
  pfy  =  &f(x,  y+1,z  ); \
  pfz  =  &f(x,  y,  z+1); \
  pfyz =  &f(x,  y+1,z+1); \
  pfzx =  &f(x+1,y,  z+1); \
  pfxy =  &f(x+1,y+1,z  ); \
  pnb  = &nb(x,  y,  z  )
# else
# define LOAD_STENCIL()    \
  pi   = &fi(x,  y,  z  ); \
  pf0  =  &f(x,  y,  z  ); \
  pfx  =  &f(x+1,y,  z  ); \
  pfy  =  &f(x,  y+1,z  ); \
  pfz  =  &f(x,  y,  z+1); \
  pfyz =  &f(x,  y+1,z+1); \
  pfzx =  &f(x+1,y,  z+1); \
  pfxy =  &f(x+1,y+1,z  )
# endif

  LOAD_STENCIL();
  
  for( ; n_voxel; n_voxel-- ) {

    // ex interpolation
    w0 = pf0->ex;
    w1 = pfy->ex;
    w2 = pfz->ex;
    w3 = pfyz->ex;
    pi->ex       = fourth*( (w3 + w0) + (w1 + w2) );
    pi->dexdy    = fourth*( (w3 - w0) + (w1 - w2) );
    pi->dexdz    = fourth*( (w3 - w0) - (w1 - w2) );
    pi->d2exdydz = fourth*( (w3 + w0) - (w1 + w2) );

    // ey interpolation coefficients
    w0 = pf0->ey;
    w1 = pfz->ey;
    w2 = pfx->ey;
    w3 = pfzx->ey;
    pi->ey       = fourth*( (w3 + w0) + (w1 + w2) );
    pi->deydz    = fourth*( (w3 - w0) + (w1 - w2) );
    pi->deydx    = fourth*( (w3 - w0) - (w1 - w2) );
    pi->d2eydzdx = fourth*( (w3 + w0) - (w1 + w2) );

    // ez interpolation coefficients
    w0 = pf0->ez;
    w1 = pfx->ez;
    w2 = pfy->ez;
    w3 = pfxy->ez;
    pi->ez       = fourth*( (w3 + w0) + (w1 + w2) );
    pi->dezdx    = fourth*( (w3 - w0) + (w1 - w2) );
    pi->dezdy    = fourth*( (w3 - w0) - (w1 - w2) );
    pi->d2ezdxdy = fourth*( (w3 + w0) - (w1 + w2) );

    // bx interpolation coefficients
    w0 = pf0->cbx;
    w1 = pfx->cbx;
    pi->cbx    = half*( w1 + w0 );
    pi->dcbxdx = half*( w1 - w0 );

    // by interpolation coefficients
    w0 = pf0->cby;
    w1 = pfy->cby;
    pi->cby    = half*( w1 + w0 );
    pi->dcbydy = half*( w1 - w0 );

    // bz interpolation coefficients
    w0 = pf0->cbz;
    w1 = pfz->cbz;
    pi->cbz    = half*( w1 + w0 );
    pi->dcbzdz = half*( w1 - w0 );

#   if HAS_SPU_INTERPOLATOR
    pi->neighbor[0] = pnb[0];
    pi->neighbor[1] = pnb[1];
    pi->neighbor[2] = pnb[2];
    pi->neighbor[3] = pnb[3];
    pi->neighbor[4] = pnb[4];
    pi->neighbor[5] = pnb[5];
    pnb += 6;
#   endif

    pi++; pf0++; pfx++; pfy++; pfz++; pfyz++; pfzx++; pfxy++;

    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      LOAD_STENCIL();
    }
  }

# undef LOAD_STENCIL

}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

#error "SPU version not hooked up yet!"

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

using namespace v4;

void
load_interpolator_pipeline_v4( load_interpolator_pipeline_args_t * args,
                               int pipeline_rank,
                               int n_pipeline ) {
  interpolator_t * ALIGNED(128) fi = args->fi;
  const field_t  * ALIGNED(128) f  = args->f;
  const grid_t   *              g  = args->g;

# if HAS_SPU_INTERPOLATOR
  const int64_t * ALIGNED(128) nb = args->g->neighbor;
  const int64_t * ALIGNED(16) pnb;
# endif

  interpolator_t * ALIGNED(16) pi;

  const field_t * ALIGNED(16) pf0;
  const field_t * ALIGNED(16) pfx,  * ALIGNED(16) pfy,  * ALIGNED(16) pfz;
  const field_t * ALIGNED(16) pfyz, * ALIGNED(16) pfzx, * ALIGNED(16) pfxy;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const v4float fourth(0.25);
  const v4float half(  0.5 );

  const v4int   sgn_1_2(  0, 1<<31, 1<<31,     0 );
  const v4int   sgn_2_3(  0,     0, 1<<31, 1<<31 );
  const v4int   sgn_1_3(  0, 1<<31,     0, 1<<31 );
  const v4int   sel_0_1( -1,    -1,     0,     0 );

  v4float w0, w1, w2, w3;

  // Process the voxels assigned to this pipeline
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );
  
# if HAS_SPU_INTERPOLATOR
# define LOAD_STENCIL()    \
  pi   = &fi(x,  y,  z  ); \
  pf0  =  &f(x,  y,  z  ); \
  pfx  =  &f(x+1,y,  z  ); \
  pfy  =  &f(x,  y+1,z  ); \
  pfz  =  &f(x,  y,  z+1); \
  pfyz =  &f(x,  y+1,z+1); \
  pfzx =  &f(x+1,y,  z+1); \
  pfxy =  &f(x+1,y+1,z  ); \
  pnb  = &nb(x,  y,  z  )
# else
# define LOAD_STENCIL()    \
  pi   = &fi(x,  y,  z  ); \
  pf0  =  &f(x,  y,  z  ); \
  pfx  =  &f(x+1,y,  z  ); \
  pfy  =  &f(x,  y+1,z  ); \
  pfz  =  &f(x,  y,  z+1); \
  pfyz =  &f(x,  y+1,z+1); \
  pfzx =  &f(x+1,y,  z+1); \
  pfxy =  &f(x+1,y+1,z  )
# endif

  LOAD_STENCIL();
  
  for( ; n_voxel; n_voxel-- ) {

    // ex interpolation coefficients 
    w0 = toggle_bits( sgn_1_2, v4float( pf0->ex) ); // [ w0 -w0 -w0 w0 ]
    w1 =                       v4float( pfy->ex);   // [ w1  w1  w1 w1 ]
    w2 = toggle_bits( sgn_1_2, v4float( pfz->ex) ); // [ w2 -w2 -w2 w2 ]
    w3 =                       v4float(pfyz->ex);   // [ w3  w3  w3 w3 ]
    store_4x1( fourth*( ( w3 + w0 ) + toggle_bits( sgn_2_3, w1 + w2 ) ),
               &pi->ex );

    // ey interpolation coefficients 
    w0 = toggle_bits( sgn_1_2, v4float( pf0->ey) ); // [ w0 -w0 -w0 w0 ]
    w1 =                       v4float( pfz->ey);   // [ w1  w1  w1 w1 ]
    w2 = toggle_bits( sgn_1_2, v4float( pfx->ey) ); // [ w2 -w2 -w2 w2 ]
    w3 =                       v4float(pfzx->ey);   // [ w3  w3  w3 w3 ]
    store_4x1( fourth*( ( w3 + w0 ) + toggle_bits( sgn_2_3, w1 + w2 ) ),
               &pi->ey );

    // ez interpolation coefficients 
    w0 = toggle_bits( sgn_1_2, v4float( pf0->ez) ); // [ w0 -w0 -w0 w0 ]
    w1 =                       v4float( pfx->ez);   // [ w1  w1  w1 w1 ]
    w2 = toggle_bits( sgn_1_2, v4float( pfy->ez) ); // [ w2 -w2 -w2 w2 ]
    w3 =                       v4float(pfxy->ez);   // [ w3  w3  w3 w3 ]
    store_4x1( fourth*( ( w3 + w0 ) + toggle_bits( sgn_2_3, w1 + w2 ) ),
               &pi->ez );

    // bx and by interpolation coefficients 
    w0  = toggle_bits( sgn_1_3,
                       merge( sel_0_1,
                              v4float(pf0->cbx),
                              v4float(pf0->cby) ) ); // [ w0x -w0x w0y -w0y ]
    w1  =              merge( sel_0_1,
                              v4float(pfx->cbx),
                              v4float(pfy->cby) );   // [ w1x  w1x w1y  w1y ]
    store_4x1( half*( w1 + w0 ), &pi->cbx );

    // bz interpolation coefficients 
    w0  = toggle_bits( sgn_1_3, v4float(pf0->cbz) ); // [ w0 -w0 d/c d/c ]
    w1  =                       v4float(pfz->cbz);   // [ w1 -w1 d/c d/c ]
    store_4x1( half*( w1 + w0 ), &pi->cbz ); // Note: Padding after bz coeff!

#   if HAS_SPU_INTERPOLATOR
    copy_4x1( pi->neighbor,   pnb   );
    copy_4x1( pi->neighbor+2, pnb+2 );
    copy_4x1( pi->neighbor+4, pnb+4 );
    pnb += 6;
#   endif

    pi++; pf0++; pfx++; pfy++; pfz++; pfyz++; pfzx++; pfxy++;

    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      LOAD_STENCIL();
    }
  }

# undef LOAD_STENCIL

}

#endif

void
load_interpolator( interpolator_t * ALIGNED(128) fi,
                   const field_t  * ALIGNED(128) f,
                   const grid_t   *              g ) {
  load_interpolator_pipeline_args_t args[1];

  if( fi==NULL ) ERROR(("Bad interpolator"));
  if( f==NULL )  ERROR(("Bad field"));
  if( g==NULL )  ERROR(("Bad grid"));

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

  args->fi = fi;
  args->f  = f;
  args->g  = g;

  EXEC_PIPELINES( load_interpolator, args, 0 );
  WAIT_PIPELINES();
}
