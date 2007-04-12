#include <field_pipelines.h>

#define fi(x,y,z) fi[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
#define f(x,y,z)  f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
load_interpolator_pipeline( load_interpolator_pipeline_args_t * args,
			    int pipeline_rank,
                            int n_pipeline ) {
  interpolator_t * ALIGNED fi = args->fi;
  const field_t  * ALIGNED f  = args->f;
  const grid_t   *         g  = args->g;
  int n_voxel;
  
  float w0, w1, w2, w3;
  int x, y, z, nx, ny, nz;
  interpolator_t *pi;
  const field_t *pf0, *pfx, *pfy, *pfz, *pfyz, *pfzx, *pfxy;

  const float fourth = 0.25;
  const float half   = 0.5;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;

# if 0 /* Original non-pipelined version */  
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
# endif

  /* Process the voxels assigned to this pipeline */
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );
                               
  pi   = &fi(x,  y,  z  );
  pf0  =  &f(x,  y,  z  );
  pfx  =  &f(x+1,y,  z  );
  pfy  =  &f(x,  y+1,z  );
  pfz  =  &f(x,  y,  z+1);
  pfyz =  &f(x,  y+1,z+1);
  pfzx =  &f(x+1,y,  z+1);
  pfxy =  &f(x+1,y+1,z  );
  
  for( ; n_voxel; n_voxel-- ) {

    /* ex interpolation */
    w0 = pf0->ex;
    w1 = pfy->ex;
    w2 = pfz->ex;
    w3 = pfyz->ex;
    pi->ex       = fourth*( (w3 + w0) + (w1 + w2) );
    pi->dexdy    = fourth*( (w3 - w0) + (w1 - w2) );
    pi->dexdz    = fourth*( (w3 - w0) - (w1 - w2) );
    pi->d2exdydz = fourth*( (w3 + w0) - (w1 + w2) );

    /* ey interpolation coefficients */
    w0 = pf0->ey;
    w1 = pfz->ey;
    w2 = pfx->ey;
    w3 = pfzx->ey;
    pi->ey       = fourth*( (w3 + w0) + (w1 + w2) );
    pi->deydz    = fourth*( (w3 - w0) + (w1 - w2) );
    pi->deydx    = fourth*( (w3 - w0) - (w1 - w2) );
    pi->d2eydzdx = fourth*( (w3 + w0) - (w1 + w2) );

    /* ez interpolation coefficients */
    w0 = pf0->ez;
    w1 = pfx->ez;
    w2 = pfy->ez;
    w3 = pfxy->ez;
    pi->ez       = fourth*( (w3 + w0) + (w1 + w2) );
    pi->dezdx    = fourth*( (w3 - w0) + (w1 - w2) );
    pi->dezdy    = fourth*( (w3 - w0) - (w1 - w2) );
    pi->d2ezdxdy = fourth*( (w3 + w0) - (w1 + w2) );

    /* bx interpolation coefficients */
    w0 = pf0->cbx;
    w1 = pfx->cbx;
    pi->cbx    = half*( w1 + w0 );
    pi->dcbxdx = half*( w1 - w0 );

    /* by interpolation coefficients */
    w0 = pf0->cby;
    w1 = pfy->cby;
    pi->cby    = half*( w1 + w0 );
    pi->dcbydy = half*( w1 - w0 );

    /* bz interpolation coefficients */
    w0 = pf0->cbz;
    w1 = pfz->cbz;
    pi->cbz    = half*( w1 + w0 );
    pi->dcbzdz = half*( w1 - w0 );

    pi++; pf0++; pfx++; pfy++; pfz++; pfyz++; pfzx++; pfxy++;
    
    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      pi   = &fi(x,  y,  z  );
      pf0  =  &f(x,  y,  z  );
      pfx  =  &f(x+1,y,  z  );
      pfy  =  &f(x,  y+1,z  );
      pfz  =  &f(x,  y,  z+1);
      pfyz =  &f(x,  y+1,z+1);
      pfzx =  &f(x+1,y,  z+1);
      pfxy =  &f(x+1,y+1,z  );
    }
  }
}
