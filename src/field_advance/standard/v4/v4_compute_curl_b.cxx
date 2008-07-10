// FIXME: IT WOULD BE NICE TO REUSE THE SCALAR PIPELINE FROM STANDARD!
// Note: This is similar to advance_e_pipeline

#define IN_sfa_v4
#include "sfa_v4_private.h"

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define UPDATE_TCAX()						          \
  f0->tcax = py*(f0->cbz*m[f0->fmatz].rmuz-fy->cbz*m[fy->fmatz].rmuz) -   \
             pz*(f0->cby*m[f0->fmaty].rmuy-fz->cby*m[fz->fmaty].rmuy)

#define UPDATE_TCAY()						          \
  f0->tcay = pz*(f0->cbx*m[f0->fmatx].rmux-fz->cbx*m[fz->fmatx].rmux) -   \
             px*(f0->cbz*m[f0->fmatz].rmuz-fx->cbz*m[fx->fmatz].rmuz)

#define UPDATE_TCAZ()						          \
  f0->tcaz = px*(f0->cby*m[f0->fmaty].rmuy-fx->cby*m[fx->fmaty].rmuy) -   \
             py*(f0->cbx*m[f0->fmatx].rmux-fy->cbx*m[fy->fmatx].rmux)

typedef struct pipeline_args {
  field_t                      * ALIGNED(128) f;
  const material_coefficient_t * ALIGNED(128) m;
  const grid_t                 *              g;
} pipeline_args_t;

static void
pipeline( pipeline_args_t * args,
          int pipeline_rank,
          int n_pipeline ) {
  field_t                      * ALIGNED(16) f = args->f;
  const material_coefficient_t * ALIGNED(16) m = args->m;
  const grid_t                 *             g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? g->cvac*g->dt*g->rdx : 0;
  const float py = (ny>1) ? g->cvac*g->dt*g->rdy : 0;
  const float pz = (nz>1) ? g->cvac*g->dt*g->rdz : 0;

  // Process the voxels assigned to this pipeline
  
  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x-1,y,  z  ); \
  fy = &f(x,  y-1,z  ); \
  fz = &f(x,  y,  z-1)
  
  LOAD_STENCIL();

  for( ; n_voxel; n_voxel-- ) {
    UPDATE_TCAX();
    UPDATE_TCAY();
    UPDATE_TCAZ(); 
    f0++; fx++;	fy++; fz++;
    
    x++;
    if( x>nx ) {
      x=2, y++;
      if( y>ny ) y=2, z++;
      LOAD_STENCIL();
    }
  }

# undef LOAD_STENCIL

}

static void
v4_pipeline( pipeline_args_t * args,
             int pipeline_rank,
             int n_pipeline ) {

  using namespace v4;

  field_t                      * ALIGNED(16) f = args->f;
  const material_coefficient_t * ALIGNED(16) m = args->m;
  const grid_t                 *             g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? g->cvac*g->dt*g->rdx : 0;
  const float py = (ny>1) ? g->cvac*g->dt*g->rdy : 0;
  const float pz = (nz>1) ? g->cvac*g->dt*g->rdz : 0;

  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);

  v4float dummy;

  v4float f0_cbx,  f0_cby,  f0_cbz;
  v4float f0_tcax, f0_tcay, f0_tcaz;
  v4float          fx_cby,  fx_cbz;
  v4float fy_cbx,           fy_cbz;
  v4float fz_cbx,  fz_cby;
  v4float m_f0_rmux, m_f0_rmuy, m_f0_rmuz;
  v4float            m_fx_rmuy, m_fx_rmuz;
  v4float m_fy_rmux,            m_fy_rmuz;
  v4float m_fz_rmux, m_fz_rmuy;

  v4float f0_cbx_rmux, f0_cby_rmuy, f0_cbz_rmuz;

  field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel quad
  field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel quad +x neighbors
  field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel quad +y neighbors
  field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel quad +z neighbors

  // Process the voxels assigned to this pipeline 
  
  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  // Process the bulk of the voxels 4 at a time
                               
# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x-1,y,  z  ); \
  fy = &f(x,  y-1,z  ); \
  fz = &f(x,  y,  z-1)

# define NEXT_STENCIL(n) \
  f0##n = f0++;          \
  fx##n = fx++;          \
  fy##n = fy++;          \
  fz##n = fz++;          \
  x++;                   \
  if( x>nx ) {           \
    x=2, y++;            \
    if( y>ny ) y=2, z++; \
    LOAD_STENCIL();      \
  }

  LOAD_STENCIL();

  for( ; n_voxel>3; n_voxel-=4 ) {
    NEXT_STENCIL(0); NEXT_STENCIL(1); NEXT_STENCIL(2); NEXT_STENCIL(3);

    load_4x3_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx, f0_cbx, f0_cby, f0_cbz );
    load_4x3_tr( &fx0->cbx, &fx1->cbx, &fx2->cbx, &fx3->cbx, dummy,  fx_cby, fx_cbz );
    load_4x3_tr( &fy0->cbx, &fy1->cbx, &fy2->cbx, &fy3->cbx, fy_cbx, dummy,  fy_cbz );
    load_4x2_tr( &fz0->cbx, &fz1->cbx, &fz2->cbx, &fz3->cbx, fz_cbx, fz_cby  /**/   );

#   define LOAD_RMU( V, D ) load_4x1_tr( &m[f##V##0->fmat##D].rmu##D, \
                                         &m[f##V##1->fmat##D].rmu##D, \
                                         &m[f##V##2->fmat##D].rmu##D, \
                                         &m[f##V##3->fmat##D].rmu##D, \
                                         m_f##V##_rmu##D )

    LOAD_RMU(0,x); LOAD_RMU(0,y); LOAD_RMU(0,z);
    /**/           LOAD_RMU(x,y); LOAD_RMU(x,z);
    LOAD_RMU(y,x);                LOAD_RMU(y,z);
    LOAD_RMU(z,x); LOAD_RMU(z,y);

#   undef LOAD_RMU
    
    f0_cbx_rmux = f0_cbx * m_f0_rmux;
    f0_cby_rmuy = f0_cby * m_f0_rmuy;
    f0_cbz_rmuz = f0_cbz * m_f0_rmuz;

    f0_tcax = fms( vpy,fnms( fy_cbz,m_fy_rmuz, f0_cbz_rmuz ),
                   vpz*fnms( fz_cby,m_fz_rmuy, f0_cby_rmuy ) );

    f0_tcay = fms( vpz,fnms( fz_cbx,m_fz_rmux, f0_cbx_rmux ),
                   vpx*fnms( fx_cbz,m_fx_rmuz, f0_cbz_rmuz ) );

    f0_tcaz = fms( vpx,fnms( fx_cby,m_fx_rmuy, f0_cby_rmuy ),
                   vpy*fnms( fy_cbx,m_fy_rmux, f0_cbx_rmux ) );

    store_4x3_tr( f0_tcax, f0_tcay, f0_tcaz, &f00->tcax, &f01->tcax, &f02->tcax, &f03->tcax );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL

}

void
v4_compute_curl_b( field_t                      * ALIGNED(16) f,
                   const material_coefficient_t * ALIGNED(16) m,
                   const grid_t                 *             g ) {
  pipeline_args_t args[1];
  
  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field"));
  if( m==NULL ) ERROR(("Bad material coefficients"));
  if( g==NULL ) ERROR(("Bad grid"));

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? g->cvac*g->dt*g->rdx : 0;
  py = (ny>1) ? g->cvac*g->dt*g->rdy : 0;
  pz = (nz>1) ? g->cvac*g->dt*g->rdz : 0;

  /***************************************************************************
   * Begin tangential B ghost setup
   ***************************************************************************/
  
  begin_remote_ghost_tang_b( f, g );
  local_ghost_tang_b( f, g );

  /***************************************************************************
   * Update interior fields
   * Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
   * Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
   * Note: ez all (1:nx+1,1:ny+1,1:nz  ) interior (1:nx,1:ny,2:nz)
   ***************************************************************************/

  // Do bulk of the interior in a single pass in the pipelines. (The
  // host handles stragglers in the interior.)

  // FIXME: CHECK IF IT IS SAFE TO DISPATCH THE PIPELINES EVEN EARLIER
  // AND COMPLETE THEM EVEN LATER.  I DON'T THINK IT IS SAFE TO
  // DISPATCH THEM EARLIER UNDER ABSORBING BOUNDARY CONDITIONS BUT I
  // AM NOT SURE.  I THINK IT IS PROBABLY SAFE TO FINISH THEM LATER.

# if 0 // Original non-pipelined version
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
	UPDATE_TCAX();
	UPDATE_TCAY();
	UPDATE_TCAZ();
	f0++; fx++; fy++; fz++;
      }
    }
  }
# endif
     
  args->f = f;
  args->m = m;
  args->g = g;

  EXEC_PIPELINES( pipeline, args, 0 );
  
  // Do left over interior ex
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fy = &f(1,y-1,z);
      fz = &f(1,y,  z-1);
      UPDATE_TCAX();
    }
  }

  // Do left over interior ey
  for( z=2; z<=nz; z++ ) {
    f0 = &f(2,1,z);
    fx = &f(1,1,z);
    fz = &f(2,1,z-1);
    for( x=2; x<=nx; x++ ) {
      UPDATE_TCAY();
      f0++;
      fx++;
      fz++;
    }
  }

  // Do left over interior ez
  for( y=2; y<=ny; y++ ) {
    f0 = &f(2,y,  1);
    fx = &f(1,y,  1);
    fy = &f(2,y-1,1);
    for( x=2; x<=nx; x++ ) {
      UPDATE_TCAZ();
      f0++;
      fx++;
      fy++;
    }
  }

  WAIT_PIPELINES();
  
  /***************************************************************************
   * Finish tangential B ghost setup
   ***************************************************************************/

  end_remote_ghost_tang_b( f, g );

  /***************************************************************************
   * Update exterior fields
   ***************************************************************************/

  // Do exterior ex
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  1);
    fy = &f(1,y-1,1);
    fz = &f(1,y,  0);
    for( x=1; x<=nx; x++ ) {
      UPDATE_TCAX();
      f0++;
      fy++;
      fz++;
    }
  }
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  nz+1);
    fy = &f(1,y-1,nz+1);
    fz = &f(1,y,  nz);
    for( x=1; x<=nx; x++ ) {
      UPDATE_TCAX();
      f0++;
      fy++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fy = &f(1,0,z);
    fz = &f(1,1,z-1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_TCAX();
      f0++;
      fy++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fy = &f(1,ny,  z);
    fz = &f(1,ny+1,z-1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_TCAX();
      f0++;
      fy++;
      fz++;
    }
  }

  // Do exterior ey
  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(0,y,z);
      fz = &f(1,y,z-1);
      UPDATE_TCAY();
    }
  }
  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,z);
      fx = &f(nx,  y,z);
      fz = &f(nx+1,y,z-1);
      UPDATE_TCAY();
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(2,y,1);
    fx = &f(1,y,1);
    fz = &f(2,y,0);
    for( x=2; x<=nx; x++ ) {
      UPDATE_TCAY();
      f0++;
      fx++;
      fz++;
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(2,y,nz+1);
    fx = &f(1,y,nz+1);
    fz = &f(2,y,nz  );
    for( x=2; x<=nx; x++ ) {
      UPDATE_TCAY();
      f0++;
      fx++;
      fz++;
    }
  }

  // Do exterior ez
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fx = &f(0,1,z);
    fy = &f(1,0,z);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_TCAZ();
      f0++;
      fx++;
      fy++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(0,ny+1,z);
    fy = &f(1,ny,  z);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_TCAZ();
      f0++;
      fx++;
      fy++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fx = &f(0,y,  z);
      fy = &f(1,y-1,z);
      UPDATE_TCAZ();
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fx = &f(nx,  y,  z);
      fy = &f(nx+1,y-1,z);
      UPDATE_TCAZ();
    }
  }
}

