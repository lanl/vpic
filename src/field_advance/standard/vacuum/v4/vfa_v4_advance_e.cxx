#define IN_vfa_v4
#define HAS_V4_PIPELINE
#include "vfa_v4_private.h"

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

#define UPDATE_EX() f0->ex += ( ( py*(f0->cbz-fy->cbz) - pz*(f0->cby-fz->cby) ) - cj*f0->jfx )
#define UPDATE_EY() f0->ey += ( ( pz*(f0->cbx-fz->cbx) - px*(f0->cbz-fx->cbz) ) - cj*f0->jfy )
#define UPDATE_EZ() f0->ez += ( ( px*(f0->cby-fx->cby) - py*(f0->cbx-fy->cbx) ) - cj*f0->jfz )

// FIXME: MERGE WITH ABOVE DIRECTORY TO ELIMINATE THIS HIDEOUESNESS

typedef struct pipeline_args {
  field_t                      * ALIGNED(128) f;
  const grid_t                 *              g;
} pipeline_args_t;

extern "C" void
vfa_advance_e_pipeline( pipeline_args_t * args,
                        int pipeline_rank,
                        int n_pipeline );

static void
vfa_advance_e_pipeline_v4( pipeline_args_t * args,
                           int pipeline_rank,
                           int n_pipeline ) {

  using namespace v4;

  field_t                      * ALIGNED(128) f = args->f;
  const grid_t                 *              g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? g->cvac*g->dt*g->rdx : 0;
  const float py = (ny>1) ? g->cvac*g->dt*g->rdy : 0;
  const float pz = (nz>1) ? g->cvac*g->dt*g->rdz : 0;
  const float cj = g->dt/g->eps0;

  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);
  const v4float vcj(cj);

  v4float save0, dummy;

  v4float f0_ex,   f0_ey,   f0_ez;
  v4float f0_cbx,  f0_cby,  f0_cbz;
  v4float f0_tcax, f0_tcay, f0_tcaz;
  v4float f0_jfx,  f0_jfy,  f0_jfz;
  v4float          fx_cby,  fx_cbz;
  v4float fy_cbx,           fy_cbz;
  v4float fz_cbx,  fz_cby;

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

    load_4x4_tr( &f00->ex,   &f01->ex,   &f02->ex,   &f03->ex,   f0_ex,   f0_ey,   f0_ez,   save0 );
    load_4x3_tr( &f00->cbx,  &f01->cbx,  &f02->cbx,  &f03->cbx,  f0_cbx,  f0_cby,  f0_cbz         );
    load_4x3_tr( &f00->jfx,  &f01->jfx,  &f02->jfx,  &f03->jfx,  f0_jfx,  f0_jfy,  f0_jfz         );

    load_4x3_tr( &fx0->cbx,  &fx1->cbx,  &fx2->cbx,  &fx3->cbx,  dummy,   fx_cby,  fx_cbz         );
    load_4x3_tr( &fy0->cbx,  &fy1->cbx,  &fy2->cbx,  &fy3->cbx,  fy_cbx,  dummy,   fy_cbz         );
    load_4x2_tr( &fz0->cbx,  &fz1->cbx,  &fz2->cbx,  &fz3->cbx,  fz_cbx,  fz_cby   /**/           );

    f0_ex += ( (vpy*(f0_cbz-fy_cbz) - vpz*(f0_cby-fz_cby)) - vcj*f0_jfx );
    f0_ey += ( (vpz*(f0_cbx-fz_cbx) - vpx*(f0_cbz-fx_cbz)) - vcj*f0_jfy );
    f0_ez += ( (vpx*(f0_cby-fx_cby) - vpy*(f0_cbx-fy_cbx)) - vcj*f0_jfz );

    // Note: Unlike load_4x3 versus load_4x4, store_4x4 is much more efficient than store_4x3!

    store_4x4_tr( f0_ex,   f0_ey,   f0_ez,   save0, &f00->ex,    &f01->ex,    &f02->ex,    &f03->ex   );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL

}

void
vfa_v4_advance_e( field_t                      * ALIGNED(128) f,
                  const material_coefficient_t * ALIGNED(128) m, // Ignored
                  const grid_t                 *              g ) {
  pipeline_args_t args[1];
  
  float damp, px, py, pz, cj;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field"));
  if( g==NULL ) ERROR(("Bad grid"));

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  damp = g->damp;
  px = (nx>1) ? (1+damp)*g->cvac*g->dt*g->rdx : 0;
  py = (ny>1) ? (1+damp)*g->cvac*g->dt*g->rdy : 0;
  pz = (nz>1) ? (1+damp)*g->cvac*g->dt*g->rdz : 0;
  cj = g->dt/g->eps0;

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

  // Do majority interior in a single pass.  The host handles
  // stragglers.

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
	UPDATE_EX();
	UPDATE_EY();
	UPDATE_EZ();
	f0++; fx++; fy++; fz++;
      }
    }
  }
# endif
     
  args->f = f;
  args->g = g;

  EXEC_PIPELINES( vfa_advance_e, args, 0 );
  
  // Do left over interior ex
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fy = &f(1,y-1,z);
      fz = &f(1,y,  z-1);
      UPDATE_EX();
    }
  }

  // Do left over interior ey
  for( z=2; z<=nz; z++ ) {
    f0 = &f(2,1,z);
    fx = &f(1,1,z);
    fz = &f(2,1,z-1);
    for( x=2; x<=nx; x++ ) {
      UPDATE_EY();
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
      UPDATE_EZ();
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
      UPDATE_EX();
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
      UPDATE_EX();
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
      UPDATE_EX();
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
      UPDATE_EX();
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
      UPDATE_EY();
    }
  }
  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,z);
      fx = &f(nx,  y,z);
      fz = &f(nx+1,y,z-1);
      UPDATE_EY();
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(2,y,1);
    fx = &f(1,y,1);
    fz = &f(2,y,0);
    for( x=2; x<=nx; x++ ) {
      UPDATE_EY();
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
      UPDATE_EY();
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
      UPDATE_EZ();
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
      UPDATE_EZ();
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
      UPDATE_EZ();
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fx = &f(nx,  y,  z);
      fy = &f(nx+1,y-1,z);
      UPDATE_EZ();
    }
  }

  local_adjust_tang_e( f, g );
}
