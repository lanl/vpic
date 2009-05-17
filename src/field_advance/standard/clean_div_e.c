#define IN_sfa
#include "sfa_private.h"

#define f(x,y,z) f[ VOXEL(x,y,z,nx,ny,nz) ]

#define MARDER_EX() \
    f0->ex += m[f0->ematx].drivex*px*(fx->div_e_err-f0->div_e_err)

#define MARDER_EY() \
    f0->ey += m[f0->ematy].drivey*py*(fy->div_e_err-f0->div_e_err)

#define MARDER_EZ() \
    f0->ez += m[f0->ematz].drivez*pz*(fz->div_e_err-f0->div_e_err)

typedef struct pipeline_args {
  field_t                      * ALIGNED(128) f;
  const material_coefficient_t * ALIGNED(128) m;
  const grid_t                 *              g;
} pipeline_args_t;

static void
pipeline( pipeline_args_t * args,
          int pipeline_rank,
          int n_pipeline ) {
  field_t                      * ALIGNED(128) f = args->f;
  const material_coefficient_t * ALIGNED(128) m = args->m;
  const grid_t                 *              g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  float alphadt, px, py, pz;

  px = (nx>1) ? g->rdx : 0;
  py = (ny>1) ? g->rdy : 0;
  pz = (nz>1) ? g->rdz : 0;
  alphadt = 0.3888889/( px*px + py*py + pz*pz );
  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;
  
  // Process voxels assigned to this pipeline

  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x+1,y,  z  ); \
  fy = &f(x,  y+1,z  ); \
  fz = &f(x,  y,  z+1)

  LOAD_STENCIL();

  for( ; n_voxel; n_voxel-- ) {
    MARDER_EX();
    MARDER_EY();
    MARDER_EZ();
    f0++; fx++; fy++; fz++;

    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      LOAD_STENCIL();
    }
  }

# undef LOAD_STENCIL

}

void
clean_div_e( field_t                      * ALIGNED(128) f,
             const material_coefficient_t * ALIGNED(128) m,
             const grid_t                 *              g ) {
  pipeline_args_t args[1];
  
  float alphadt, px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field"));
  if( m==NULL ) ERROR(("Bad material coefficients"));
  if( g==NULL ) ERROR(("Bad grid"));

  // Do majority of field components in single pass on the pipelines.
  // The host handles stragglers.

# if 0 // Original non-pipelined version
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fx = &f(2,y,  z);
      fy = &f(1,y+1,z);
      fz = &f(1,y,z+1);
      for( x=1; x<=nx; x++ ) {
	MARDER_EX();
	MARDER_EY();
	MARDER_EZ();
	f0++; fx++; fy++; fz++;
      }
    }
  }
# endif

  args->f = f;
  args->m = m;
  args->g = g;

  EXEC_PIPELINES( pipeline, args, 0 );
  
  // Do left over field components on the host

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? g->rdx : 0;
  py = (ny>1) ? g->rdy : 0;
  pz = (nz>1) ? g->rdz : 0;
  alphadt = 0.3888889/( px*px + py*py + pz*pz );
  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;
  
  // Do left over ex
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,nz+1);
    fx = &f(2,y,nz+1);
    for( x=1; x<=nx; x++ ) {
      MARDER_EX();
      f0++; fx++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(2,ny+1,z);
    for( x=1; x<=nx; x++ ) {
      MARDER_EX();
      f0++; fx++;
    }
  }

  // Do left over ey
  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fy = &f(nx+1,y+1,z);
      MARDER_EY();
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,  nz+1);
    fy = &f(1,y+1,nz+1);
    for( x=1; x<=nx; x++ ) {
      MARDER_EY();
      f0++; fy++;
    }
  }

  // Do left over ez
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fz = &f(1,ny+1,z+1);
    for( x=1; x<=nx+1; x++ ) {
      MARDER_EZ();
      f0++; fz++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,z);
      fz = &f(nx+1,y,z+1);
      MARDER_EZ();
    }
  }

  WAIT_PIPELINES(); // FIXME: FINISH EVEN LATER??

  local_adjust_tang_e(f,g);
}
