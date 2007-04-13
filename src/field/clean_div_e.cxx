#include <field.h>

// FIXME: This is at the very edge of being worthwhile to vectorize.
// For now, it is not.

#undef V4_ACCELERATION

#ifndef V4_ACCELERATION
#define CLEAN_DIV_E_PIPELINE (pipeline_func_t)clean_div_e_pipeline
#else
#define CLEAN_DIV_E_PIPELINE (pipeline_func_t)clean_div_e_pipeline_v4
#endif

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define MARDER_EX() \
    f0->ex += m[f0->ematx].drivex*px*(fx->div_e_err-f0->div_e_err)

#define MARDER_EY() \
    f0->ey += m[f0->ematy].drivey*py*(fy->div_e_err-f0->div_e_err)

#define MARDER_EZ() \
    f0->ez += m[f0->ematz].drivez*pz*(fz->div_e_err-f0->div_e_err)

typedef struct clean_div_e_pipeline_args {
  field_t                      * ALIGNED f;
  const material_coefficient_t * ALIGNED m;
  const grid_t                 *         g;
} clean_div_e_pipeline_args_t;

static void
clean_div_e_pipeline( clean_div_e_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline ) {
  field_t                      * ALIGNED f = args->f;
  const material_coefficient_t * ALIGNED m = args->m;
  const grid_t                 *         g = args->g;

  field_t * ALIGNED f0, * ALIGNED fx, * ALIGNED fy, * ALIGNED fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  float alphadt, px, py, pz;

  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;
  alphadt = 0.3888889/( px*px + py*py + pz*pz );
  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;
  
  /* Process voxels assigned to this pipeline */

  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
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

#ifdef V4_ACCELERATION
#error "V4 version not implemented"
#endif

void
clean_div_e( field_t * ALIGNED f,
             const material_coefficient_t * ALIGNED m,
             const grid_t * g ) {
  clean_div_e_pipeline_args_t args[1];
  
  float alphadt, px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field"));
  if( m==NULL ) ERROR(("Bad material coefficients"));
  if( g==NULL ) ERROR(("Bad grid"));

  /* Do majority of field components in single pass on the pipelines.
     The host handles stragglers. */

# if 0 /* Original non-pipelined version */
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

  PMETHOD.dispatch( CLEAN_DIV_E_PIPELINE, args, 0 );
  clean_div_e_pipeline( args, PMETHOD.n_pipeline, PMETHOD.n_pipeline );
  
  /* Do left over field components on the host */

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;
  alphadt = 0.3888889/( px*px + py*py + pz*pz );
  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;
  
  /* Do left over ex */
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

  /* Do left over ey */
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

  /* Do left over ez */
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

  PMETHOD.wait(); /* FIXME: FINISH EVEN LATER?? */

  local_adjust_tang_e(f,g);
}
