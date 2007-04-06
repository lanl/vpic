#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define MARDER_EX() \
    f0->ex += m[f0->ematx].drivex*px*(fx->div_e_err-f0->div_e_err)

#define MARDER_EY() \
    f0->ey += m[f0->ematy].drivey*py*(fy->div_e_err-f0->div_e_err);

#define MARDER_EZ() \
    f0->ez += m[f0->ematz].drivez*pz*(fz->div_e_err-f0->div_e_err);

void
clean_div_e_pipeline( clean_div_e_pipeline_args_t * args,
                      int pipeline_rank ) {
  field_t                      * ALIGNED f = args->f;
  const material_coefficient_t * ALIGNED m = args->m;
  const grid_t                 *         g = args->g;
  int n_voxel;

  float alphadt, px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

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

  /* Process voxels assigned to this pipeline */

  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  f0 = &f(x,  y,  z  );
  fx = &f(x+1,y,  z  );
  fy = &f(x,  y+1,z  );
  fz = &f(x,  y,  z+1);

  for( ; n_voxel; n_voxel-- ) {
    MARDER_EX();
    MARDER_EY();
    MARDER_EZ();
    f0++; fx++; fy++; fz++;

    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      f0 = &f(x,  y,  z  );
      fx = &f(x+1,y,  z  );
      fy = &f(x,  y+1,z  );
      fz = &f(x,  y,  z+1);
    }
  }
}

