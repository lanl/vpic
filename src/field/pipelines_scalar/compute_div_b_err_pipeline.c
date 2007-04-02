#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void compute_div_b_err_pipeline( compute_div_b_err_pipeline_args_t * args,
                                 int pipeline_rank ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;
  int n_voxel;
  
  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;

# if 0 /* Original non-pipelined version */
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(2,y,z);
      fy = &f(1,y+1,z);
      fz = &f(1,y,z+1);
      for( x=1; x<=nx; x++ ) {
	f0->div_b_err = px*( fx->cbx - f0->cbx ) +
	                py*( fy->cby - f0->cby ) +
                        pz*( fz->cbz - f0->cbz );
	f0++; fx++; fy++; fz++;
      }
    }
  }
# endif

  /* Process the voxels assigned to this pipeline */
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  f0 = &f(x,  y,  z  );
  fx = &f(x+1,y,  z  );
  fy = &f(x,  y+1,z  );
  fz = &f(x,  y,  z+1);

  for( ; n_voxel; n_voxel-- ) {
    f0->div_b_err = px*( fx->cbx - f0->cbx ) +
	            py*( fy->cby - f0->cby ) +
                    pz*( fz->cbz - f0->cbz );
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
