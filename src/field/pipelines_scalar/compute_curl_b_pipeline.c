#include <field_pipelines.h>

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
             
void
compute_curl_b_pipeline( compute_curl_b_pipeline_args_t * args,
                         int pipeline_rank ) {
  field_t                      * ALIGNED f = args->f;
  const material_coefficient_t * ALIGNED m = args->m;
  const grid_t                 *         g = args->g;
  int n_voxel;

  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? g->cvac*g->dt/g->dx : 0;
  py = (ny>1) ? g->cvac*g->dt/g->dy : 0;
  pz = (nz>1) ? g->cvac*g->dt/g->dz : 0;

# if 0 /* Original non-pipelined version */
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

  /* Process the voxels assigned to this pipeline */
  
  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );
                               
  f0 = &f(x,  y,  z  );
  fx = &f(x-1,y,  z  );
  fy = &f(x,  y-1,z  );
  fz = &f(x,  y,  z-1);
  
  for( ; n_voxel; n_voxel-- ) {
    UPDATE_TCAX();
    UPDATE_TCAY();
    UPDATE_TCAZ(); 
    f0++; fx++;	fy++; fz++;
    
    x++;
    if( x>nx ) {
      x=2, y++;
      if( y>ny ) y=2, z++;
      f0 = &f(x,  y,  z  );
      fx = &f(x-1,y,  z  );
      fy = &f(x,  y-1,z  );
      fz = &f(x,  y,  z-1);
    }
  }
}
