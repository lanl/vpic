#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define MARDER_CBX() f0->cbx += px*( f0->div_b_err - fx->div_b_err )
#define MARDER_CBY() f0->cby += py*( f0->div_b_err - fy->div_b_err )
#define MARDER_CBZ() f0->cbz += pz*( f0->div_b_err - fz->div_b_err )

void
clean_div_b_pipeline( clean_div_b_pipeline_args_t * args,
                      int pipeline_rank ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;
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
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
	MARDER_CBX();
	MARDER_CBY();
	MARDER_CBZ();
	f0++; fx++; fy++; fz++;
      }
    }
  }
# endif

  /* Process voxels assigned to this pipeline */
  
  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  f0 = &f(x,  y,  z  );
  fx = &f(x-1,y,  z  );
  fy = &f(x,  y-1,z  );
  fz = &f(x,  y,  z-1);
  
  for( ; n_voxel; n_voxel-- ) {
    MARDER_CBX();
    MARDER_CBY();
    MARDER_CBZ();
    f0++; fx++; fy++; fz++;
    
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
