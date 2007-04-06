#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
 
#define UPDATE_CBX() f0->cbx -= py*( fy->ez-f0->ez ) - pz*( fz->ey-f0->ey )
#define UPDATE_CBY() f0->cby -= pz*( fz->ex-f0->ex ) - px*( fx->ez-f0->ez )
#define UPDATE_CBZ() f0->cbz -= px*( fx->ey-f0->ey ) - py*( fy->ex-f0->ex )

void
advance_b_pipeline( advance_b_pipeline_args_t * args,
		    int pipeline_rank ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;
  float frac               = args->frac;
  int n_voxel;

  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? frac*g->cvac*g->dt/g->dx : 0;
  py = (ny>1) ? frac*g->cvac*g->dt/g->dy : 0;
  pz = (nz>1) ? frac*g->cvac*g->dt/g->dz : 0;

# if 0 /* Original non-pipeline version */
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(2,y,z);
      fy = &f(1,y+1,z);
      fz = &f(1,y,z+1);
      for( x=1; x<=nx; x++ ) {
	UPDATE_CBX();
	UPDATE_CBY();
	UPDATE_CBZ();
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
    UPDATE_CBX();
    UPDATE_CBY();
    UPDATE_CBZ();
    f0++; fx++; fy++; fz++;    

    x++;
    if( x>nx ) {
      x=1,  y++;
      if( y>ny ) y=1, z++;
      f0 = &f(x,  y,  z  );
      fx = &f(x+1,y,  z  );
      fy = &f(x,  y+1,z  );
      fz = &f(x,  y,  z+1);
    }
  }
}
