#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define UPDATE_DERR_E() f0->div_e_err = m[f0->nmat].nonconductive * \
   ( px*( m[f0->ematx].epsx*f0->ex - m[fx->ematx].epsx*fx->ex ) +   \
     py*( m[f0->ematy].epsy*f0->ey - m[fy->ematy].epsy*fy->ey ) +   \
     pz*( m[f0->ematz].epsz*f0->ez - m[fz->ematz].epsz*fz->ez ) -   \
     cj*( f0->rhof + f0->rhob ) )

void
compute_div_e_err_pipeline( compute_div_e_err_pipeline_args_t * args,
                            int pipeline_rank,
                            int n_pipeline ) {
  field_t                      * ALIGNED f = args->f;
  const material_coefficient_t * ALIGNED m = args->m;
  const grid_t                 *         g = args->g;
  int n_voxel;

  float px, py, pz, cj;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;
  cj = 1./g->eps0;

# if 0 /* Original non-pipelined version */
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
	UPDATE_DERR_E();
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
    UPDATE_DERR_E();
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

