#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define UPDATE_EX()						            \
  f0->tcax = ( py*(f0->cbz*m[f0->fmatz].rmuz-fy->cbz*m[fy->fmatz].rmuz) -   \
               pz*(f0->cby*m[f0->fmaty].rmuy-fz->cby*m[fz->fmaty].rmuy) ) - \
    damp*f0->tcax;                                                          \
  f0->ex = m[f0->ematx].decayx*f0->ex +                                     \
           m[f0->ematx].drivex*( f0->tcax - cj*f0->jfx )
#define UPDATE_EY()						            \
  f0->tcay = ( pz*(f0->cbx*m[f0->fmatx].rmux-fz->cbx*m[fz->fmatx].rmux) -   \
               px*(f0->cbz*m[f0->fmatz].rmuz-fx->cbz*m[fx->fmatz].rmuz) ) - \
    damp*f0->tcay;                                                          \
  f0->ey = m[f0->ematy].decayy*f0->ey +                                     \
           m[f0->ematy].drivey*( f0->tcay - cj*f0->jfy )
#define UPDATE_EZ()						            \
  f0->tcaz = ( px*(f0->cby*m[f0->fmaty].rmuy-fx->cby*m[fx->fmaty].rmuy) -   \
               py*(f0->cbx*m[f0->fmatx].rmux-fy->cbx*m[fy->fmatx].rmux) ) - \
    damp*f0->tcaz;                                                          \
  f0->ez = m[f0->ematz].decayz*f0->ez +                                     \
           m[f0->ematz].drivez*( f0->tcaz - cj*f0->jfz )

void
advance_e_pipeline( advance_e_pipeline_args_t * args,
                    int pipeline_rank,
                    int n_pipeline ) {
  field_t                      * ALIGNED f = args->f;
  const material_coefficient_t * ALIGNED m = args->m;
  const grid_t                 *         g = args->g;
  int n_voxel;

  float damp, px, py, pz, cj;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  damp = g->damp;
  px = (nx>1) ? (1+damp)*g->cvac*g->dt/g->dx : 0;
  py = (ny>1) ? (1+damp)*g->cvac*g->dt/g->dy : 0;
  pz = (nz>1) ? (1+damp)*g->cvac*g->dt/g->dz : 0;
  cj = g->dt/g->eps0;

# if 0 /* Original non-pipelined version */
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

  /* Process the voxels assigned to this pipeline */
  
  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );
                               
  f0 = &f(x,  y,  z  );
  fx = &f(x-1,y,  z  );
  fy = &f(x,  y-1,z  );
  fz = &f(x,  y,  z-1);
  
  for( ; n_voxel; n_voxel-- ) {
    UPDATE_EX();
    UPDATE_EY();
    UPDATE_EZ(); 
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
