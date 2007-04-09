#include <v4.h>
#ifdef V4_ACCELERATION
using namespace v4;

#if 0 // Original non-pipelined non-vectorized version 
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
#endif

#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
compute_div_b_err_pipeline_v4( compute_div_b_err_pipeline_args_t * args,
                               int pipeline_rank ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;

  int x, y, z, n_voxel;
  field_t *f0, *fx, *fy, *fz;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? 1./g->dx : 0;
  const float py = (ny>1) ? 1./g->dy : 0;
  const float pz = (nz>1) ? 1./g->dz : 0;

  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);

  v4float f0_cbx, f0_cby, f0_cbz; // Voxel quad magnetic fields
  v4float f0_div_b_err;           // Voxel quad div b errs
  v4float fx_cbx;                 // Voxel quad +x neighbor x magnetic fields
  v4float fy_cby;                 // Voxel quad +y neighbor y magnetic fields
  v4float fz_cbz;                 // Voxel quad +z neighbor z magnetic fields

  field_t * ALIGNED f00, * ALIGNED f01, * ALIGNED f02, * ALIGNED f03; // Voxel quad
  field_t * ALIGNED fx0, * ALIGNED fx1, * ALIGNED fx2, * ALIGNED fx3; // Voxel quad +x neighbors
  field_t * ALIGNED fy0, * ALIGNED fy1, * ALIGNED fy2, * ALIGNED fy3; // Voxel quad +x neighbors
  field_t * ALIGNED fz0, * ALIGNED fz1, * ALIGNED fz2, * ALIGNED fz3; // Voxel quad +x neighbors

  // Process the voxels assigned to this pipeline 
  
  n_voxel = distribute_voxels_v4( 1,nx, 1,ny, 1,nz,
                                  pipeline_rank, n_pipeline,
                                  &x, &y, &z );

  // Process bulk of voxels 4 at a time

# define LOAD_PTRS()    \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x+1,y,  z  ); \
  fy = &f(x,  y+1,z  ); \
  fz = &f(x,  y,  z+1)

# define ADVANCE_PTRS(n) \
  f0##n = f0++;          \
  fx##n = fx++;          \
  fy##n = fy++;          \
  fz##n = fz++;          \
  x++;                   \
  if( x>nx ) {           \
    x=1, y++;            \
    if( y>ny ) y=1, z++; \
    LOAD_PTRS();         \
  }

  LOAD_PTRS();
  ADVANCE_PTRS(0); ADVANCE_PTRS(1); ADVANCE_PTRS(2); ADVANCE_PTRS(3);

  for( ; n_voxel>3; n_voxel-=4 ) {

    load_4x1_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx, f0_cbx );
    load_4x2_tr( &f00->cby, &f01->cby, &f02->cby, &f03->cby, f0_cby, f0_cbz );

    load_4x1_tr( &fx0->cbx, &fx1->cbx, &fx2->cbx, &fx3->cbx, fx_cbx );
    load_4x1_tr( &fy0->cby, &fy1->cby, &fy2->cby, &fy3->cby, fy_cby );
    load_4x1_tr( &fz0->cbz, &fz1->cbz, &fz2->cbz, &fz3->cbz, fz_cbz );

    f0_div_b_err = fma( vpx,fx_cbx-f0_cbx, fma( vpy,fy_cby-f0_cby, vpz*(fz_cbz-f0_cbz) ) );

    store_4x1_tr( f0_div_b_err, &f00->div_b_err, &f01->div_b_err, &f02->div_b_err, &f03->div_b_err );

    ADVANCE_PTRS(0); ADVANCE_PTRS(1); ADVANCE_PTRS(2); ADVANCE_PTRS(3);

  }

  // Process remaining straggler voxels
  // Note that the stragglers have already been loaded above

  // FIXME: CONSIDER A MODEL WHERE THE HOST PROCESSES STRAGGLERS

  for( ; n_voxel; n_voxel-- ) {

    f00->div_b_err = px*(fx0->cbx-f00->cbx) + ( py*(fy0->cby-f00->cby) + pz*(fz0->cbz-f00->cbz) );

    f00 = f01; fx0 = fx1; fy0 = fy1; fz0 = fz1;
    f01 = f02; fx1 = fx2; fy1 = fy2; fz1 = fz2;

  }

}

#endif

