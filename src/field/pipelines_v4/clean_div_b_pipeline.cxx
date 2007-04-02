#if 0 // Original non-pipelined non-vectorized version
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
        f0->cbx += px*( f0->div_b_err - fx->div_b_err )
        f0->cby += py*( f0->div_b_err - fy->div_b_err )
        f0->cbz += pz*( f0->div_b_err - fz->div_b_err )
	f0++; fx++; fy++; fz++;
      }
    }
  }
#endif

#include <field_pipelines.h>
#include <v4.hxx>

using namespace v4;

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
clean_div_b_pipeline( clean_div_b_pipeline_args_t * args,
                      int pipeline_rank ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;

  field_t *f0, *fx, *fy, *fz;
  int x, y, z, n_voxel;
  
  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  float px, py, pz, alphadt;

  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;
  alphadt = 0.3888889/( px*px + py*py + pz*pz );
  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;

  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);

  v4float f0_cbx, f0_cby, f0_cbz; // Voxel quad magnetic fields
  v4float f0_div_b_err;           // Voxel quad div b errs
  v4float fx_div_b_err;           // Voxel quad -x neighbor div b err
  v4float fy_div_b_err;           // Voxel quad -y neighbor div b err
  v4float fz_div_b_err;           // Voxel quad -z neighbor div b err

  field_t * ALIGNED f00, * ALIGNED f01, * ALIGNED f02, * ALIGNED f03; // Voxel quad
  field_t * ALIGNED fx0, * ALIGNED fx1, * ALIGNED fx2, * ALIGNED fx3; // Voxel quad +x neighbors
  field_t * ALIGNED fy0, * ALIGNED fy1, * ALIGNED fy2, * ALIGNED fy3; // Voxel quad +x neighbors
  field_t * ALIGNED fz0, * ALIGNED fz1, * ALIGNED fz2, * ALIGNED fz3; // Voxel quad +x neighbors

  // Process voxels assigned to this pipeline 
  
  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  // Process bulk of voxels 4 at a time

# define LOAD_PTRS()    \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x-1,y,  z  ); \
  fy = &f(x,  y-1,z  ); \
  fz = &f(x,  y,  z-1)

# define ADVANCE_PTRS(n) \
  f0##n = f0++;          \
  fx##n = fx++;          \
  fy##n = fy++;          \
  fz##n = fz++;          \
  x++;                   \
  if( x>nx ) {           \
    x=2, y++;            \
    if( y>ny ) y=2, z++; \
    LOAD_PTRS();         \
  }

  LOAD_PTRS();
  ADVANCE_PTRS(0); ADVANCE_PTRS(1); ADVANCE_PTRS(2); ADVANCE_PTRS(3);

  for( ; n_voxel>3; n_voxel-=4 ) {

    load_4x1_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx, f0_cbx );
    load_4x2_tr( &f00->cby, &f01->cby, &f02->cby, &f03->cby, f0_cby, f0_cbz );

    load_4x1_tr( &f00->div_b_err, &f01->div_b_err, &f02->div_b_err, &f03->div_b_err, f0_div_b_err );
    load_4x1_tr( &fx0->div_b_err, &fx1->div_b_err, &fx2->div_b_err, &fx3->div_b_err, fx_div_b_err );
    load_4x1_tr( &fy0->div_b_err, &fy1->div_b_err, &fy2->div_b_err, &fy3->div_b_err, fy_div_b_err );
    load_4x1_tr( &fz0->div_b_err, &fz1->div_b_err, &fz2->div_b_err, &fz3->div_b_err, fz_div_b_err );

    f0_cbx = fma( f0_div_b_err-fx_div_b_err, px, f0_cbx );
    f0_cby = fma( f0_div_b_err-fy_div_b_err, py, f0_cby );
    f0_cbz = fma( f0_div_b_err-fz_div_b_err, pz, f0_cbz );

    store_4x1_tr( f0_cbx,         &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx );
    store_4x2_tr( f0_cby, f0_cbz, &f00->cby, &f01->cby, &f02->cby, &f03->cby );

    ADVANCE_PTRS(0); ADVANCE_PTRS(1); ADVANCE_PTRS(2); ADVANCE_PTRS(3);

  }
  
  // Process remaining straggler voxels
  // Note that the stragglers have already been loaded above

  // FIXME: CONSIDER A MODEL WHERE THE HOST PROCESSES STRAGGLERS

  for( ; n_voxel; n_voxel-- ) {

    f00->cbx += px*( f00->div_b_err-fx0->div_b_err );
    f00->cby += py*( f00->div_b_err-fy0->div_b_err );
    f00->cbz += pz*( f00->div_b_err-fz0->div_b_err );

    f00 = f01; fx0 = fx1; fy0 = fy1; fz0 = fz1;
    f01 = f02; fx1 = fx2; fy1 = fy2; fz1 = fz2;

  }

}

