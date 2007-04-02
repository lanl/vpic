#if 0 // Original non-pipelined non-vectorized version 
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(2,y,z);
      fy = &f(1,y+1,z);
      fz = &f(1,y,z+1);
      for( x=1; x<=nx; x++ ) {
        f0->cbx -= py*( fy->ez-f0->ez ) - pz*( fz->ey-f0->ey )
        f0->cby -= pz*( fz->ex-f0->ex ) - px*( fx->ez-f0->ez )
        f0->cbz -= px*( fx->ey-f0->ey ) - py*( fy->ex-f0->ex )
	f0++; fx++; fy++; fz++;
      }
    }
  }
#endif

#include <field_pipelines.h>
#include CONCAT3(<,V4VERSION,>)

using namespace v4;

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
 
void
advance_b_pipeline( advance_b_pipeline_args_t * args,
		    int pipeline_rank ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;
  float frac               = args->frac;

  field_t * ALIGNED f0, * ALIGNED fx, * ALIGNED fy, * ALIGNED fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? frac*g->cvac*g->dt/g->dx : 0;
  const float py = (ny>1) ? frac*g->cvac*g->dt/g->dy : 0;
  const float pz = (nz>1) ? frac*g->cvac*g->dt/g->dz : 0;

  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);

  v4float f0_ex,  f0_ey,  f0_ez;  // Voxel quad electric fields
  v4float f0_cbx, f0_cby, f0_cbz; // Voxel quad magnetic fields
  v4float fx_ey, fx_ez;           // Voxel quad +x neighbor fields
  v4float fy_ez, fy_ex;           // Voxel quad +y neighbor fields
  v4float fz_ex, fz_ey;           // Voxel quad +z neighbor fields
  v4float dummy;

  field_t * ALIGNED f00, * ALIGNED f01, * ALIGNED f02, * ALIGNED f03; // Voxel quad
  field_t * ALIGNED fx0, * ALIGNED fx1, * ALIGNED fx2, * ALIGNED fx3; // Voxel quad +x neighbors
  field_t * ALIGNED fy0, * ALIGNED fy1, * ALIGNED fy2, * ALIGNED fy3; // Voxel quad +y neighbors
  field_t * ALIGNED fz0, * ALIGNED fz1, * ALIGNED fz2, * ALIGNED fz3; // Voxel quad +z neighbors
  
  // Process the voxels assigned to this pipeline 
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );
  
  // Process the bulk of the voxels 4 at a time

# define LOAD_PTRS()     \
  f0 = &f(x,  y,  z  );  \
  fx = &f(x+1,y,  z  );  \
  fy = &f(x,  y+1,z  );  \
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

    load_4x4_tr(  &f00->ex,  &f01->ex,  &f02->ex,  &f03->ex,  f0_ex,  f0_ey,  f0_ez,
                  /**/                                        f0_cbx                 );
    load_4x2_tr(  &f00->cby, &f01->cby, &f02->cby, &f03->cby,         f0_cby, f0_cbz );

    load_4x3_tr(  &fx0->ex,  &fx1->ex,  &fx2->ex,  &fx3->ex,  dummy, fx_ey, fx_ez );
    load_4x3_tr(  &fy0->ex,  &fy1->ex,  &fy2->ex,  &fy3->ex,  fy_ex, dummy, fy_ez );
    load_4x2_tr(  &fz0->ex,  &fz1->ex,  &fz2->ex,  &fz3->ex,  fz_ex, fz_ey  /**/  );

    f0_cbx -= fms( vpy,( fy_ez-f0_ez ),  vpz*( fz_ey-f0_ey ) );
    f0_cby -= fms( vpz,( fz_ex-f0_ex ),  vpx*( fx_ez-f0_ez ) );
    f0_cbz -= fms( vpx,( fx_ey-f0_ey ),  vpy*( fy_ex-f0_ex ) );

    store_4x1_tr( f0_cbx,         &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx );
    store_4x2_tr( f0_cby, f0_cbz, &f00->cby, &f01->cby, &f02->cby, &f03->cby );

    ADVANCE_PTRS(0); ADVANCE_PTRS(1); ADVANCE_PTRS(2); ADVANCE_PTRS(3);

  }

  // Process straggler voxels with scalar processing
  // Note that the pointers needed for straggler voxel processing were already
  // loaded implicitly above.

  // FIXME: CONSIDER A MODEL WHERE THE HOST PROCESSES STRAGGLERS.

  for( ; n_voxel; n_voxel-- ) {

    f00->cbx -= ( py*( fy0->ez-f00->ez ) - pz*( fz0->ey-f00->ey ) );
    f00->cby -= ( pz*( fz0->ex-f00->ex ) - px*( fx0->ez-f00->ez ) );
    f00->cbz -= ( px*( fx0->ey-f00->ey ) - py*( fy0->ex-f00->ex ) );

    f00 = f01; fx0 = fx1; fy0 = fy1; fz0 = fz1;
    f01 = f02; fx1 = fx2; fy1 = fy2; fz1 = fz2;

  }

}

