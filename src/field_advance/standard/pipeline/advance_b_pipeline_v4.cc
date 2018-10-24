#define IN_sfa
#define IN_advance_b_pipeline

#include "advance_b_pipeline.h"

#include "../sfa_private.h"

#if defined(V4_ACCELERATION)

using namespace v4;

void
advance_b_pipeline_v4( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline )
{
  DECLARE_STENCIL();

  int n_voxel;

  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  const v4float vpx( px );
  const v4float vpy( py );
  const v4float vpz( pz );

  v4float f0_ex,  f0_ey,  f0_ez;  // Voxel quad electric fields
  v4float f0_cbx, f0_cby, f0_cbz; // Voxel quad magnetic fields
  v4float fx_ey, fx_ez;           // Voxel quad +x neighbor fields
  v4float fy_ez, fy_ex;           // Voxel quad +y neighbor fields
  v4float fz_ex, fz_ey;           // Voxel quad +z neighbor fields
  v4float dummy;

  field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel quad

  field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel quad +x neighbors

  field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel quad +y neighbors

  field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel quad +z neighbors

  // Process the bulk of the voxels 4 at a time

  INIT_STENCIL();

  for( ; n_voxel > 3; n_voxel -= 4 )
  {
    f00 = f0; fx0 = fx; fy0 = fy; fz0 = fz; NEXT_STENCIL();
    f01 = f0; fx1 = fx; fy1 = fy; fz1 = fz; NEXT_STENCIL();
    f02 = f0; fx2 = fx; fy2 = fy; fz2 = fz; NEXT_STENCIL();
    f03 = f0; fx3 = fx; fy3 = fy; fz3 = fz; NEXT_STENCIL();

    load_4x3_tr( &f00->ex, &f01->ex, &f02->ex, &f03->ex,
		 f0_ex, f0_ey, f0_ez );

    load_4x3_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx,
		 f0_cbx, f0_cby, f0_cbz );

    load_4x3_tr( &fx0->ex, &fx1->ex, &fx2->ex, &fx3->ex,
		 dummy, fx_ey, fx_ez );

    load_4x3_tr( &fy0->ex, &fy1->ex, &fy2->ex, &fy3->ex,
		 fy_ex, dummy, fy_ez );

    load_4x2_tr( &fz0->ex, &fz1->ex, &fz2->ex, &fz3->ex,
		 fz_ex, fz_ey );

    f0_cbx += fnms( vpy, ( fy_ez - f0_ez ), vpz*( fz_ey - f0_ey ) );
    f0_cby += fnms( vpz, ( fz_ex - f0_ex ), vpx*( fx_ez - f0_ez ) );
    f0_cbz += fnms( vpx, ( fx_ey - f0_ey ), vpy*( fy_ex - f0_ex ) );

    store_4x3_tr( f0_cbx, f0_cby, f0_cbz,
		  &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx );
  }
}

#else

void
advance_b_pipeline_v4( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline )
{
  // No v4 implementation.
  ERROR( ( "No advance_b_pipeline_v4 implementation." ) );
}

#endif
