#define IN_sfa

#include "sfa_private.h"

#if defined(V8_ACCELERATION)

using namespace v8;

void
vacuum_advance_e_pipeline_v8( pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline )
{
  DECLARE_STENCIL();

  int n_voxel;

  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  const v8float vdecayx( decayx ), vdrivex( drivex );
  const v8float vdecayy( decayy ), vdrivey( drivey );
  const v8float vdecayz( decayz ), vdrivez( drivez );
  const v8float vdamp( damp );
  const v8float vpx_muz( px_muz ), vpx_muy( px_muy );
  const v8float vpy_mux( py_mux ), vpy_muz( py_muz );
  const v8float vpz_muy( pz_muy ), vpz_mux( pz_mux );
  const v8float vcj( cj );

  v8float save0, save1, dummy;

  v8float f0_ex,   f0_ey,   f0_ez;
  v8float f0_cbx,  f0_cby,  f0_cbz;
  v8float f0_tcax, f0_tcay, f0_tcaz;
  v8float f0_jfx,  f0_jfy,  f0_jfz;
  v8float          fx_cby,  fx_cbz;
  v8float fy_cbx,           fy_cbz;
  v8float fz_cbx,  fz_cby;

  field_t * ALIGNED(32) f00, * ALIGNED(32) f01, * ALIGNED(32) f02, * ALIGNED(32) f03; // Voxel block
  field_t * ALIGNED(32) f04, * ALIGNED(32) f05, * ALIGNED(32) f06, * ALIGNED(32) f07; // Voxel block

  field_t * ALIGNED(32) fx0, * ALIGNED(32) fx1, * ALIGNED(32) fx2, * ALIGNED(32) fx3; // Voxel block +x neighbors
  field_t * ALIGNED(32) fx4, * ALIGNED(32) fx5, * ALIGNED(32) fx6, * ALIGNED(32) fx7; // Voxel block +x neighbors

  field_t * ALIGNED(32) fy0, * ALIGNED(32) fy1, * ALIGNED(32) fy2, * ALIGNED(32) fy3; // Voxel block +y neighbors
  field_t * ALIGNED(32) fy4, * ALIGNED(32) fy5, * ALIGNED(32) fy6, * ALIGNED(32) fy7; // Voxel block +y neighbors

  field_t * ALIGNED(32) fz0, * ALIGNED(32) fz1, * ALIGNED(32) fz2, * ALIGNED(32) fz3; // Voxel block +z neighbors
  field_t * ALIGNED(32) fz4, * ALIGNED(32) fz5, * ALIGNED(32) fz6, * ALIGNED(32) fz7; // Voxel block +z neighbors

  // Process the bulk of the voxels 8 at a time

  INIT_STENCIL();

  for( ; n_voxel > 7; n_voxel -= 8 )
  {
    f00 = f0; fx0 = fx; fy0 = fy; fz0 = fz; NEXT_STENCIL();
    f01 = f0; fx1 = fx; fy1 = fy; fz1 = fz; NEXT_STENCIL();
    f02 = f0; fx2 = fx; fy2 = fy; fz2 = fz; NEXT_STENCIL();
    f03 = f0; fx3 = fx; fy3 = fy; fz3 = fz; NEXT_STENCIL();
    f04 = f0; fx4 = fx; fy4 = fy; fz4 = fz; NEXT_STENCIL();
    f05 = f0; fx5 = fx; fy5 = fy; fz5 = fz; NEXT_STENCIL();
    f06 = f0; fx6 = fx; fy6 = fy; fz6 = fz; NEXT_STENCIL();
    f07 = f0; fx7 = fx; fy7 = fy; fz7 = fz; NEXT_STENCIL();

    //------------------------------------------------------------------------//
    // Load field data.
    //------------------------------------------------------------------------//

    load_8x4_tr( &f00->ex, &f01->ex, &f02->ex, &f03->ex,
                 &f04->ex, &f05->ex, &f06->ex, &f07->ex,
                 f0_ex, f0_ey, f0_ez, save0 );

    load_8x3_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx,
                 &f04->cbx, &f05->cbx, &f06->cbx, &f07->cbx,
                 f0_cbx, f0_cby, f0_cbz );

    load_8x4_tr( &f00->tcax, &f01->tcax, &f02->tcax, &f03->tcax,
                 &f04->tcax, &f05->tcax, &f06->tcax, &f07->tcax,
                 f0_tcax, f0_tcay, f0_tcaz, save1 );

    load_8x3_tr( &f00->jfx, &f01->jfx, &f02->jfx, &f03->jfx,
                 &f04->jfx, &f05->jfx, &f06->jfx, &f07->jfx,
                 f0_jfx, f0_jfy, f0_jfz );

    load_8x3_tr( &fx0->cbx, &fx1->cbx, &fx2->cbx, &fx3->cbx,
                 &fx4->cbx, &fx5->cbx, &fx6->cbx, &fx7->cbx,
                 dummy, fx_cby, fx_cbz );

    load_8x3_tr( &fy0->cbx, &fy1->cbx, &fy2->cbx, &fy3->cbx,
                 &fy4->cbx, &fy5->cbx, &fy6->cbx, &fy7->cbx,
                 fy_cbx, dummy, fy_cbz );

    load_8x2_tr( &fz0->cbx, &fz1->cbx, &fz2->cbx, &fz3->cbx,
                 &fz4->cbx, &fz5->cbx, &fz6->cbx, &fz7->cbx,
                 fz_cbx, fz_cby );

    f0_tcax = fnms( vdamp,
                    f0_tcax,
                    fms( vpy_muz,
                         ( f0_cbz - fy_cbz ),
                         vpz_muy * ( f0_cby - fz_cby ) ) );

    f0_tcay = fnms( vdamp,
                    f0_tcay,
                    fms( vpz_mux,
                         ( f0_cbx - fz_cbx ),
                         vpx_muz * ( f0_cbz - fx_cbz ) ) );

    f0_tcaz = fnms( vdamp,
                    f0_tcaz,
                    fms( vpx_muy,
                         ( f0_cby - fx_cby ),
                         vpy_mux * ( f0_cbx - fy_cbx ) ) );

    f0_ex   = fma( vdecayx, f0_ex, vdrivex * fnms( vcj, f0_jfx, f0_tcax ) );

    f0_ey   = fma( vdecayy, f0_ey, vdrivey * fnms( vcj, f0_jfy, f0_tcay ) );

    f0_ez   = fma( vdecayz, f0_ez, vdrivez * fnms( vcj, f0_jfz, f0_tcaz ) );

    //------------------------------------------------------------------------//
    // Note:
    //------------------------------------------------------------------------//
    // Unlike load_8x3 versus load_8x4, store_8x4 is much more efficient
    // than store_8x3.
    //------------------------------------------------------------------------//

    store_8x4_tr( f0_ex, f0_ey, f0_ez, save0,
                  &f00->ex, &f01->ex, &f02->ex, &f03->ex,
                  &f04->ex, &f05->ex, &f06->ex, &f07->ex );

    store_8x4_tr( f0_tcax, f0_tcay, f0_tcaz, save1,
                  &f00->tcax, &f01->tcax, &f02->tcax, &f03->tcax,
                  &f04->tcax, &f05->tcax, &f06->tcax, &f07->tcax );
  }
}

#else

void
vacuum_advance_e_pipeline_v8( pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline )
{
  // No v8 implementation.
  ERROR( ( "No vacuum_advance_e_pipeline_v8 implementation." ) );
}

#endif
