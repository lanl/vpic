using namespace v8;

void
advance_b_pipeline_v8( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline )
{
  DECLARE_STENCIL();

  int n_voxel;

  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  const v8float vpx( px );
  const v8float vpy( py );
  const v8float vpz( pz );

  v8float f0_ex,  f0_ey,  f0_ez;  // Voxel block electric fields
  v8float f0_cbx, f0_cby, f0_cbz; // Voxel block magnetic fields
  v8float fx_ey, fx_ez;           // Voxel block +x neighbor fields
  v8float fy_ez, fy_ex;           // Voxel block +y neighbor fields
  v8float fz_ex, fz_ey;           // Voxel block +z neighbor fields
  v8float dummy;

  field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel block
  field_t * ALIGNED(16) f04, * ALIGNED(16) f05, * ALIGNED(16) f06, * ALIGNED(16) f07; // Voxel block

  field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel block +x neighbors
  field_t * ALIGNED(16) fx4, * ALIGNED(16) fx5, * ALIGNED(16) fx6, * ALIGNED(16) fx7; // Voxel block +x neighbors

  field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel block +y neighbors
  field_t * ALIGNED(16) fy4, * ALIGNED(16) fy5, * ALIGNED(16) fy6, * ALIGNED(16) fy7; // Voxel block +y neighbors

  field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel block +z neighbors
  field_t * ALIGNED(16) fz4, * ALIGNED(16) fz5, * ALIGNED(16) fz6, * ALIGNED(16) fz7; // Voxel block +z neighbors

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

    load_8x3_tr( &f00->ex, &f01->ex, &f02->ex, &f03->ex,
		 &f04->ex, &f05->ex, &f06->ex, &f07->ex,
		 f0_ex, f0_ey, f0_ez );

    load_8x3_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx,
		 &f04->cbx, &f05->cbx, &f06->cbx, &f07->cbx,
		 f0_cbx, f0_cby, f0_cbz );

    load_8x3_tr( &fx0->ex, &fx1->ex, &fx2->ex, &fx3->ex,
		 &fx4->ex, &fx5->ex, &fx6->ex, &fx7->ex,
		 dummy, fx_ey, fx_ez );

    load_8x3_tr( &fy0->ex, &fy1->ex, &fy2->ex, &fy3->ex,
		 &fy4->ex, &fy5->ex, &fy6->ex, &fy7->ex,
		 fy_ex, dummy, fy_ez );

    load_8x2_tr( &fz0->ex, &fz1->ex, &fz2->ex, &fz3->ex,
		 &fz4->ex, &fz5->ex, &fz6->ex, &fz7->ex,
		 fz_ex, fz_ey );

    f0_cbx += fnms( vpy, ( fy_ez - f0_ez ), vpz*( fz_ey - f0_ey ) );
    f0_cby += fnms( vpz, ( fz_ex - f0_ex ), vpx*( fx_ez - f0_ez ) );
    f0_cbz += fnms( vpx, ( fx_ey - f0_ey ), vpy*( fy_ex - f0_ex ) );

    store_8x3_tr( f0_cbx, f0_cby, f0_cbz,
		  &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx,
		  &f04->cbx, &f05->cbx, &f06->cbx, &f07->cbx );
  }
}
