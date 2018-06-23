using namespace v8;

void
clean_div_b_pipeline_v8( pipeline_args_t * args,
                         int pipeline_rank,
                         int n_pipeline )
{
  field_t      * ALIGNED(128) f = args->f;
  const grid_t *              g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;
  
  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  float px, py, pz, alphadt;

  px = ( nx > 1 ) ? g->rdx : 0;
  py = ( ny > 1 ) ? g->rdy : 0;
  pz = ( nz > 1 ) ? g->rdz : 0;

  alphadt = 0.3888889/( px*px + py*py + pz*pz );

  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;

  const v8float vpx(px);
  const v8float vpy(py);
  const v8float vpz(pz);

  v8float f0_cbx, f0_cby, f0_cbz; // Voxel quad magnetic fields
  v8float f0_div_b_err;           // Voxel quad div b errs
  v8float fx_div_b_err;           // Voxel quad -x neighbor div b err
  v8float fy_div_b_err;           // Voxel quad -y neighbor div b err
  v8float fz_div_b_err;           // Voxel quad -z neighbor div b err

  field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel quad
  field_t * ALIGNED(16) f04, * ALIGNED(16) f05, * ALIGNED(16) f06, * ALIGNED(16) f07; // Voxel quad

  field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel quad +x neighbors
  field_t * ALIGNED(16) fx4, * ALIGNED(16) fx5, * ALIGNED(16) fx6, * ALIGNED(16) fx7; // Voxel quad +x neighbors

  field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel quad +y neighbors
  field_t * ALIGNED(16) fy4, * ALIGNED(16) fy5, * ALIGNED(16) fy6, * ALIGNED(16) fy7; // Voxel quad +y neighbors

  field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel quad +z neighbors
  field_t * ALIGNED(16) fz4, * ALIGNED(16) fz5, * ALIGNED(16) fz6, * ALIGNED(16) fz7; // Voxel quad +z neighbors

  // Process voxels assigned to this pipeline 
  
  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  // Process bulk of voxels 8 at a time

# define LOAD_STENCIL()     \
  f0 = &f( x,   y,   z   ); \
  fx = &f( x-1, y,   z   ); \
  fy = &f( x,   y-1, z   ); \
  fz = &f( x,   y,   z-1 )

# define NEXT_STENCIL(n)      \
  f0##n = f0++;               \
  fx##n = fx++;               \
  fy##n = fy++;               \
  fz##n = fz++;               \
  x++;                        \
  if ( x > nx )               \
  {			      \
                  x = 2, y++; \
    if ( y > ny ) y = 2, z++; \
    LOAD_STENCIL();           \
  }

  LOAD_STENCIL();

  for( ; n_voxel > 7; n_voxel -= 8 )
  {
    NEXT_STENCIL(0);
    NEXT_STENCIL(1);
    NEXT_STENCIL(2);
    NEXT_STENCIL(3);
    NEXT_STENCIL(4);
    NEXT_STENCIL(5);
    NEXT_STENCIL(6);
    NEXT_STENCIL(7);

    load_8x4_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx,
		 &f04->cbx, &f05->cbx, &f06->cbx, &f07->cbx,
		 f0_cbx, f0_cby, f0_cbz, f0_div_b_err );

    fx_div_b_err = v8float( fx0->div_b_err, fx1->div_b_err, fx2->div_b_err, fx3->div_b_err,
			    fx4->div_b_err, fx5->div_b_err, fx6->div_b_err, fx7->div_b_err );

    fy_div_b_err = v8float( fy0->div_b_err, fy1->div_b_err, fy2->div_b_err, fy3->div_b_err,
			    fy4->div_b_err, fy5->div_b_err, fy6->div_b_err, fy7->div_b_err );

    fz_div_b_err = v8float( fz0->div_b_err, fz1->div_b_err, fz2->div_b_err, fz3->div_b_err,
			    fz4->div_b_err, fz5->div_b_err, fz6->div_b_err, fz7->div_b_err );

    f0_cbx = fma( f0_div_b_err - fx_div_b_err, px, f0_cbx );
    f0_cby = fma( f0_div_b_err - fy_div_b_err, py, f0_cby );
    f0_cbz = fma( f0_div_b_err - fz_div_b_err, pz, f0_cbz );

    store_8x4_tr( f0_cbx, f0_cby, f0_cbz, f0_div_b_err,
		  &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx,
		  &f04->cbx, &f05->cbx, &f06->cbx, &f07->cbx );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL
}