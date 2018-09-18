using namespace v4;

void
clean_div_b_pipeline_v4( pipeline_args_t * args,
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

  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);

  v4float f0_cbx, f0_cby, f0_cbz; // Voxel block magnetic fields
  v4float f0_div_b_err;           // Voxel block div b errs
  v4float fx_div_b_err;           // Voxel block -x neighbor div b err
  v4float fy_div_b_err;           // Voxel block -y neighbor div b err
  v4float fz_div_b_err;           // Voxel block -z neighbor div b err

  field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel block

  field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel block +x neighbors

  field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel block +y neighbors

  field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel block +z neighbors

  // Process voxels assigned to this pipeline 
  
  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  // Process bulk of voxels 4 at a time

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

  for( ; n_voxel > 3; n_voxel -= 4 )
  {
    NEXT_STENCIL(0);
    NEXT_STENCIL(1);
    NEXT_STENCIL(2);
    NEXT_STENCIL(3);

    load_4x4_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx,
                 f0_cbx, f0_cby, f0_cbz, f0_div_b_err );

    fx_div_b_err = v4float( fx0->div_b_err, fx1->div_b_err, fx2->div_b_err, fx3->div_b_err );

    fy_div_b_err = v4float( fy0->div_b_err, fy1->div_b_err, fy2->div_b_err, fy3->div_b_err );

    fz_div_b_err = v4float( fz0->div_b_err, fz1->div_b_err, fz2->div_b_err, fz3->div_b_err );

    f0_cbx = fma( f0_div_b_err - fx_div_b_err, px, f0_cbx );
    f0_cby = fma( f0_div_b_err - fy_div_b_err, py, f0_cby );
    f0_cbz = fma( f0_div_b_err - fz_div_b_err, pz, f0_cbz );

    store_4x4_tr( f0_cbx, f0_cby, f0_cbz, f0_div_b_err,
                  &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL
}
