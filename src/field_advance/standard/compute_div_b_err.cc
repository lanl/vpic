#define IN_sfa
#include "sfa_private.h"

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

typedef struct pipeline_args {
  field_t      * ALIGNED(128) f;
  const grid_t *              g;
} pipeline_args_t;

void
compute_div_b_err_pipeline( pipeline_args_t * args,
                            int pipeline_rank,
                            int n_pipeline ) {
  field_t      * ALIGNED(128) f = args->f;
  const grid_t *              g = args->g;
  
  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? g->rdx : 0;
  const float py = (ny>1) ? g->rdy : 0;
  const float pz = (nz>1) ? g->rdz : 0;

  // Process the voxels assigned to this pipeline
  
  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x+1,y,  z  ); \
  fy = &f(x,  y+1,z  ); \
  fz = &f(x,  y,  z+1)

  LOAD_STENCIL();

  for( ; n_voxel; n_voxel-- ) {
    f0->div_b_err = px*( fx->cbx - f0->cbx ) +
                    py*( fy->cby - f0->cby ) +
                    pz*( fz->cbz - f0->cbz );
    f0++; fx++; fy++; fz++;
    
    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      LOAD_STENCIL();
    }
  }

# undef LOAD_STENCIL

}

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

using namespace v4;

static void
compute_div_b_err_pipeline_v4( pipeline_args_t * args,
                               int pipeline_rank,
                               int n_pipeline ) {
  field_t      * ALIGNED(128) f = args->f;
  const grid_t *              g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? g->rdx : 0;
  const float py = (ny>1) ? g->rdy : 0;
  const float pz = (nz>1) ? g->rdz : 0;

  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);

  v4float f0_cbx, f0_cby, f0_cbz; // Voxel quad magnetic fields
  v4float f0_div_b_err;           // Voxel quad div b errs
  v4float fx_cbx;                 // Voxel quad +x neighbor x magnetic fields
  v4float fy_cby;                 // Voxel quad +y neighbor y magnetic fields
  v4float fz_cbz;                 // Voxel quad +z neighbor z magnetic fields

  field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel quad
  field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel quad +x neighbors
  field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel quad +x neighbors
  field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel quad +x neighbors

  // Process the voxels assigned to this pipeline 
  
  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  // Process bulk of voxels 4 at a time

# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x+1,y,  z  ); \
  fy = &f(x,  y+1,z  ); \
  fz = &f(x,  y,  z+1)

# define NEXT_STENCIL(n) \
  f0##n = f0++;          \
  fx##n = fx++;          \
  fy##n = fy++;          \
  fz##n = fz++;          \
  x++;                   \
  if( x>nx ) {           \
    x=1, y++;            \
    if( y>ny ) y=1, z++; \
    LOAD_STENCIL();      \
  }

  LOAD_STENCIL();

  for( ; n_voxel>3; n_voxel-=4 ) {
    NEXT_STENCIL(0); NEXT_STENCIL(1); NEXT_STENCIL(2); NEXT_STENCIL(3);

    load_4x3_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx, f0_cbx, f0_cby, f0_cbz );

    fx_cbx = v4float( fx0->cbx, fx1->cbx, fx2->cbx, fx3->cbx );
    fy_cby = v4float( fy0->cby, fy1->cby, fy2->cby, fy3->cby );
    fz_cbz = v4float( fz0->cbz, fz1->cbz, fz2->cbz, fz3->cbz );

    f0_div_b_err = fma( vpx,fx_cbx-f0_cbx, fma( vpy,fy_cby-f0_cby, vpz*(fz_cbz-f0_cbz) ) );

    store_4x1_tr( f0_div_b_err, &f00->div_b_err, &f01->div_b_err, &f02->div_b_err, &f03->div_b_err );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL

}

#endif

#if defined(V8_ACCELERATION) && defined(HAS_V8_PIPELINE)

using namespace v8;

static void
compute_div_b_err_pipeline_v8( pipeline_args_t * args,
                               int pipeline_rank,
                               int n_pipeline ) {
  field_t      * ALIGNED(128) f = args->f;
  const grid_t *              g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? g->rdx : 0;
  const float py = (ny>1) ? g->rdy : 0;
  const float pz = (nz>1) ? g->rdz : 0;

  const v8float vpx(px);
  const v8float vpy(py);
  const v8float vpz(pz);

  v8float f0_cbx, f0_cby, f0_cbz; // Voxel quad magnetic fields
  v8float f0_div_b_err;           // Voxel quad div b errs
  v8float fx_cbx;                 // Voxel quad +x neighbor x magnetic fields
  v8float fy_cby;                 // Voxel quad +y neighbor y magnetic fields
  v8float fz_cbz;                 // Voxel quad +z neighbor z magnetic fields

  field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel quad
  field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel quad +x neighbors
  field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel quad +x neighbors
  field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel quad +x neighbors

  // Process the voxels assigned to this pipeline 
  
  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  // Process bulk of voxels 4 at a time

# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x+1,y,  z  ); \
  fy = &f(x,  y+1,z  ); \
  fz = &f(x,  y,  z+1)

# define NEXT_STENCIL(n) \
  f0##n = f0++;          \
  fx##n = fx++;          \
  fy##n = fy++;          \
  fz##n = fz++;          \
  x++;                   \
  if( x>nx ) {           \
    x=1, y++;            \
    if( y>ny ) y=1, z++; \
    LOAD_STENCIL();      \
  }

  LOAD_STENCIL();

  for( ; n_voxel>3; n_voxel-=4 ) {
    NEXT_STENCIL(0); NEXT_STENCIL(1); NEXT_STENCIL(2); NEXT_STENCIL(3);

    load_8x3_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx, f0_cbx, f0_cby, f0_cbz );

    fx_cbx = v8float( fx0->cbx, fx1->cbx, fx2->cbx, fx3->cbx );
    fy_cby = v8float( fy0->cby, fy1->cby, fy2->cby, fy3->cby );
    fz_cbz = v8float( fz0->cbz, fz1->cbz, fz2->cbz, fz3->cbz );

    f0_div_b_err = fma( vpx,fx_cbx-f0_cbx, fma( vpy,fy_cby-f0_cby, vpz*(fz_cbz-f0_cbz) ) );

    store_8x1_tr( f0_div_b_err, &f00->div_b_err, &f01->div_b_err, &f02->div_b_err, &f03->div_b_err );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL

}

#endif

void
compute_div_b_err( field_array_t * RESTRICT fa ) {
  pipeline_args_t args[1];

  if( !fa ) ERROR(( "Bad args" ));
  
# if 0 // Original non-pipelined version
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
# endif

  args->f = fa->f;
  args->g = fa->g;
  EXEC_PIPELINES( compute_div_b_err, args, 0 );
  WAIT_PIPELINES();
}
