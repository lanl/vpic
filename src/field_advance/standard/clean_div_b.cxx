#define IN_sfa
#define HAS_V4_PIPELINE
#include "sfa_private.h"

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

#define MARDER_CBX() f0->cbx += px*( f0->div_b_err - fx->div_b_err )
#define MARDER_CBY() f0->cby += py*( f0->div_b_err - fy->div_b_err )
#define MARDER_CBZ() f0->cbz += pz*( f0->div_b_err - fz->div_b_err )

typedef struct pipeline_args {
  field_t      * ALIGNED(128) f;
  const grid_t *              g;
} pipeline_args_t;

void
clean_div_b_pipeline( pipeline_args_t * args,
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

  float px, py, pz, alphadt;

  px = (nx>1) ? g->rdx : 0;
  py = (ny>1) ? g->rdy : 0;
  pz = (nz>1) ? g->rdz : 0;
  alphadt = 0.3888889/( px*px + py*py + pz*pz );
  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;

  // Process voxels assigned to this pipeline
  
  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x-1,y,  z  ); \
  fy = &f(x,  y-1,z  ); \
  fz = &f(x,  y,  z-1)

  LOAD_STENCIL();
  
  for( ; n_voxel; n_voxel-- ) {
    MARDER_CBX();
    MARDER_CBY();
    MARDER_CBZ();
    f0++; fx++; fy++; fz++;
    
    x++;
    if( x>nx ) {
      x=2, y++;
      if( y>ny ) y=2, z++;
      LOAD_STENCIL();
    }      
  }

# undef LOAD_STENCIL

}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

// SPU pipeline is defined in a different compilation unit

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

using namespace v4;

void
clean_div_b_pipeline_v4( pipeline_args_t * args,
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

  float px, py, pz, alphadt;

  px = (nx>1) ? g->rdx : 0;
  py = (ny>1) ? g->rdy : 0;
  pz = (nz>1) ? g->rdz : 0;
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

  field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel quad
  field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel quad +x neighbors
  field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel quad +x neighbors
  field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel quad +x neighbors

  // Process voxels assigned to this pipeline 
  
  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  // Process bulk of voxels 4 at a time

# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x-1,y,  z  ); \
  fy = &f(x,  y-1,z  ); \
  fz = &f(x,  y,  z-1)

# define NEXT_STENCIL(n) \
  f0##n = f0++;          \
  fx##n = fx++;          \
  fy##n = fy++;          \
  fz##n = fz++;          \
  x++;                   \
  if( x>nx ) {           \
    x=2, y++;            \
    if( y>ny ) y=2, z++; \
    LOAD_STENCIL();      \
  }

  LOAD_STENCIL();

  for( ; n_voxel>3; n_voxel-=4 ) {
    NEXT_STENCIL(0); NEXT_STENCIL(1); NEXT_STENCIL(2); NEXT_STENCIL(3);

    load_4x4_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx, f0_cbx, f0_cby, f0_cbz, f0_div_b_err );

    fx_div_b_err = v4float( fx0->div_b_err, fx1->div_b_err, fx2->div_b_err, fx3->div_b_err );
    fy_div_b_err = v4float( fy0->div_b_err, fy1->div_b_err, fy2->div_b_err, fy3->div_b_err );
    fz_div_b_err = v4float( fz0->div_b_err, fz1->div_b_err, fz2->div_b_err, fz3->div_b_err );

    f0_cbx = fma( f0_div_b_err-fx_div_b_err, px, f0_cbx );
    f0_cby = fma( f0_div_b_err-fy_div_b_err, py, f0_cby );
    f0_cbz = fma( f0_div_b_err-fz_div_b_err, pz, f0_cbz );

    store_4x4_tr( f0_cbx, f0_cby, f0_cbz, f0_div_b_err, &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL

}

#endif

void
clean_div_b( field_array_t * fa ) {
  pipeline_args_t args[1];
  
  field_t * f, * f0, * fx, * fy, * fz;
  const grid_t * g;
  float alphadt, px, py, pz;
  int x, y, z, nx, ny, nz;

  if( !fa ) ERROR(( "Bad args" ));
  f = fa->f;
  g = fa->g;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? g->rdx : 0;
  py = (ny>1) ? g->rdy : 0;
  pz = (nz>1) ? g->rdz : 0;
  alphadt = 0.3888889/( px*px + py*py + pz*pz );
  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;

  // Have pipelines do Marder pass in interior.  The host handles
  // stragglers.

# if 0 // Original non-pipelined version
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
	MARDER_CBX();
	MARDER_CBY();
	MARDER_CBZ();
	f0++; fx++; fy++; fz++;
      }
    }
  }
# endif

  // Begin setting derr ghosts
  begin_remote_ghost_div_b( f, g );
  local_ghost_div_b( f, g);

  // Have pipelines do interior of the local domain
  args->f = f;
  args->g = g;
  EXEC_PIPELINES( clean_div_b, args, 0 );
  
  // Do left over interior bx
  for( y=1; y<=ny; y++ ) {
    f0 = &f(2,y,1);
    fx = &f(1,y,1);
    for( x=2; x<=nx; x++ ) {
      MARDER_CBX();
      f0++;
      fx++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(2,1,z);
    fx = &f(1,1,z);
    for( x=2; x<=nx; x++ ) {
      MARDER_CBX();
      f0++;
      fx++;
    }
  }

  // Left over interior by
  for( z=1; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fy = &f(1,y-1,z);
      MARDER_CBY();
    }
  }
  for( y=2; y<=ny; y++ ) {
    f0 = &f(2,y,  1);
    fy = &f(2,y-1,1);
    for( x=2; x<=nx; x++ ) {
      MARDER_CBY();
      f0++;
      fy++;
    }
  }

  // Left over interior bz
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fz = &f(1,1,z-1);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBZ();
      f0++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fz = &f(1,y,z-1);
      MARDER_CBZ();
    }
  }

  // Finish setting derr ghosts
  
  end_remote_ghost_div_b( f, g );

  // Do Marder pass in exterior

  // Exterior bx
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(0,y,z);
      MARDER_CBX();
    }
  }
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,z);
      fx = &f(nx,  y,z);
      MARDER_CBX();
    }
  }

  // Exterior by
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fy = &f(1,0,z);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBY();
      f0++;
      fy++;
    }
  }
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fy = &f(1,ny,  z);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBY();
      f0++;
      fy++;
    }
  }

  // Exterior bz
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,1);
    fz = &f(1,y,0);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBZ();
      f0++;
      fz++;
    }
  }
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,nz+1);
    fz = &f(1,y,nz);
    for( x=1; x<=nx; x++ ) {
      MARDER_CBZ();
      f0++;
      fz++;
    }
  }

  // Wait for pipelines to finish up cleaning div_b in interior
  
  WAIT_PIPELINES();
  
  local_adjust_norm_b(f,g);
}
