#define IN_sfa
#include "sfa_private.h"

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define MARDER_CBX() f0->cbx += px*( f0->div_b_err - fx->div_b_err )
#define MARDER_CBY() f0->cby += py*( f0->div_b_err - fy->div_b_err )
#define MARDER_CBZ() f0->cbz += pz*( f0->div_b_err - fz->div_b_err )

typedef struct pipeline_args {
  field_t      * ALIGNED(128) f;
  const grid_t *              g;
} pipeline_args_t;

static void
pipeline( pipeline_args_t * args,
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
  
  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

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

void
clean_div_b( field_t      * ALIGNED(128) f,
             const grid_t *              g ) {
  pipeline_args_t args[1];
  
  float alphadt, px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field"));
  if( g==NULL ) ERROR(("Bad grid"));

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
  EXEC_PIPELINES( pipeline, args, 0 );
  
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
