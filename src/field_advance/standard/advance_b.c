#define IN_sfa
#include "sfa_private.h"

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]
 
// WTF!  Under -ffast-math, gcc-4.1.1 thinks it is okay to treat the
// below as
//   f0->cbx = ( f0->cbx + py*( blah ) ) - pz*( blah )
// even with explicit parenthesis are in there!  Oh my ...
// -fno-unsafe-math-optimizations must be used

#define UPDATE_CBX() f0->cbx -= ( py*( fy->ez-f0->ez ) - pz*( fz->ey-f0->ey ) )
#define UPDATE_CBY() f0->cby -= ( pz*( fz->ex-f0->ex ) - px*( fx->ez-f0->ez ) )
#define UPDATE_CBZ() f0->cbz -= ( px*( fx->ey-f0->ey ) - py*( fy->ex-f0->ex ) )

typedef struct pipeline_args {
  field_t      * ALIGNED(128) f;
  const grid_t *              g;
  float frac;
} pipeline_args_t;

static void
pipeline( pipeline_args_t * args,
          int pipeline_rank,
          int n_pipeline ) {
  field_t      * ALIGNED(128) f    = args->f;
  const grid_t *              g    = args->g;
  const float                 frac = args->frac;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? frac*g->cvac*g->dt*g->rdx : 0;
  const float py = (ny>1) ? frac*g->cvac*g->dt*g->rdy : 0;
  const float pz = (nz>1) ? frac*g->cvac*g->dt*g->rdz : 0;

  // Process the voxels assigned to this pipeline
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

# define LOAD_STENCIL() \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x+1,y,  z  ); \
  fy = &f(x,  y+1,z  ); \
  fz = &f(x,  y,  z+1)

  LOAD_STENCIL();

  for( ; n_voxel; n_voxel-- ) {
    UPDATE_CBX();
    UPDATE_CBY();
    UPDATE_CBZ();
    f0++; fx++; fy++; fz++;    

    x++;
    if( x>nx ) {
      x=1,  y++;
      if( y>ny ) y=1, z++;
      LOAD_STENCIL();
    }
  }

# undef LOAD_STENCIL

}

void
advance_b( field_t      * ALIGNED(128) f,
           const grid_t *              g,
           float                       frac ) {
  pipeline_args_t args[1];
  
  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field")); 
  if( g==NULL ) ERROR(("Bad grid"));

  // Do the bulk of the magnetic fields in the pipelines.  The host
  // handles stragglers.

# if 0 // Original non-pipeline version
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(1,y,z);
      fx = &f(2,y,z);
      fy = &f(1,y+1,z);
      fz = &f(1,y,z+1);
      for( x=1; x<=nx; x++ ) {
	UPDATE_CBX();
	UPDATE_CBY();
	UPDATE_CBZ();
	f0++; fx++; fy++; fz++;
      }
    }
  }
# endif
  
  args->f    = f;
  args->g    = g;
  args->frac = frac;

  EXEC_PIPELINES( pipeline, args, 0 );
  
  // While the pipelines are busy, do surface fields
  
  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? frac*g->cvac*g->dt*g->rdx : 0;
  py = (ny>1) ? frac*g->cvac*g->dt*g->rdy : 0;
  pz = (nz>1) ? frac*g->cvac*g->dt*g->rdz : 0;

  // Do left over bx
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fy = &f(nx+1,y+1,z);
      fz = &f(nx+1,y,  z+1);
      UPDATE_CBX();
    }
  }

  // Do left over by
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(2,ny+1,z);
    fz = &f(1,ny+1,z+1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_CBY();
      f0++;
      fx++;
      fz++;
    }
  }

  // Do left over bz
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,  nz+1);
    fx = &f(2,y,  nz+1);
    fy = &f(1,y+1,nz+1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_CBZ();
      f0++;
      fx++;
      fy++;
    }
  }

  local_adjust_norm_b(f,g);
  
  WAIT_PIPELINES(); // FIXME: Check how late this can be done!
}
