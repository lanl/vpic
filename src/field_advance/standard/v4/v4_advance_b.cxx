#define IN_sfa_v4
#define HAS_V4_PIPELINE
#include "sfa_v4_private.h"

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]
 
// WTF!  Under -ffast-math, gcc-4.1.1 thinks it is okay to treat the
// below as
//   f0->cbx = ( f0->cbx + py*( blah ) ) - pz*( blah )
// even with explicit parenthesis are in there!  Oh my ...
// -fno-unsafe-math-optimizations must be used

#define UPDATE_CBX() f0->cbx -= ( py*( fy->ez-f0->ez ) - pz*( fz->ey-f0->ey ) )
#define UPDATE_CBY() f0->cby -= ( pz*( fz->ex-f0->ex ) - px*( fx->ez-f0->ez ) )
#define UPDATE_CBZ() f0->cbz -= ( px*( fx->ey-f0->ey ) - py*( fy->ex-f0->ex ) )

// FIXME: MERGE WITH ADVANCE_B IN PREVIOUS DIRECTORY TO ELIMINATE HIDEOUS
// HACKS LIKE THIS!

typedef struct pipeline_args {
  field_t      * ALIGNED(128) f;
  const grid_t *              g;
  float frac;
} pipeline_args_t;

extern "C" void
advance_b_pipeline( pipeline_args_t * args,
                    int pipeline_rank,
                    int n_pipeline );

static void
advance_b_pipeline_v4( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {

  using namespace v4;

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

  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);

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
  
  // Process the voxels assigned to this pipeline 
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );
  
  // Process the bulk of the voxels 4 at a time

# define LOAD_STENCIL()  \
  f0 = &f(x,  y,  z  );  \
  fx = &f(x+1,y,  z  );  \
  fy = &f(x,  y+1,z  );  \
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

    load_4x3_tr(  &f00->ex,  &f01->ex,  &f02->ex,  &f03->ex,  f0_ex,  f0_ey,  f0_ez  );
    load_4x3_tr(  &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx, f0_cbx, f0_cby, f0_cbz );

    load_4x3_tr(  &fx0->ex,  &fx1->ex,  &fx2->ex,  &fx3->ex,  dummy,  fx_ey,  fx_ez  );
    load_4x3_tr(  &fy0->ex,  &fy1->ex,  &fy2->ex,  &fy3->ex,  fy_ex,  dummy,  fy_ez  );
    load_4x2_tr(  &fz0->ex,  &fz1->ex,  &fz2->ex,  &fz3->ex,  fz_ex,  fz_ey   /**/   );

    f0_cbx += fnms( vpy,( fy_ez-f0_ez ),  vpz*( fz_ey-f0_ey ) );
    f0_cby += fnms( vpz,( fz_ex-f0_ex ),  vpx*( fx_ez-f0_ez ) );
    f0_cbz += fnms( vpx,( fx_ey-f0_ey ),  vpy*( fy_ex-f0_ex ) );

    store_4x3_tr( f0_cbx, f0_cby, f0_cbz, &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL

}

void
v4_advance_b( field_t      * ALIGNED(128) f,
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

  EXEC_PIPELINES( advance_b, args, 0 );
  
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
