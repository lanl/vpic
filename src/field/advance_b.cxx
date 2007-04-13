#include <field.h>

#ifndef V4_ACCELERATION
#define ADVANCE_B_PIPELINE (pipeline_func_t)advance_b_pipeline
#else
#define ADVANCE_B_PIPELINE (pipeline_func_t)advance_b_pipeline_v4
#endif

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
 
/* WTF!  Under -ffast-math, gcc-4.1.1 thinks it is okay to treat the
   below as
     f0->cbx = ( f0->cbx + py*( blah ) ) - pz*( blah )
   even with explicit parenthesis are in there!  Oh my ...
   -ffast-math must not be used */

#define UPDATE_CBX() f0->cbx -= ( py*( fy->ez-f0->ez ) - pz*( fz->ey-f0->ey ) )
#define UPDATE_CBY() f0->cby -= ( pz*( fz->ex-f0->ex ) - px*( fx->ez-f0->ez ) )
#define UPDATE_CBZ() f0->cbz -= ( px*( fx->ey-f0->ey ) - py*( fy->ex-f0->ex ) )

typedef struct advance_b_pipeline_args {
  field_t      * ALIGNED f;
  const grid_t *         g;
  float frac;
} advance_b_pipeline_args_t;

static void
advance_b_pipeline( advance_b_pipeline_args_t * args,
		    int pipeline_rank,
                    int n_pipeline ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;
  const float frac         = args->frac;

  field_t * ALIGNED f0, * ALIGNED fx, * ALIGNED fy, * ALIGNED fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? frac*g->cvac*g->dt/g->dx : 0;
  const float py = (ny>1) ? frac*g->cvac*g->dt/g->dy : 0;
  const float pz = (nz>1) ? frac*g->cvac*g->dt/g->dz : 0;

  /* Process the voxels assigned to this pipeline */
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
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

#ifdef V4_ACCELERATION

using namespace v4;
 
static void
advance_b_pipeline_v4( advance_b_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {
  field_t      * ALIGNED f = args->f;
  const grid_t *         g = args->g;
  const float frac         = args->frac;

  field_t * ALIGNED f0, * ALIGNED fx, * ALIGNED fy, * ALIGNED fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? frac*g->cvac*g->dt/g->dx : 0;
  const float py = (ny>1) ? frac*g->cvac*g->dt/g->dy : 0;
  const float pz = (nz>1) ? frac*g->cvac*g->dt/g->dz : 0;

  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);

  v4float f0_ex,  f0_ey,  f0_ez;  // Voxel quad electric fields
  v4float f0_cbx, f0_cby, f0_cbz; // Voxel quad magnetic fields
  v4float fx_ey, fx_ez;           // Voxel quad +x neighbor fields
  v4float fy_ez, fy_ex;           // Voxel quad +y neighbor fields
  v4float fz_ex, fz_ey;           // Voxel quad +z neighbor fields
  v4float dummy;

  field_t * ALIGNED f00, * ALIGNED f01, * ALIGNED f02, * ALIGNED f03; // Voxel quad
  field_t * ALIGNED fx0, * ALIGNED fx1, * ALIGNED fx2, * ALIGNED fx3; // Voxel quad +x neighbors
  field_t * ALIGNED fy0, * ALIGNED fy1, * ALIGNED fy2, * ALIGNED fy3; // Voxel quad +y neighbors
  field_t * ALIGNED fz0, * ALIGNED fz1, * ALIGNED fz2, * ALIGNED fz3; // Voxel quad +z neighbors
  
  // Process the voxels assigned to this pipeline 
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
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

    load_4x4_tr(  &f00->ex,  &f01->ex,  &f02->ex,  &f03->ex,  f0_ex,  f0_ey,  f0_ez,
                  /**/                                        f0_cbx                 );
    load_4x2_tr(  &f00->cby, &f01->cby, &f02->cby, &f03->cby,         f0_cby, f0_cbz );

    load_4x3_tr(  &fx0->ex,  &fx1->ex,  &fx2->ex,  &fx3->ex,  dummy, fx_ey, fx_ez );
    load_4x3_tr(  &fy0->ex,  &fy1->ex,  &fy2->ex,  &fy3->ex,  fy_ex, dummy, fy_ez );
    load_4x2_tr(  &fz0->ex,  &fz1->ex,  &fz2->ex,  &fz3->ex,  fz_ex, fz_ey  /**/  );

    f0_cbx -= fms( vpy,( fy_ez-f0_ez ),  vpz*( fz_ey-f0_ey ) );
    f0_cby -= fms( vpz,( fz_ex-f0_ex ),  vpx*( fx_ez-f0_ez ) );
    f0_cbz -= fms( vpx,( fx_ey-f0_ey ),  vpy*( fy_ex-f0_ex ) );

    store_4x1_tr( f0_cbx,         &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx );
    store_4x2_tr( f0_cby, f0_cbz, &f00->cby, &f01->cby, &f02->cby, &f03->cby );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL

}

#endif

void
advance_b( field_t * ALIGNED f,
           const grid_t * g,
           float frac ) {
  advance_b_pipeline_args_t args[1];
  
  float px, py, pz;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field")); 
  if( g==NULL ) ERROR(("Bad grid"));

  /* Do the bulk of the magnetic fields in the pipelines.  The host
     handles stragglers. */

# if 0 /* Original non-pipeline version */
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

  PMETHOD.dispatch( ADVANCE_B_PIPELINE, args, 0 );
  advance_b_pipeline( args, PMETHOD.n_pipeline, PMETHOD.n_pipeline );
  
  /* While the pipelines are busy, do surface fields */
  
  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? frac*g->cvac*g->dt/g->dx : 0;
  py = (ny>1) ? frac*g->cvac*g->dt/g->dy : 0;
  pz = (nz>1) ? frac*g->cvac*g->dt/g->dz : 0;

  /* Do left over bx */
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fy = &f(nx+1,y+1,z);
      fz = &f(nx+1,y,  z+1);
      UPDATE_CBX();
    }
  }

  /* Do left over by */
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

  /* Do left over bz */
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
  
  PMETHOD.wait();
}
