// Note: This is virtually identical to compute_div_e_err
#define IN_sfa
#include "sfa_private.h"

typedef struct pipeline_args {
  /**/  field_t      * ALIGNED(128) f;
  const sfa_params_t *              p;
  const grid_t       *              g;
} pipeline_args_t;

#define DECLARE_STENCIL()                                       \
  /**/  field_t                * ALIGNED(128) f = args->f;      \
  const material_coefficient_t * ALIGNED(128) m = args->p->mc;  \
  const grid_t                 *              g = args->g;      \
  const int nx = g->nx, ny = g->ny, nz = g->nz;                 \
                                                                \
  const float px = (nx>1) ? g->eps0*g->rdx : 0;                 \
  const float py = (ny>1) ? g->eps0*g->rdy : 0;                 \
  const float pz = (nz>1) ? g->eps0*g->rdz : 0;                 \
                                                                \
  field_t * ALIGNED(16) f0;                                     \
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz; \
  int x, y, z

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

#define INIT_STENCIL()  \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x-1,y,  z  ); \
  fy = &f(x,  y-1,z  ); \
  fz = &f(x,  y,  z-1)

#define NEXT_STENCIL()                \
  f0++; fx++; fy++; fz++; x++;        \
  if( x>nx ) {                        \
    /**/       y++;            x = 2; \
    if( y>ny ) z++; if( y>ny ) y = 2; \
    INIT_STENCIL();                   \
  }

#define UPDATE_DERR_E() f0->rhob = m[f0->nmat].nonconductive *   \
  ( px*( m[f0->ematx].epsx*f0->ex - m[fx->ematx].epsx*fx->ex ) + \
    py*( m[f0->ematy].epsy*f0->ey - m[fy->ematy].epsy*fy->ey ) + \
    pz*( m[f0->ematz].epsz*f0->ez - m[fz->ematz].epsz*fz->ez ) - \
    f0->rhof )

void
compute_rhob_pipeline( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {
  DECLARE_STENCIL();

  int n_voxel;
  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  INIT_STENCIL();
  for( ; n_voxel; n_voxel-- ) {
    UPDATE_DERR_E();
    NEXT_STENCIL();
  }
}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

// SPU pipeline is defined in a different compilation unit

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "Not implemented"

#endif

void
compute_rhob( field_array_t * RESTRICT fa ) {
  if( !fa ) ERROR(( "Bad args" ));

  // Have pipelines compute the interior of local domain (the host
  // handles stragglers in the interior)

  // Begin setting normal e ghosts

  begin_remote_ghost_norm_e( fa->f, fa->g );
  local_ghost_norm_e( fa->f, fa->g );

  // Have pipelines compute interior of local domain

  pipeline_args_t args[1];  
  args->f = fa->f;
  args->p = (sfa_params_t *)fa->params;
  args->g = fa->g;
  EXEC_PIPELINES( compute_rhob, args, 0 );

  // While pipelines are busy, have host compute the exterior
  // of the local domain

  DECLARE_STENCIL();

  // Finish setting normal e ghosts
  end_remote_ghost_norm_e( fa->f, fa->g );

  // z faces, x edges, y edges and all corners
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  1);
    fx = &f(0,y,  1);
    fy = &f(1,y-1,1);
    fz = &f(1,y,  0);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_DERR_E();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
  for( y=1; y<=ny+1; y++ ) {
    f0 = &f(1,y,  nz+1);
    fx = &f(0,y,  nz+1);
    fy = &f(1,y-1,nz+1);
    fz = &f(1,y,  nz);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_DERR_E();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  // y faces, z edges
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,1,z);
    fx = &f(0,1,z);
    fy = &f(1,0,z);
    fz = &f(1,1,z-1);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_DERR_E();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }
  for( z=2; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(0,ny+1,z);
    fy = &f(1,ny,  z);
    fz = &f(1,ny+1,z-1);
    for( x=1; x<=nx+1; x++ ) {
      UPDATE_DERR_E();
      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  // x faces
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(1,y,  z);
      fx = &f(0,y,  z);
      fy = &f(1,y-1,z);
      fz = &f(1,y,  z-1);
      UPDATE_DERR_E();
      f0 = &f(nx+1,y,  z);
      fx = &f(nx,  y,  z);
      fy = &f(nx+1,y-1,z);
      fz = &f(nx+1,y,  z-1);
      UPDATE_DERR_E();
    }
  }

  // Finish up setting interior

  WAIT_PIPELINES();

  local_adjust_rhob( fa->f, fa->g );
}
