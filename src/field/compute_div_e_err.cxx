// Note: This is virtually identical to compute_rhob
#define IN_field_pipeline
#include <field_pipelines.h>

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define UPDATE_DERR_E() f0->div_e_err = m[f0->nmat].nonconductive * \
  ( px*( m[f0->ematx].epsx*f0->ex - m[fx->ematx].epsx*fx->ex ) +    \
    py*( m[f0->ematy].epsy*f0->ey - m[fy->ematy].epsy*fy->ey ) +    \
    pz*( m[f0->ematz].epsz*f0->ez - m[fz->ematz].epsz*fz->ez ) -    \
    cj*( f0->rhof + f0->rhob ) )

static void
compute_div_e_err_pipeline( compute_div_e_err_pipeline_args_t * args,
                            int pipeline_rank,
                            int n_pipeline ) {
  field_t                      * ALIGNED(16) f = args->f;
  const material_coefficient_t * ALIGNED(16) m = args->m;
  const grid_t                 *             g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? 1./g->dx : 0;
  const float py = (ny>1) ? 1./g->dy : 0;
  const float pz = (nz>1) ? 1./g->dz : 0;
  const float cj = 1./g->eps0;

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
    UPDATE_DERR_E();
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

#error "SPU version not hooked up yet!"

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

// FIXME: This is probably not worth vectorizing.  Namely, this
// operation is not executed very often and, in either horizontal or
// vertical SIMD, would require predominantly scalar gather / scatter
// operations to assemble the desired vectors.

#error "V4 version not hooked up yet!"

#endif

void
compute_div_e_err( field_t                      * ALIGNED(16) f,
                   const material_coefficient_t * ALIGNED(16) m,
                   const grid_t                 *             g ) {
  compute_div_e_err_pipeline_args_t args[1];  

  float px, py, pz, cj;
  field_t *f0, *fx, *fy, *fz;
  int x, y, z, nx, ny, nz;

  if( f==NULL ) ERROR(("Bad field"));
  if( m==NULL ) ERROR(("Bad material coefficients"));
  if( g==NULL ) ERROR(("Bad grid"));

  // Have pipelines compute the interior of local domain (the host
  // handles stragglers in the interior)

  // FIXME: CHECK IF THIS CAN BE STARTED THIS EARLY

# if 0 // Original non-pipelined version
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
	UPDATE_DERR_E();
	f0++; fx++; fy++; fz++;
      }
    }
  }
# endif

  args->f = f;
  args->m = m;
  args->g = g;

  EXEC_PIPELINES( compute_div_e_err, args, 0 );

  // Have host compute the exterior of the local domain

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  px = (nx>1) ? 1./g->dx : 0;
  py = (ny>1) ? 1./g->dy : 0;
  pz = (nz>1) ? 1./g->dz : 0;
  cj = 1./g->eps0;

  // Begin setting normal e ghosts

  begin_remote_ghost_norm_e( f, g );

  local_ghost_norm_e( f, g );

  // Finish setting normal e ghosts

  end_remote_ghost_norm_e( f, g );

  // Compute divergence error in exterior

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

  WAIT_PIPELINES(); // FIXME: CHECK EXACTLY HOW LATE THIS CAN BE DONE

  local_adjust_div_e(f,g);
}
