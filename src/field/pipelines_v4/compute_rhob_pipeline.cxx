// Note: This is virtually identical to compute_div_e_err_pipeline

#if 0 // Original non-pipelined non-vectorized version 
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {
        f0->rhob = m[f0->nmat].nonconductive *
          ( px*( m[f0->ematx].epsx*f0->ex - m[fx->ematx].epsx*fx->ex ) + 
            py*( m[f0->ematy].epsy*f0->ey - m[fy->ematy].epsy*fy->ey ) + 
            pz*( m[f0->ematz].epsz*f0->ez - m[fz->ematz].epsz*fz->ez ) - 
            f0->rhof );
	f0++; fx++; fy++; fz++;
      }
    }
  }
#endif

#include <field_pipelines.h>
#include <v4.h>

using namespace v4;

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
compute_rhob_pipeline( compute_rhob_pipeline_args_t * args,
                       int pipeline_rank ) {
  field_t                      * ALIGNED f = args->f;
  const material_coefficient_t * ALIGNED m = args->m;
  const grid_t                 *         g = args->g;

  field_t *f0, *fx, *fy, *fz;
  int x, y, z, n_voxel;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float px = (nx>1) ? g->eps0/g->dx : 0;
  const float py = (ny>1) ? g->eps0/g->dy : 0;
  const float pz = (nz>1) ? g->eps0/g->dz : 0;

  // Process voxels assigned to this pipeline 

  n_voxel = distribute_voxels( 2,nx, 2,ny, 2,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

# define LOAD_PTRS()    \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x-1,y,  z  ); \
  fy = &f(x,  y-1,z  ); \
  fz = &f(x,  y,  z-1)

# define ADVANCE_PTRS()   \
  f0++; fx++; fy++; fz++; \
  x++;                    \
  if( x>nx ) {            \
    x=2, y++;             \
    if( y>ny ) y=2, z++;  \
    LOAD_PTRS();          \
  }

  LOAD_PTRS();

  for( ; n_voxel; n_voxel-- ) {

    // FIXME: I stared at this for a long time and came to the
    // conclusion that this is not worth vectorizing even on the
    // SPUs.  Namely, this operation is not executed very often
    // and, in either horizontal or vertical SIMD, would require
    // predominantly scalar gather / scatter operations to
    // assemble the desired vectors.
   
    f0->rhob = m[f0->nmat].nonconductive *
      ( px*( m[f0->ematx].epsx*f0->ex - m[fx->ematx].epsx*fx->ex ) + 
        py*( m[f0->ematy].epsy*f0->ey - m[fy->ematy].epsy*fy->ey ) + 
        pz*( m[f0->ematz].epsz*f0->ez - m[fz->ematz].epsz*fz->ez ) - 
        f0->rhof );

    ADVANCE_PTRS();

  }

}

