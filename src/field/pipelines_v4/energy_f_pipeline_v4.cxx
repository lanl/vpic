#include <field_pipelines.h>

#ifdef V4_ACCELERATION
using namespace v4;

#if 0 // Original non-pipelined non-vectorized version
  en_ex = ey_ey = en_ez = en_bx = en_by = en_bz = 0;
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0  = &f(1,y,  z  );
      fx  = &f(2,y,  z  );
      fy  = &f(1,y+1,z  );
      fz  = &f(1,y,  z+1);
      fyz = &f(1,y+1,z+1);
      fzx = &f(2,y,  z+1);
      fxy = &f(2,y+1,z  );
      for( x=1; x<=nx; x++ ) {
      
         // FIXME: CHECK IF THIS IS THE CORRECT LAGRANGIAN DEFINITION 
         
         en_ex += 0.25*( m[ f0->ematx].epsx* f0->ex * f0->ex +
                         m[ fy->ematx].epsx* fy->ex * fy->ex +
                         m[ fz->ematx].epsx* fz->ex * fz->ex + 
                         m[fyz->ematx].epsx*fyz->ex *fyz->ex );
         en_ey += 0.25*( m[ f0->ematy].epsy* f0->ey * f0->ey +
                         m[ fz->ematy].epsy* fz->ey * fz->ey +
                         m[ fx->ematy].epsy* fx->ey * fx->ey + 
                         m[fzx->ematy].epsy*fzx->ey *fzx->ey );
         en_ez += 0.25*( m[ f0->ematz].epsz* f0->ez * f0->ez +
                         m[ fx->ematz].epsz* fx->ez * fx->ez +
                         m[ fy->ematz].epsz* fy->ez * fy->ez + 
                         m[fxy->ematz].epsz*fxy->ez *fxy->ez );
         
         en_bx += 0.5 *( m[ f0->fmatx].rmux* f0->cbx* f0->cbx +
                         m[ fx->fmatx].rmux* fx->cbx* fx->cbx );
         en_by += 0.5 *( m[ f0->fmaty].rmuy* f0->cby* f0->cby +
                         m[ fy->fmaty].rmuy* fy->cby* fy->cby );
         en_bz += 0.5 *( m[ f0->fmatz].rmuz* f0->cbz* f0->cbz +
                         m[ fz->fmatz].rmuz* fz->cbz* fz->cbz );
                            
         f0++; fx++; fy++; fz++; fyz++; fzx++; fxy++;

      }
    }
  }
#endif

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
energy_f_pipeline_v4( energy_f_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline ) {
  const field_t                * ALIGNED f = args->f;
  const material_coefficient_t * ALIGNED m = args->m;
  const grid_t                 * ALIGNED g = args->g;

  int x, y, z, n_voxel;
  
  const field_t *f0, *fx, *fy, *fz, *fyz, *fzx, *fxy;
  double en_ex, en_ey, en_ez, en_bx, en_by, en_bz;
  
  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  // Process voxels assigned to this pipeline 
  
  n_voxel = distribute_voxels_v4( 1,nx, 1,ny, 1,nz,
                                  pipeline_rank, n_pipeline,
                                  &x, &y, &z );

# define LOAD_PTRS()                                                   \
  f0  = &f(x,  y,  z  );                                               \
  fx  = &f(x+1,y,  z  ); fy  = &f(x,  y+1,z  ); fz  = &f(x,  y,  z+1); \
  fyz = &f(x,  y+1,z+1); fzx = &f(x+1,y,  z+1); fxy = &f(x+1,y+1,z  )

# define ADVANCE_PTRS()                        \
  f0++; fx++; fy++; fz++; fyz++; fzx++; fxy++; \
  x++;                                         \
  if( x>nx ) {                                 \
    x=1, y++;                                  \
    if( y>ny ) y=1, z++;                       \
    LOAD_PTRS();                               \
  }

  LOAD_PTRS();

  en_ex = en_ey = en_ez = en_bx = en_by = en_bz = 0;
  for( ; n_voxel; n_voxel-- ) {

    // FIXME: CHECK IF THIS IS THE CORRECT LAGRANGIAN DEFINITION 

    // FIXME: I stared at this for a long time and, for a great
    // many reasons, this is not worth vectorizing on the SPUs.
    // Namely, it is not called often enough, it uses V4 unfriendly
    // mixed precision arithmetic, horizontal implementations are
    // almost entirely vector gather operations, efficient vertical
    // implementations will have different round-off behavior from
    // its scalar counterpart, etc.

    en_ex += 0.25*( m[ f0->ematx].epsx* f0->ex * f0->ex +
                    m[ fy->ematx].epsx* fy->ex * fy->ex +
                    m[ fz->ematx].epsx* fz->ex * fz->ex + 
                    m[fyz->ematx].epsx*fyz->ex *fyz->ex );
    en_ey += 0.25*( m[ f0->ematy].epsy* f0->ey * f0->ey +
                    m[ fz->ematy].epsy* fz->ey * fz->ey +
                    m[ fx->ematy].epsy* fx->ey * fx->ey + 
                    m[fzx->ematy].epsy*fzx->ey *fzx->ey );                            
    en_ez += 0.25*( m[ f0->ematz].epsz* f0->ez * f0->ez +
                    m[ fx->ematz].epsz* fx->ez * fx->ez +
                    m[ fy->ematz].epsz* fy->ez * fy->ez + 
                    m[fxy->ematz].epsz*fxy->ez *fxy->ez );

    en_bx += 0.5 *( m[ f0->fmatx].rmux* f0->cbx* f0->cbx +
                    m[ fx->fmatx].rmux* fx->cbx* fx->cbx );
    en_by += 0.5 *( m[ f0->fmaty].rmuy* f0->cby* f0->cby +
                    m[ fy->fmaty].rmuy* fy->cby* fy->cby );
    en_bz += 0.5 *( m[ f0->fmatz].rmuz* f0->cbz* f0->cbz +
                    m[ fz->fmatz].rmuz* fz->cbz* fz->cbz );

    ADVANCE_PTRS();

  }

  args->en[pipeline_rank][0] = en_ex;
  args->en[pipeline_rank][1] = en_ey;
  args->en[pipeline_rank][2] = en_ez;
  args->en[pipeline_rank][3] = en_bx;
  args->en[pipeline_rank][4] = en_by;
  args->en[pipeline_rank][5] = en_bz;

}

#endif
