#define IN_vfa
#include "vfa_private.h"

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

typedef struct pipeline_args {
  const field_t                * ALIGNED(128) f;
  const grid_t                 *              g;
  double en[MAX_PIPELINE+1][6];
} pipeline_args_t;

static void
vfa_energy_f_pipeline( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {
  const field_t * ALIGNED(128) f = args->f;
  const grid_t  *              g = args->g;
    
  const field_t * ALIGNED(16) f0;
  const field_t * ALIGNED(16) fx,  * ALIGNED(16) fy,  * ALIGNED(16) fz;
  const field_t * ALIGNED(16) fyz, * ALIGNED(16) fzx, * ALIGNED(16) fxy;
  int x, y, z, n_voxel;
  
  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  double en_ex, en_ey, en_ez, en_bx, en_by, en_bz;
  
  // Process voxels assigned to this pipeline
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );
  
# define LOAD_STENCIL()  \
  f0  = &f(x,  y,  z  ); \
  fx  = &f(x+1,y,  z  ); \
  fy  = &f(x,  y+1,z  ); \
  fz  = &f(x,  y,  z+1); \
  fyz = &f(x,  y+1,z+1); \
  fzx = &f(x+1,y,  z+1); \
  fxy = &f(x+1,y+1,z  )
  
  LOAD_STENCIL();

  // FOR VACUUM REGIONS, THIS COULD BE ACCELERATED

  en_ex = en_ey = en_ez = en_bx = en_by = en_bz = 0;
  for( ; n_voxel; n_voxel-- ) {
    // FIXME: CHECK IF THIS IS THE CORRECT LAGRANGIAN DEFINITION

    en_ex += 0.25*( f0->ex * f0->ex +
                    fy->ex * fy->ex +
                    fz->ex * fz->ex + 
                    fyz->ex*fyz->ex );
    en_ey += 0.25*( f0->ey * f0->ey +
                    fz->ey * fz->ey +
                    fx->ey * fx->ey + 
                    fzx->ey*fzx->ey );
    en_ez += 0.25*( f0->ez * f0->ez +
                    fx->ez * fx->ez +
                    fy->ez * fy->ez + 
                    fxy->ez*fxy->ez );

    en_bx += 0.5 *( f0->cbx* f0->cbx +
                    fx->cbx* fx->cbx );
    en_by += 0.5 *( f0->cby* f0->cby +
                    fy->cby* fy->cby );
    en_bz += 0.5 *( f0->cbz* f0->cbz +
                    fz->cbz* fz->cbz );

    f0++; fx++; fy++; fz++; fyz++; fzx++; fxy++;
    
    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      LOAD_STENCIL();
    }
  }

# undef LOAD_STENCIL

  args->en[pipeline_rank][0] = en_ex;
  args->en[pipeline_rank][1] = en_ey;
  args->en[pipeline_rank][2] = en_ez;
  args->en[pipeline_rank][3] = en_bx;
  args->en[pipeline_rank][4] = en_by;
  args->en[pipeline_rank][5] = en_bz;
}

void
vfa_energy_f( double                       *              global,
              const field_t                * ALIGNED(128) f,
              const material_coefficient_t * ALIGNED(128) m,
              const grid_t                 * ALIGNED(128) g ) {
  pipeline_args_t args[1];
  double v0;
  int p;
  
  if( global==NULL ) ERROR(("Bad energy"));
  if( f==NULL )      ERROR(("Bad field"));
  if( g==NULL )      ERROR(("Bad grid"));

# if 0 // Original non-pipelined version
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
         
         en_ex += 0.25*(  f0->ex * f0->ex +
                          fy->ex * fy->ex +
                          fz->ex * fz->ex + 
                         fyz->ex *fyz->ex );

         en_ey += 0.25*(  f0->ey * f0->ey +
                          fz->ey * fz->ey +
                          fx->ey * fx->ey + 
                         fzx->ey *fzx->ey );
                            
         en_ez += 0.25*(  f0->ez * f0->ez +
                          fx->ez * fx->ez +
                          fy->ez * fy->ez + 
                         fxy->ez *fxy->ez );
         
         en_bx += 0.5 *(  f0->cbx* f0->cbx +
                          fx->cbx* fx->cbx );

         en_by += 0.5 *(  f0->cby* f0->cby +
                          fy->cby* fy->cby );

         en_bz += 0.5 *(  f0->cbz* f0->cbz +
                          fz->cbz* fz->cbz );
                            
         f0++; fx++; fy++; fz++; fyz++; fzx++; fxy++;
      }
    }
  }
# endif

  // Have each pipelines work on a portion of the local voxels
  
  args->f = f;
  args->g = g;

  EXEC_PIPELINES( vfa_energy_f, args, 0 );
  WAIT_PIPELINES();

  // Reduce results from each pipelines
  
  for( p=1; p<=N_PIPELINE; p++ ) {
    args->en[0][0] += args->en[p][0]; args->en[0][1] += args->en[p][1];
    args->en[0][2] += args->en[p][2]; args->en[0][3] += args->en[p][3];
    args->en[0][4] += args->en[p][4]; args->en[0][5] += args->en[p][5];
  }
    
  // Convert to physical units and reduce results between nodes
  
  v0 = 0.5*g->eps0*g->dx*g->dy*g->dz;
  args->en[0][0] *= v0; args->en[0][1] *= v0;
  args->en[0][2] *= v0; args->en[0][3] *= v0;
  args->en[0][4] *= v0; args->en[0][5] *= v0;

  // Reduce results between nodes

  mp_allsum_d( args->en[0], global, 6, g->mp );
}
