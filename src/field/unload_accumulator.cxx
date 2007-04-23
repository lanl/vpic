#include <field.h>

#ifndef V4_ACCELERATION
#define UNLOAD_ACCUMULATOR_PIPELINE (pipeline_func_t)unload_accumulator_pipeline
#else
#define UNLOAD_ACCUMULATOR_PIPELINE (pipeline_func_t)unload_accumulator_pipeline_v4
#endif

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
#define a(x,y,z) a[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

typedef struct unload_accumulator_pipeline_args {
  field_t             * ALIGNED(16)  f;
  const accumulator_t * ALIGNED(128) a;
  const grid_t        *              g;
} unload_accumulator_pipeline_args_t;

static void
unload_accumulator_pipeline( unload_accumulator_pipeline_args_t * args,
			     int pipeline_rank,
                             int n_pipeline ) {
  field_t             * ALIGNED(16)  f = args->f;
  const accumulator_t * ALIGNED(128) a = args->a;
  const grid_t        *              g = args->g;
  
  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx,  * ALIGNED(16) fy,  * ALIGNED(16) fz;
  field_t * ALIGNED(16) fyz, * ALIGNED(16) fzx, * ALIGNED(16) fxy;
  const float * ALIGNED(16) pa;
  int x, y, z, n_voxel;
  
  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float cx = 0.25 / (g->dt*g->dy*g->dz);
  const float cy = 0.25 / (g->dt*g->dz*g->dx);
  const float cz = 0.25 / (g->dt*g->dx*g->dy);

  // Process the voxels assigned to this pipeline
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

# define LOAD_STENCIL()        \
  pa  = &a(x,  y,  z  ).jx[0]; \
  f0  = &f(x,  y,  z  );       \
  fx  = &f(x+1,y,  z  );       \
  fy  = &f(x,  y+1,z  );       \
  fz  = &f(x,  y,  z+1);       \
  fyz = &f(x,  y+1,z+1);       \
  fzx = &f(x+1,y,  z+1);       \
  fxy = &f(x+1,y+1,z  )

  LOAD_STENCIL();

  for( ; n_voxel; n_voxel-- ) {
    f0->jfx  += cx*pa[0];  // f(x,y,  z  ).jfx += a(x,y,z).jx[0]
    fy->jfx  += cx*pa[1];  // f(x,y+1,z  ).jfx += a(x,y,z).jx[1]
    fz->jfx  += cx*pa[2];  // f(x,y,  z+1).jfx += a(x,y,z).jx[2]
    fyz->jfx += cx*pa[3];  // f(x,y+1,z+1).jfx += a(x,y,z).jx[3]

    f0->jfy  += cy*pa[4];  // f(x,  y,z  ).jfy += a(x,y,z).jy[0]
    fz->jfy  += cy*pa[5];  // f(x,  y,z+1).jfy += a(x,y,z).jy[1]
    fx->jfy  += cy*pa[6];  // f(x+1,y,z  ).jfy += a(x,y,z).jy[2]
    fzx->jfy += cy*pa[7];  // f(x+1,y,z+1).jfy += a(x,y,z).jy[3]

    f0->jfz  += cz*pa[8];  // f(x,  y,  z).jfz += a(x,y,z).jz[0]
    fx->jfz  += cz*pa[9];  // f(x+1,y,  z).jfz += a(x,y,z).jz[1]
    fy->jfz  += cz*pa[10]; // f(x,  y+1,z).jfz += a(x,y,z).jz[2]
    fxy->jfz += cz*pa[11]; // f(x+1,y+1,z).jfz += a(x,y,z).jz[3]

    pa+=12; f0++; fx++; fy++; fz++; fyz++; fzx++; fxy++;
    
    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      LOAD_STENCIL();
    }
  }

# undef LOAD_STENCIL

}

#ifdef V4_ACCELERATION

using namespace v4;

static void
unload_accumulator_pipeline_v4( unload_accumulator_pipeline_args_t * args,
                                int pipeline_rank,
                                int n_pipeline ) {
  field_t             * ALIGNED(16)  f = args->f;
  const accumulator_t * ALIGNED(128) a = args->a;
  const grid_t        *              g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx,  * ALIGNED(16) fy,  * ALIGNED(16) fz;
  field_t * ALIGNED(16) fyz, * ALIGNED(16) fzx, * ALIGNED(16) fxy;
  const float * ALIGNED(16) pa;
  int x, y, z, n_voxel;
  
  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const v4float cx( 0.25 / (g->dt*g->dy*g->dz) );
  const v4float cy( 0.25 / (g->dt*g->dz*g->dx) );
  const v4float cz( 0.25 / (g->dt*g->dx*g->dy) );

  v4float jfx, jfy, jfz, ajfx, ajfy, ajfz;

  // Process the voxels assigned to this pipeline 
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

# define LOAD_STENCIL()        \
  pa  = &a(x,  y,  z  ).jx[0]; \
  f0  = &f(x,  y,  z  );       \
  fx  = &f(x+1,y,  z  );       \
  fy  = &f(x,  y+1,z  );       \
  fz  = &f(x,  y,  z+1);       \
  fyz = &f(x,  y+1,z+1);       \
  fzx = &f(x+1,y,  z+1);       \
  fxy = &f(x+1,y+1,z  )

  LOAD_STENCIL();

  for( ; n_voxel; n_voxel-- ) {

    load_4x1_tr( &f0->jfx, &fy->jfx, &fz->jfx, &fyz->jfx, jfx ); load_4x1( pa,   ajfx );
    load_4x1_tr( &f0->jfy, &fz->jfy, &fx->jfy, &fzx->jfy, jfy ); load_4x1( pa+4, ajfy );
    load_4x1_tr( &f0->jfz, &fx->jfz, &fy->jfz, &fxy->jfz, jfz ); load_4x1( pa+8, ajfz );

    store_4x1_tr( fma( ajfx, cx, jfx ), &f0->jfx, &fy->jfx, &fz->jfx, &fyz->jfx );
    store_4x1_tr( fma( ajfy, cy, jfy ), &f0->jfy, &fz->jfy, &fx->jfy, &fzx->jfy );
    store_4x1_tr( fma( ajfz, cz, jfz ), &f0->jfz, &fx->jfz, &fy->jfz, &fxy->jfz );

    pa+=12; f0++; fx++; fy++; fz++; fyz++; fzx++; fxy++;
    
    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      LOAD_STENCIL();
    }

  }

# undef LOAD_STENCIL

}

#endif

void
unload_accumulator( field_t             * ALIGNED(16)  f, 
                    const accumulator_t * ALIGNED(128) a,
                    const grid_t        *              g ) {
  unload_accumulator_pipeline_args_t args[1];
  
  if( f==NULL ) ERROR(("Bad field"));
  if( a==NULL ) ERROR(("Bad accumulator"));
  if( g==NULL ) ERROR(("Bad grid"));

# if 0 // Original non-pipelined version
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {

      pa  = &a(1,y,z).jx[0];
      f0  = &f(1,y,  z  );
      fx  = &f(2,y,  z  );
      fy  = &f(1,y+1,z  );
      fz  = &f(1,y,  z+1);
      fyz = &f(1,y+1,z+1);
      fzx = &f(2,y,  z+1);
      fxy = &f(2,y+1,z  );

      for( x=1; x<=nx; x++ ) {

        f0->jfx  += cx*pa[0];  // f(x,y,  z  ).jfx += a(x,y,z).jx[0]
        fy->jfx  += cx*pa[1];  // f(x,y+1,z  ).jfx += a(x,y,z).jx[1]
        fz->jfx  += cx*pa[2];  // f(x,y,  z+1).jfx += a(x,y,z).jx[2]
        fyz->jfx += cx*pa[3];  // f(x,y+1,z+1).jfx += a(x,y,z).jx[3]

        f0->jfy  += cy*pa[4];  // f(x,  y,z  ).jfy += a(x,y,z).jy[0]
        fz->jfy  += cy*pa[5];  // f(x,  y,z+1).jfy += a(x,y,z).jy[1]
        fx->jfy  += cy*pa[6];  // f(x+1,y,z  ).jfy += a(x,y,z).jy[2]
        fzx->jfy += cy*pa[7];  // f(x+1,y,z+1).jfy += a(x,y,z).jy[3]

        f0->jfz  += cz*pa[8];  // f(x,  y,  z).jfz += a(x,y,z).jz[0]
        fx->jfz  += cz*pa[9];  // f(x+1,y,  z).jfz += a(x,y,z).jz[1]
        fy->jfz  += cz*pa[10]; // f(x,  y+1,z).jfz += a(x,y,z).jz[2]
        fxy->jfz += cz*pa[11]; // f(x+1,y+1,z).jfz += a(x,y,z).jz[3]

        pa += 12; f0++; fx++; fy++; fz++; fyz++; fzx++; fxy++;
      }
    }
  }
# endif

  args->f = f;
  args->a = a;
  args->g = g;

  PSTYLE.dispatch( UNLOAD_ACCUMULATOR_PIPELINE, args, 0 );
  unload_accumulator_pipeline( args, PSTYLE.n_pipeline, PSTYLE.n_pipeline );
  PSTYLE.wait();
}
