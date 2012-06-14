// FIXME: This function assumes that the accumlator ghost values are
// zero.  Further, assumes that the ghost values of jfx, jfy, jfz are
// meaningless.  This might be changed to a more robust but slightly
// slower implementation in the near future.

#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include <sf_interface_private.h>

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]
#define a(x,y,z) a[ VOXEL(x,y,z, nx,ny,nz) ]

void
unload_accumulator_pipeline( unload_accumulator_pipeline_args_t * args,
			     int pipeline_rank,
                             int n_pipeline ) {
  field_t             * ALIGNED(128) f = args->f;
  const accumulator_t * ALIGNED(128) a = args->a;
  
  const accumulator_t * ALIGNED(16) a0;
  const accumulator_t * ALIGNED(16) ax,  * ALIGNED(16) ay,  * ALIGNED(16) az;
  const accumulator_t * ALIGNED(16) ayz, * ALIGNED(16) azx, * ALIGNED(16) axy;
  field_t * ALIGNED(16) f0;
  int x, y, z, n_voxel;
  
  const int nx = args->nx;
  const int ny = args->ny;
  const int nz = args->nz;

  const float cx = args->cx;
  const float cy = args->cy;
  const float cz = args->cz;

  // Process the voxels assigned to this pipeline
  
  if( pipeline_rank==n_pipeline ) return; // No need for straggler cleanup
  DISTRIBUTE_VOXELS( 1,nx+1, 1,ny+1, 1,nz+1, 1,
                     pipeline_rank, n_pipeline, x, y, z, n_voxel );

# define LOAD_STENCIL()                                                 \
  f0  = &f(x,  y,  z  );                                                \
  a0  = &a(x,  y,  z  );                                                \
  ax  = &a(x-1,y,  z  ); ay  = &a(x,  y-1,z  ); az  = &a(x,  y,  z-1);  \
  ayz = &a(x,  y-1,z-1); azx = &a(x-1,y,  z-1); axy = &a(x-1,y-1,z  )

  LOAD_STENCIL();

  for( ; n_voxel; n_voxel-- ) {

    f0->jfx += cx*( a0->jx[0] + ay->jx[1] + az->jx[2] + ayz->jx[3] );
    f0->jfy += cy*( a0->jy[0] + az->jy[1] + ax->jy[2] + azx->jy[3] );
    f0->jfz += cz*( a0->jz[0] + ax->jz[1] + ay->jz[2] + axy->jz[3] );

    f0++; a0++; ax++; ay++; az++; ayz++; azx++; axy++;

    x++;
    if( x>nx+1 ) {
      x=1, y++;
      if( y>ny+1 ) y=1, z++;
      LOAD_STENCIL();
    }

  }

# undef LOAD_STENCIL

}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

// SPU pipeline defined in a different compile unit

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "V4 version not hooked up yet!"

#endif

void
unload_accumulator_array( /**/  field_array_t       * RESTRICT fa,
                          const accumulator_array_t * RESTRICT aa ) {
  unload_accumulator_pipeline_args_t args[1];

  if( !fa || !aa || fa->g!=aa->g ) ERROR(( "Bad args" ));

# if 0 // Original non-pipelined version

  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny+1; y++ ) {

      x   = 1;
      f0  = &f(x,  y,  z  );
      a0  = &a(x,  y,  z  );
      ax  = &a(x-1,y,  z  ); ay  = &a(x,  y-1,z  ); az  = &a(x,  y,  z-1);
      ayz = &a(x,  y-1,z-1); azx = &a(x-1,y,  z-1); axy = &a(x-1,y-1,z  );

      for( x=1; x<=nx+1; x++ ) {

        f0->jfx += cx*( a0->jx[0] + ay->jx[1] + az->jx[2] + ayz->jx[3] );
        f0->jfy += cy*( a0->jy[0] + az->jy[1] + ax->jy[2] + azx->jy[3] );
        f0->jfz += cz*( a0->jz[0] + ax->jz[1] + ay->jz[2] + axy->jz[3] );

        f0++; a0++; ax++; ay++; az++; ayz++; azx++; axy++;

      }
    }
  }

# endif

  args->f  = fa->f;
  args->a  = aa->a;
  args->nx = fa->g->nx;
  args->ny = fa->g->ny;
  args->nz = fa->g->nz;
  args->cx = 0.25*fa->g->rdy*fa->g->rdz/fa->g->dt;
  args->cy = 0.25*fa->g->rdz*fa->g->rdx/fa->g->dt;
  args->cz = 0.25*fa->g->rdx*fa->g->rdy/fa->g->dt;

  EXEC_PIPELINES( unload_accumulator, args, 0 );
  WAIT_PIPELINES();
}
