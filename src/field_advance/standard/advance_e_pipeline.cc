#define IN_sfa

#define HAS_V4_PIPELINE
#define HAS_V8_PIPELINE
#define HAS_V16_PIPELINE

#include "sfa_private.h"

#include "../../util/pipelines/pipelines_exec.h"

typedef struct pipeline_args
{
  field_t            * ALIGNED(128) f;
  const sfa_params_t *              p;
  const grid_t       *              g;
} pipeline_args_t;

#define DECLARE_STENCIL()                                        \
        field_t                * ALIGNED(128) f = args->f;       \
  const material_coefficient_t * ALIGNED(128) m = args->p->mc;   \
  const grid_t                 *              g = args->g;       \
  const int nx = g->nx, ny = g->ny, nz = g->nz;                  \
                                                                 \
  const float damp = args->p->damp;                              \
  const float px   = (nx>1) ? (1+damp)*g->cvac*g->dt*g->rdx : 0; \
  const float py   = (ny>1) ? (1+damp)*g->cvac*g->dt*g->rdy : 0; \
  const float pz   = (nz>1) ? (1+damp)*g->cvac*g->dt*g->rdz : 0; \
  const float cj   = g->dt/g->eps0;                              \
                                                                 \
  field_t * ALIGNED(16) f0;                                      \
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;  \
  int x, y, z

#define f(x,y,z) f[ VOXEL( x, y, z, nx, ny, nz ) ]

#define INIT_STENCIL()  \
  f0 = &f( x,   y,   z   ); \
  fx = &f( x-1, y,   z   ); \
  fy = &f( x,   y-1, z   ); \
  fz = &f( x,   y,   z-1 )

#define NEXT_STENCIL()                      \
  f0++; fx++; fy++; fz++; x++;              \
  if ( x > nx )                             \
  {                                         \
                  y++;               x = 2; \
    if ( y > ny ) z++; if ( y > ny ) y = 2; \
    INIT_STENCIL();                         \
  }

#define UPDATE_EX()                                         \
  f0->tcax = ( py * ( f0->cbz * m[f0->fmatz].rmuz -         \
		      fy->cbz * m[fy->fmatz].rmuz ) -       \
               pz * ( f0->cby * m[f0->fmaty].rmuy -         \
		      fz->cby * m[fz->fmaty].rmuy ) ) -     \
             damp * f0->tcax;                               \
  f0->ex   = m[f0->ematx].decayx * f0->ex +                 \
             m[f0->ematx].drivex * ( f0->tcax - cj * f0->jfx )

#define UPDATE_EY()                                         \
  f0->tcay = ( pz * ( f0->cbx * m[f0->fmatx].rmux -         \
		      fz->cbx * m[fz->fmatx].rmux ) -       \
               px * ( f0->cbz * m[f0->fmatz].rmuz -         \
		      fx->cbz * m[fx->fmatz].rmuz ) ) -     \
             damp * f0->tcay;                               \
  f0->ey   = m[f0->ematy].decayy * f0->ey +                 \
             m[f0->ematy].drivey * ( f0->tcay - cj * f0->jfy )

#define UPDATE_EZ()                                         \
  f0->tcaz = ( px * ( f0->cby * m[f0->fmaty].rmuy -         \
		      fx->cby * m[fx->fmaty].rmuy) -        \
               py * ( f0->cbx * m[f0->fmatx].rmux -         \
		      fy->cbx * m[fy->fmatx].rmux ) ) -     \
             damp * f0->tcaz;                               \
  f0->ez   = m[f0->ematz].decayz * f0->ez +                 \
             m[f0->ematz].drivez * ( f0->tcaz - cj * f0->jfz )

//----------------------------------------------------------------------------//
// Reference implementation for an advance_e pipeline function which does not
// make use of explicit calls to vector intrinsic functions.
//----------------------------------------------------------------------------//

void
advance_e_pipeline_scalar( pipeline_args_t * args,
                           int pipeline_rank,
                           int n_pipeline )
{
  DECLARE_STENCIL();

  int n_voxel;

  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  INIT_STENCIL();

  for( ; n_voxel; n_voxel-- )
  {
    UPDATE_EX();
    UPDATE_EY();
    UPDATE_EZ();

    NEXT_STENCIL();
  }
}

//----------------------------------------------------------------------------//
// If using v4, include an implementation for advance_e_pipeline_v4.
//----------------------------------------------------------------------------//

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#include "advance_e_pipeline_v4.cc"

#endif

//----------------------------------------------------------------------------//
// If using v8, include an implementation for advance_e_pipeline_v8.
//----------------------------------------------------------------------------//

#if defined(V8_ACCELERATION) && defined(HAS_V8_PIPELINE)

#include "advance_e_pipeline_v8.cc"

#endif

//----------------------------------------------------------------------------//
// If using v16, include an implementation for advance_e_pipeline_v16.
//----------------------------------------------------------------------------//

#if defined(V16_ACCELERATION) && defined(HAS_V16_PIPELINE)

#include "advance_e_pipeline_v16.cc"

#endif

//----------------------------------------------------------------------------//
// Top level function to select and call the proper advance_e pipeline
// function.
//----------------------------------------------------------------------------//

void
advance_e_pipeline( field_array_t * RESTRICT fa,
                    float                    frac )
{
  if ( !fa  )
  {
    ERROR( ( "Bad args" ) );
  }

  if ( frac != 1 )
  {
    ERROR( ( "standard advance_e does not support frac!=1 yet" ) );
  }

  /***************************************************************************
   * Begin tangential B ghost setup
   ***************************************************************************/
  
  begin_remote_ghost_tang_b( fa->f, fa->g );

  local_ghost_tang_b( fa->f, fa->g );

  /***************************************************************************
   * Update interior fields
   * Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
   * Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
   * Note: ez all (1:nx+1,1:ny+1,1:nz  ) interior (1:nx,1:ny,2:nz)
   ***************************************************************************/

  // Do majority interior in a single pass.  The host handles
  // stragglers.

  pipeline_args_t args[1];
  args->f = fa->f;
  args->p = (sfa_params_t *)fa->params;
  args->g = fa->g;

  EXEC_PIPELINES( advance_e, args, 0 );
  
  // While the pipelines are busy, do non-bulk interior fields

  DECLARE_STENCIL();

  // Do left over interior ex
  for( z = 2; z <= nz; z++ )
  {
    for( y = 2; y <= ny; y++ )
    {
      f0 = &f( 1, y,   z   );
      fy = &f( 1, y-1, z   );
      fz = &f( 1, y,   z-1 );

      UPDATE_EX();
    }
  }

  // Do left over interior ey
  for( z = 2; z <= nz; z++ )
  {
    f0 = &f( 2, 1, z   );
    fx = &f( 1, 1, z   );
    fz = &f( 2, 1, z-1 );

    for( x = 2; x <= nx; x++ )
    {
      UPDATE_EY();

      f0++;
      fx++;
      fz++;
    }
  }

  // Do left over interior ez
  for( y = 2; y <= ny; y++ )
  {
    f0 = &f( 2, y,   1 );
    fx = &f( 1, y,   1 );
    fy = &f( 2, y-1, 1 );

    for( x = 2; x <= nx; x++ )
    {
      UPDATE_EZ();

      f0++;
      fx++;
      fy++;
    }
  }

  WAIT_PIPELINES();
  
  /***************************************************************************
   * Finish tangential B ghost setup
   ***************************************************************************/

  end_remote_ghost_tang_b( fa->f, fa->g );

  /***************************************************************************
   * Update exterior fields
   ***************************************************************************/

  // Do exterior ex
  for( y = 1; y <= ny+1; y++ )
  {
    f0 = &f( 1, y,   1 );
    fy = &f( 1, y-1, 1 );
    fz = &f( 1, y,   0 );

    for( x = 1; x <= nx; x++ )
    {
      UPDATE_EX();

      f0++;
      fy++;
      fz++;
    }
  }

  for( y = 1; y <= ny+1; y++ )
  {
    f0 = &f( 1, y,   nz+1 );
    fy = &f( 1, y-1, nz+1 );
    fz = &f( 1, y,   nz   );

    for( x = 1; x <= nx; x++ )
    {
      UPDATE_EX();

      f0++;
      fy++;
      fz++;
    }
  }

  for( z = 2; z <= nz; z++ )
  {
    f0 = &f( 1, 1, z   );
    fy = &f( 1, 0, z   );
    fz = &f( 1, 1, z-1 );

    for( x = 1; x <= nx; x++ )
    {
      UPDATE_EX();

      f0++;
      fy++;
      fz++;
    }
  }

  for( z = 2; z <= nz; z++ )
  {
    f0 = &f( 1, ny+1, z   );
    fy = &f( 1, ny,   z   );
    fz = &f( 1, ny+1, z-1 );

    for( x = 1; x <= nx; x++ )
    {
      UPDATE_EX();

      f0++;
      fy++;
      fz++;
    }
  }

  // Do exterior ey
  for( z = 1; z <= nz+1; z++ )
  {
    for( y = 1; y <= ny; y++ )
    {
      f0 = &f( 1, y, z   );
      fx = &f( 0, y, z   );
      fz = &f( 1, y, z-1 );

      UPDATE_EY();
    }
  }

  for( z = 1; z <= nz+1; z++ )
  {
    for( y = 1; y <= ny; y++ )
    {
      f0 = &f( nx+1, y, z   );
      fx = &f( nx,   y, z   );
      fz = &f( nx+1, y, z-1 );

      UPDATE_EY();
    }
  }

  for( y = 1; y <= ny; y++ )
  {
    f0 = &f( 2, y, 1 );
    fx = &f( 1, y, 1 );
    fz = &f( 2, y, 0 );

    for( x = 2; x <= nx; x++ )
    {
      UPDATE_EY();

      f0++;
      fx++;
      fz++;
    }
  }

  for( y = 1; y <= ny; y++ )
  {
    f0 = &f( 2, y, nz+1 );
    fx = &f( 1, y, nz+1 );
    fz = &f( 2, y, nz   );

    for( x = 2; x <= nx; x++ )
    {
      UPDATE_EY();

      f0++;
      fx++;
      fz++;
    }
  }

  // Do exterior ez
  for( z = 1; z <= nz; z++ )
  {
    f0 = &f( 1, 1, z );
    fx = &f( 0, 1, z );
    fy = &f( 1, 0, z );

    for( x = 1; x <= nx+1; x++ )
    {
      UPDATE_EZ();

      f0++;
      fx++;
      fy++;
    }
  }

  for( z = 1; z <= nz; z++ )
  {
    f0 = &f( 1, ny+1, z );
    fx = &f( 0, ny+1, z );
    fy = &f( 1, ny,   z );

    for( x = 1; x <= nx+1; x++ )
    {
      UPDATE_EZ();

      f0++;
      fx++;
      fy++;
    }
  }

  for( z = 1; z <= nz; z++ )
  {
    for( y = 2; y <= ny; y++ )
    {
      f0 = &f( 1, y,   z );
      fx = &f( 0, y,   z );
      fy = &f( 1, y-1, z );

      UPDATE_EZ();
    }
  }

  for( z = 1; z <= nz; z++ )
  {
    for( y = 2; y <= ny; y++ )
    {
      f0 = &f( nx+1, y,   z );
      fx = &f( nx,   y,   z );
      fy = &f( nx+1, y-1, z );

      UPDATE_EZ();
    }
  }

  local_adjust_tang_e( fa->f, fa->g );
}
