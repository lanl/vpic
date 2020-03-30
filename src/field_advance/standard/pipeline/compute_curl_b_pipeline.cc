#define IN_sfa
#define IN_compute_curl_b_pipeline

#define HAS_V4_PIPELINE
#define HAS_V8_PIPELINE
#define HAS_V16_PIPELINE

#include "compute_curl_b_pipeline.h"

#include "../sfa_private.h"

#include "../../../util/pipelines/pipelines_exec.h"

//----------------------------------------------------------------------------//
// Reference implementation for a compute_curl_b pipeline function which does
// not make use of explicit calls to vector intrinsic functions.
//----------------------------------------------------------------------------//

void
compute_curl_b_pipeline_scalar( pipeline_args_t * args,
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
// Top level function to select and call the proper compute_curl_b pipeline
// function.
//----------------------------------------------------------------------------//

void
compute_curl_b_pipeline( field_array_t * RESTRICT fa )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  //--------------------------------------------------------------------------//
  // Begin tangential B ghost setup
  //--------------------------------------------------------------------------//

  begin_remote_ghost_tang_b( fa->f, fa->g );

  local_ghost_tang_b( fa->f, fa->g );

  //--------------------------------------------------------------------------//
  // Update interior fields
  //--------------------------------------------------------------------------//
  // Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
  // Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
  // Note: ez all (1:nx+1,1:ny+1,1:nz  ) interior (1:nx,1:ny,2:nz)
  //--------------------------------------------------------------------------//

  // Do majority interior in a single pass.  The host handles
  // stragglers.

  pipeline_args_t args[1];

  args->f = fa->f;
  args->p = (sfa_params_t *)fa->params;
  args->g = fa->g;

  EXEC_PIPELINES( compute_curl_b, args, 0 );

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

  //--------------------------------------------------------------------------//
  // Finish tangential B ghost setup
  //--------------------------------------------------------------------------//

  end_remote_ghost_tang_b( fa->f, fa->g );

  //--------------------------------------------------------------------------//
  // Update exterior fields
  //--------------------------------------------------------------------------//

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
