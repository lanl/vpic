// Note: This is virtually identical to vacuum_compute_rhob

#define IN_sfa
#define IN_vacuum_compute_div_e_err_pipeline

#include "vacuum_compute_div_e_err_pipeline.h"

#include "../sfa_private.h"

#include "../../../util/pipelines/pipelines_exec.h"

void
vacuum_compute_div_e_err_pipeline_scalar( pipeline_args_t * args,
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
    UPDATE_DERR_E();

    NEXT_STENCIL();
  }
}

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "Not implemented"

#endif

void
vacuum_compute_div_e_err_pipeline( field_array_t * RESTRICT fa )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Have pipelines compute the interior of local domain (the host
  // handles stragglers in the interior)

  // Begin setting normal e ghosts

  begin_remote_ghost_norm_e( fa->f, fa->g );

  local_ghost_norm_e( fa->f, fa->g );

  // Have pipelines compute interior of local domain

  pipeline_args_t args[1];

  args->f = fa->f;
  args->p = (sfa_params_t *) fa->params;
  args->g = fa->g;

  EXEC_PIPELINES( vacuum_compute_div_e_err, args, 0 );

  // While pipelines are busy, have host compute the exterior
  // of the local domain

  DECLARE_STENCIL();

  // Finish setting normal e ghosts
  end_remote_ghost_norm_e( fa->f, fa->g );

  // z faces, x edges, y edges and all corners
  for( y = 1; y <= ny+1; y++ )
  {
    f0 = &f( 1, y,   1 );
    fx = &f( 0, y,   1 );
    fy = &f( 1, y-1, 1 );
    fz = &f( 1, y,   0 );

    for( x = 1; x <= nx+1; x++ )
    {
      UPDATE_DERR_E();

      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  for( y = 1; y <= ny+1; y++ )
  {
    f0 = &f( 1, y,   nz+1 );
    fx = &f( 0, y,   nz+1 );
    fy = &f( 1, y-1, nz+1 );
    fz = &f( 1, y,   nz   );

    for( x = 1; x <= nx+1; x++ )
    {
      UPDATE_DERR_E();

      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  // y faces, z edges
  for( z = 2; z <= nz; z++ )
  {
    f0 = &f( 1, 1, z   );
    fx = &f( 0, 1, z   );
    fy = &f( 1, 0, z   );
    fz = &f( 1, 1, z-1 );

    for( x = 1; x <= nx+1; x++ )
    {
      UPDATE_DERR_E();

      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  for( z = 2; z <= nz; z++ )
  {
    f0 = &f( 1, ny+1, z   );
    fx = &f( 0, ny+1, z   );
    fy = &f( 1, ny,   z   );
    fz = &f( 1, ny+1, z-1 );

    for( x = 1; x <= nx+1; x++ )
    {
      UPDATE_DERR_E();

      f0++;
      fx++;
      fy++;
      fz++;
    }
  }

  // x faces
  for( z = 2; z <= nz; z++ )
  {
    for( y = 2; y <= ny; y++ )
    {
      f0 = &f( 1, y,   z   );
      fx = &f( 0, y,   z   );
      fy = &f( 1, y-1, z   );
      fz = &f( 1, y,   z-1 );

      UPDATE_DERR_E();

      f0 = &f( nx+1, y,   z   );
      fx = &f( nx,   y,   z   );
      fy = &f( nx+1, y-1, z   );
      fz = &f( nx+1, y,   z-1 );

      UPDATE_DERR_E();
    }
  }

  // Finish up setting interior

  WAIT_PIPELINES();

  local_adjust_div_e( fa->f, fa->g );
}
