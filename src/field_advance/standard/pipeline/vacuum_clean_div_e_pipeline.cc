#define IN_sfa
#define IN_vacuum_clean_div_e_pipeline

#include "vacuum_clean_div_e_pipeline.h"

#include "../sfa_private.h"

#include "../../../util/pipelines/pipelines_exec.h"

static void
vacuum_clean_div_e_pipeline_scalar( pipeline_args_t * args,
                                    int pipeline_rank,
                                    int n_pipeline )
{
  DECLARE_STENCIL();
  
  int n_voxel;

  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  INIT_STENCIL();
  for( ; n_voxel; n_voxel-- )
  {
    MARDER_EX();
    MARDER_EY();
    MARDER_EZ();

    NEXT_STENCIL();
  }
}

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "Not implemented"

#endif

void
vacuum_clean_div_e_pipeline( field_array_t * fa )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Do majority of field components in single pass on the pipelines.
  // The host handles stragglers.

  pipeline_args_t args[1];

  args->f = fa->f;
  args->p = (sfa_params_t *) fa->params;
  args->g = fa->g;

  EXEC_PIPELINES( vacuum_clean_div_e, args, 0 );

  // While pipelines are busy, do left overs on the host

  DECLARE_STENCIL();

  // Do left over ex
  for( y = 1; y <= ny+1; y++ )
  {
    f0 = &f( 1, y, nz+1 );
    fx = &f( 2, y, nz+1 );

    for( x = 1; x <= nx; x++ )
    {
      MARDER_EX();

      f0++;
      fx++;
    }
  }

  for( z = 1; z <= nz; z++ )
  {
    f0 = &f( 1, ny+1, z );
    fx = &f( 2, ny+1, z );

    for( x = 1; x <= nx; x++ )
    {
      MARDER_EX();

      f0++;
      fx++;
    }
  }

  // Do left over ey
  for( z = 1; z <= nz+1; z++ )
  {
    for( y = 1; y <= ny; y++ )
    {
      f0 = &f( nx+1, y,   z );
      fy = &f( nx+1, y+1, z );

      MARDER_EY();
    }
  }

  for( y = 1; y <= ny; y++ )
  {
    f0 = &f( 1, y,   nz+1 );
    fy = &f( 1, y+1, nz+1 );

    for( x = 1; x <= nx; x++ )
    {
      MARDER_EY();

      f0++;
      fy++;
    }
  }

  // Do left over ez
  for( z = 1; z <= nz; z++ )
  {
    f0 = &f( 1, ny+1, z   );
    fz = &f( 1, ny+1, z+1 );

    for( x = 1; x <= nx+1; x++ )
    {
      MARDER_EZ();

      f0++;
      fz++;
    }
  }

  for( z = 1; z <= nz; z++ )
  {
    for( y = 1; y <= ny; y++ )
    {
      f0 = &f( nx+1, y, z   );
      fz = &f( nx+1, y, z+1 );

      MARDER_EZ();
    }
  }

  WAIT_PIPELINES();

  local_adjust_tang_e( fa->f, fa->g );
}
