#define IN_sfa
#define IN_advance_b_pipeline

#define HAS_V4_PIPELINE
#define HAS_V8_PIPELINE
#define HAS_V16_PIPELINE

#include "advance_b_pipeline.h"

#include "../sfa_private.h"

#include "../../../util/pipelines/pipelines_exec.h"

//----------------------------------------------------------------------------//
// Reference implementation for an advance_b pipeline function which does not
// make use of explicit calls to vector intrinsic functions.
//----------------------------------------------------------------------------//

void
advance_b_pipeline_scalar( pipeline_args_t * args,
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
    UPDATE_CBX();
    UPDATE_CBY();
    UPDATE_CBZ();

    NEXT_STENCIL();
  }

# undef LOAD_STENCIL
}

//----------------------------------------------------------------------------//
// Top level function to select and call the proper advance_b pipeline
// function.
//----------------------------------------------------------------------------//

void
advance_b_pipeline( field_array_t * RESTRICT fa,
                    float _frac )
{
  if ( !fa )
  {
    ERROR( ( "Bad args" ) );
  }
  
  // Do the bulk of the magnetic fields in the pipelines.  The host
  // handles stragglers.

  pipeline_args_t args[1];

  args->f    = fa->f;
  args->g    = fa->g;
  args->frac = _frac;

  EXEC_PIPELINES( advance_b, args, 0 );

  // While the pipelines are busy, do surface fields

  DECLARE_STENCIL();

  // Do left over bx
  for( z = 1; z <= nz; z++ )
  {
    for( y = 1; y <= ny; y++ )
    {
      f0 = &f( nx+1, y,   z   );
      fy = &f( nx+1, y+1, z   );
      fz = &f( nx+1, y,   z+1 );

      UPDATE_CBX();
    }
  }

  // Do left over by
  for( z = 1; z <= nz; z++ )
  {
    f0 = &f( 1, ny+1, z   );
    fx = &f( 2, ny+1, z   );
    fz = &f( 1, ny+1, z+1 );

    for( x = 1; x <= nx; x++ )
    {
      UPDATE_CBY();

      f0++;
      fx++;
      fz++;
    }
  }

  // Do left over bz
  for( y = 1; y <= ny; y++ )
  {
    f0 = &f( 1, y,   nz+1 );
    fx = &f( 2, y,   nz+1 );
    fy = &f( 1, y+1, nz+1 );

    for( x = 1; x <= nx; x++ )
    {
      UPDATE_CBZ();

      f0++;
      fx++;
      fy++;
    }
  }

  WAIT_PIPELINES();

  local_adjust_norm_b( f, g );
}
