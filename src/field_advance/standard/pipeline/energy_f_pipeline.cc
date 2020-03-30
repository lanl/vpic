// FIXME: USE THE DISCRETIZED VARIATIONAL PRINCIPLE DEFINITION OF ENERGY

#define IN_sfa
#define IN_energy_f_pipeline

#include "energy_f_pipeline.h"

#include "../sfa_private.h"

#include "../../../util/pipelines/pipelines_exec.h"

void
energy_f_pipeline_scalar( pipeline_args_t * args,
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
    REDUCE_EN();
    NEXT_STENCIL();
  }

  args->en[ pipeline_rank ][ 0 ] = en_ex;
  args->en[ pipeline_rank ][ 1 ] = en_ey;
  args->en[ pipeline_rank ][ 2 ] = en_ez;
  args->en[ pipeline_rank ][ 3 ] = en_bx;
  args->en[ pipeline_rank ][ 4 ] = en_by;
  args->en[ pipeline_rank ][ 5 ] = en_bz;
}

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "Not implemented"

#endif

void
energy_f_pipeline( double * global,
                   const field_array_t * RESTRICT fa )
{
  if ( !global || !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Have each pipeline and the host handle a portion of the
  // local voxels

  pipeline_args_t args[1];

  args->f = fa->f;
  args->p = (sfa_params_t *) fa->params;
  args->g = fa->g;

  EXEC_PIPELINES( energy_f, args, 0 );

  WAIT_PIPELINES();

  // Reduce results from each pipelines

  int p;
  for( p = 1; p <= N_PIPELINE; p++ )
  {
    args->en[ 0 ][ 0 ] += args->en[ p ][ 0 ];
    args->en[ 0 ][ 1 ] += args->en[ p ][ 1 ];
    args->en[ 0 ][ 2 ] += args->en[ p ][ 2 ];
    args->en[ 0 ][ 3 ] += args->en[ p ][ 3 ];
    args->en[ 0 ][ 4 ] += args->en[ p ][ 4 ];
    args->en[ 0 ][ 5 ] += args->en[ p ][ 5 ];
  }

  // Convert to physical units and reduce results between nodes

  double v0 = 0.5 * fa->g->eps0 * fa->g->dV;

  args->en[ 0 ][ 0 ] *= v0;
  args->en[ 0 ][ 1 ] *= v0;
  args->en[ 0 ][ 2 ] *= v0;
  args->en[ 0 ][ 3 ] *= v0;
  args->en[ 0 ][ 4 ] *= v0;
  args->en[ 0 ][ 5 ] *= v0;

  // Reduce results between nodes

  mp_allsum_d( args->en[0], global, 6 );
}
