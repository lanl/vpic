#include <field_pipelines.h>

void energy_f( double * global,
               const field_t * ALIGNED f,
               const material_coefficient_t * ALIGNED m,
               const grid_t * ALIGNED g ) {
  energy_f_pipeline_args_t args[1];
  pipeline_request_t request[1];
  double v0;
  int p;
  
  if( global==NULL ) { ERROR(("Bad energy"));                return; }
  if( f==NULL )      { ERROR(("Bad field"));                 return; }
  if( m==NULL )      { ERROR(("Bad material coefficients")); return; }
  if( g==NULL )      { ERROR(("Bad grid"));                  return; }

  /* Have each pipelines work on a portion of the local voxels */
  
  args->f = f;
  args->m = m;
  args->g = g;
  dispatch_pipelines( energy_f_pipeline, args, 0,request );
  wait_for_pipelines( request );

  /* Reduce results from each pipelines */
  
  for( p=1; p<n_pipeline; p++ ) {
    args->en[0][0] += args->en[p][0]; args->en[0][1] += args->en[p][1];
    args->en[0][2] += args->en[p][2]; args->en[0][3] += args->en[p][3];
    args->en[0][4] += args->en[p][4]; args->en[0][5] += args->en[p][5];
  }
    
  /* Convert to physical units and reduce results between nodes */
  
  v0 = 0.5*g->eps0*g->dx*g->dy*g->dz;
  args->en[0][0] *= v0; args->en[0][1] *= v0;
  args->en[0][2] *= v0; args->en[0][3] *= v0;
  args->en[0][4] *= v0; args->en[0][5] *= v0;

  /* Reduce results between nodes */

  mp_allsum_d( args->en[0], global, 6, g->mp );
}

