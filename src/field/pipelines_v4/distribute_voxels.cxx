#include <v4.h>
#ifdef V4_ACCELERATION
using namespace v4;

#include <field_pipelines.h>
   
int
distribute_voxels_v4( int x0, int x1,     /* range of x-indices (inclusive) */
                      int y0, int y1,     /* range of y-indices (inclusive) */
                      int z0, int z1,     /* range of z-indices (inclusive) */
                      int job, int n_job, /* job ... on [0,n_job-1] */
                      int * _x, int * _y, int * _z ) {
  double n_target;
  int x, y, z, nvx, nvy, nvz, n_voxel;
  
  nvx = x1-x0+1;
  nvy = y1-y0+1;
  nvz = z1-z0+1;
  
  n_target = (double)( nvx*nvy*nvz ) / (double)n_job;
  x        = (int)( n_target*(double)( job     ) + 0.5 );
  n_voxel  = (int)( n_target*(double)( job + 1 ) + 0.5 ) - x;
  
  /**/           /* x = (x-x0) + nvx*( (y-y0) + nvy*(z-z0) ) */
  y  = x/nvx;    /* y =                (y-y0) + nvy*(z-z0)   */
  z  = y/nvy;    /* z =                             (z-z0)   */
  x -= y*nvx;    /* x = (x-x0)                               */
  y -= z*nvy;    /* y =                (y-y0)                */
  
  x += x0;
  y += y0;
  z += z0;

  *_x = x;
  *_y = y;
  *_z = z;

  return n_voxel;
}

#endif

