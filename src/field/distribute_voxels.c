#include <field.h>

int
distribute_voxels( int x0, int x1,     /* range of x-indices (inclusive) */
                   int y0, int y1,     /* range of y-indices (inclusive) */
                   int z0, int z1,     /* range of z-indices (inclusive) */
                   int job, int n_job, /* job ... on [0,n_job-1] */
                   int * _x, int * _y, int * _z ) {
  const int nvx = x1-x0+1;
  const int nvy = y1-y0+1;
  const int nvz = z1-z0+1;

  int x, y, z, n_voxel;

  n_voxel = nvx*nvy*nvz;

  if( job==n_job ) {

    /* The host processes the final incomplete bundle */

    x  = n_voxel;
    n_voxel &= 3;
    x -= n_voxel;

  } else {

    /* Pipelines process a roughly equal share of complete voxel bundles */

    double n_target = (double)( n_voxel>>2 ) / (double)n_job;
    x               = 4*(int)( n_target*(double)( job     ) + 0.5 );
    n_voxel         = 4*(int)( n_target*(double)( job + 1 ) + 0.5 ) - x;

  }

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

