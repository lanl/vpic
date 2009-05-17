#include "grid.h"

int
distribute_voxels( int x0, int x1,     // range of x-indices (inclusive)
                   int y0, int y1,     // range of y-indices (inclusive)
                   int z0, int z1,     // range of z-indices (inclusive)
                   int bundle,         // number of voxels in a bundle
                   int job, int n_job, // job ... on [0,n_job-1]
                   int * x, int * y, int * z ) {
  int v, n_voxel;
  DISTRIBUTE_VOXELS( x0,x1, y0,y1, z0,z1, bundle, job,n_job,
                     v,*x,*y,*z,n_voxel );
  return n_voxel;
}
