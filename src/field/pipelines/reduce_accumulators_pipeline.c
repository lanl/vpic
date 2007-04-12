#include <field_pipelines.h>

#define a(x,y,z) a[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
reduce_accumulators_pipeline( reduce_accumulators_pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline ) {
  accumulator_t * ALIGNED a = args->a;
  const grid_t  *         g = args->g;
  int n_voxel;
  
  float * ALIGNED pa, * ALIGNED paa;
  int x, y, z, nx, ny, nz, stride;
  int p;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;
  
  n_voxel = distribute_voxels( 1,nx, 1,ny, 1,nz,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

  stride = 12*(nx+2)*(ny+2)*(nz+2);

  pa = &a(x,y,z).jx[0];
  
  for( ; n_voxel; n_voxel-- ) {
  
    paa = pa + stride;
    for( p=0; p<n_pipeline; p++ ) { 
      pa[0]  += paa[0];
      pa[1]  += paa[1];
      pa[2]  += paa[2];
      pa[3]  += paa[3];
      pa[4]  += paa[4];
      pa[5]  += paa[5];
      pa[6]  += paa[6];
      pa[7]  += paa[7];
      pa[8]  += paa[8];
      pa[9]  += paa[9];
      pa[10] += paa[10];
      pa[11] += paa[11];
      paa    += stride;
    }
    pa += 12;
    
    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
      pa = &a(x,y,z).jx[0];
    }
  }
}
