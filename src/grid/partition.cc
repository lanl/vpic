/*
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "grid.h"

#define RANK_TO_INDEX(rank,ix,iy,iz) do {               \
    int _ix, _iy, _iz;                                  \
    _ix  = (rank);  /* ix = ix + gpx*( iy + gpy*iz ) */ \
    _iy  = _ix/gpx; /* iy = iy + gpy*iz */              \
    _ix -= _iy*gpx; /* ix = ix */                       \
    _iz  = _iy/gpy; /* iz = iz */                       \
    _iy -= _iz*gpy; /* iy = iy */                       \
    (ix) = _ix;                                         \
    (iy) = _iy;                                         \
    (iz) = _iz;                                         \
  } while(0)

#define INDEX_TO_RANK(ix,iy,iz,rank) do {            \
    int _ix = (ix), _iy = (iy), _iz = (iz);          \
    /* Wrap processor index periodically */          \
    while(_ix>=gpx) _ix-=gpx; while(_ix<0) _ix+=gpx; \
    while(_iy>=gpy) _iy-=gpy; while(_iy<0) _iy+=gpy; \
    while(_iz>=gpz) _iz-=gpz; while(_iz<0) _iz+=gpz; \
    /* Compute the rank */                           \
    (rank) = _ix + gpx*( _iy + gpy*_iz );            \
  } while(0)

void
partition_periodic_box( grid_t * g,
                        double gx0, double gy0, double gz0,
                        double gx1, double gy1, double gz1,
                        int gnx, int gny, int gnz,
                        int gpx, int gpy, int gpz ) {
  double f;
  int rank, px, py, pz;

  // Make sure the grid can be setup

  if( !g ) ERROR(( "NULL grid" ));

  if( gpx<1 || gpy<1 || gpz<1 || gpx*gpy*gpz!=world_size )
    ERROR(( "Bad domain decompostion (%ix%ix%i)", gpx, gpy, gpz ));

  if( gnx<1 || gny<1 || gnz<1 || gnx%gpx!=0 || gny%gpy!=0 || gnz%gpz!=0 )
    ERROR(( "Bad resolution (%ix%ix%i) for domain decomposition",
            gnx, gny, gnz, gpx, gpy, gpz ));

  // Setup basic variables
  RANK_TO_INDEX( world_rank, px,py,pz );

  // Capture global processor decomposition
  g->gnx = gnx;
  g->gny = gny;
  g->gnz = gnz;

  g->dx = (gx1-gx0)/(double)gnx;
  g->dy = (gy1-gy0)/(double)gny;
  g->dz = (gz1-gz0)/(double)gnz;
  g->dV = ((gx1-gx0)/(double)gnx)*
          ((gy1-gy0)/(double)gny)*
          ((gz1-gz0)/(double)gnz);

  g->rdx =  (double)gnx/(gx1-gx0);
  g->rdy =  (double)gny/(gy1-gy0);
  g->rdz =  (double)gnz/(gz1-gz0);
  g->r8V = ((double)gnx/(gx1-gx0))*
           ((double)gny/(gy1-gy0))*
           ((double)gnz/(gz1-gz0))*0.125;

  int nx0 = (gnx * px) / gpx;
  int ny0 = (gny * py) / gpy;
  int nz0 = (gnz * pz) / gpz;

  int nx1 = (gnx * (px+1)) / gpx;
  int ny1 = (gny * (py+1)) / gpy;
  int nz1 = (gnz * (pz+1)) / gpz;

  g->nx0 = nx0;
  g->ny0 = ny0;
  g->nz0 = nz0;

  f = (double)nx0 / (double)gnx; g->x0 = gx0*(1-f) + gx1*f;
  f = (double)ny0 / (double)gny; g->y0 = gy0*(1-f) + gy1*f;
  f = (double)nz0 / (double)gnz; g->z0 = gz0*(1-f) + gz1*f;

  f = (double)nx1 / (double)gnx; g->x1 = gx0*(1-f) + gx1*f;
  f = (double)ny1 / (double)gny; g->y1 = gy0*(1-f) + gy1*f;
  f = (double)nz1 / (double)gnz; g->z1 = gz0*(1-f) + gz1*f;

  // Size the local grid
  size_grid(g,nx1-nx0,ny1-ny0,nz1-nz0);

  // Join the grid to neighbors
  INDEX_TO_RANK(px-1,py,  pz,  rank); join_grid(g,BOUNDARY((-1), 0, 0),rank);
  INDEX_TO_RANK(px,  py-1,pz,  rank); join_grid(g,BOUNDARY( 0,(-1), 0),rank);
  INDEX_TO_RANK(px,  py,  pz-1,rank); join_grid(g,BOUNDARY( 0, 0,(-1)),rank);
  INDEX_TO_RANK(px+1,py,  pz,  rank); join_grid(g,BOUNDARY( 1, 0, 0),rank);
  INDEX_TO_RANK(px,  py+1,pz,  rank); join_grid(g,BOUNDARY( 0, 1, 0),rank);
  INDEX_TO_RANK(px,  py,  pz+1,rank); join_grid(g,BOUNDARY( 0, 0, 1),rank);
}

void
partition_absorbing_box( grid_t * g,
                         double gx0, double gy0, double gz0,
                         double gx1, double gy1, double gz1,
                         int gnx, int gny, int gnz,
                         int gpx, int gpy, int gpz,
                         int pbc ) {
  int px, py, pz;

  partition_periodic_box( g,
                          gx0, gy0, gz0,
                          gx1, gy1, gz1,
                          gnx, gny, gnz,
                          gpx, gpy, gpz );

  // Override periodic boundary conditions

  RANK_TO_INDEX( world_rank, px,py,pz );

  if( px==0 && gnx>1 ) {
    set_fbc(g,BOUNDARY((-1),0,0),absorb_fields);
    set_pbc(g,BOUNDARY((-1),0,0),pbc);
  }

  if( px==gpx-1 && gnx>1 ) {
    set_fbc(g,BOUNDARY( 1,0,0),absorb_fields);
    set_pbc(g,BOUNDARY( 1,0,0),pbc);
  }

  if( py==0 && gny>1 ) {
    set_fbc(g,BOUNDARY(0,(-1),0),absorb_fields);
    set_pbc(g,BOUNDARY(0,(-1),0),pbc);
  }

  if( py==gpy-1 && gny>1 ) {
    set_fbc(g,BOUNDARY(0, 1,0),absorb_fields);
    set_pbc(g,BOUNDARY(0, 1,0),pbc);
  }

  if( pz==0 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0,(-1)),absorb_fields);
    set_pbc(g,BOUNDARY(0,0,(-1)),pbc);
  }

  if( pz==gpz-1 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0, 1),absorb_fields);
    set_pbc(g,BOUNDARY(0,0, 1),pbc);
  }
}

// FIXME: HANDLE 1D and 2D SIMULATIONS IN PARTITION_METAL_BOX
// FIXME: ALLOW USER TO SPECIFIC PBC TO USE ON BOX

void
partition_metal_box( grid_t * g,
                     double gx0, double gy0, double gz0,
                     double gx1, double gy1, double gz1,
                     int gnx, int gny, int gnz,
                     int gpx, int gpy, int gpz ) {
  int px, py, pz;

  partition_periodic_box( g,
                          gx0, gy0, gz0,
                          gx1, gy1, gz1,
                          gnx, gny, gnz,
                          gpx, gpy, gpz );

  // Override periodic boundary conditions

  RANK_TO_INDEX( world_rank, px,py,pz );

  if( px==0 && gnx>1 ) {
    set_fbc(g,BOUNDARY((-1),0,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY((-1),0,0),reflect_particles);
  }

  if( px==gpx-1 && gnx>1 ) {
    set_fbc(g,BOUNDARY(1,0,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(1,0,0),reflect_particles);
  }

  if( py==0 && gny>1 ) {
    set_fbc(g,BOUNDARY(0,(-1),0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,(-1),0),reflect_particles);
  }

  if( py==gpy-1 && gny>1 ) {
    set_fbc(g,BOUNDARY(0,1,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,1,0),reflect_particles);
  }

  if( pz==0 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0,(-1)),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,0,(-1)),reflect_particles);
  }

  if( pz==gpz-1 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0,1),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,0,1),reflect_particles);
  }
}

void
partition_nonuniform_periodic_box( grid_t *g,
                                   double gx0, double gy0, double gz0,
                                   double gx1, double gy1, double gz1,
                                   int gnx, int gny, int gnz,
                                   int gpx, int gpy, int gpz,
                                   Cost& c ) {
  double f;
  int rank, px, py, pz;

  // Make sure the grid can be setup

  if( !g ) ERROR(( "NULL grid" ));

  if( gpx<1 || gpy<1 || gpz<1 || gpx*gpy*gpz!=world_size )
    ERROR(( "Bad domain decompostion (%ix%ix%i)", gpx, gpy, gpz ));

  if( gnx<1 || gny<1 || gnz<1 || gnx%gpx!=0 || gny%gpy!=0 || gnz%gpz!=0 )
    ERROR(( "Bad resolution (%ix%ix%i) for domain decomposition",
            gnx, gny, gnz, gpx, gpy, gpz ));

  // Setup basic variables
  RANK_TO_INDEX( world_rank, px,py,pz );

  // Capture global processor decomposition
  g->gnx = gnx;
  g->gny = gny;
  g->gnz = gnz;

  g->dx = (gx1-gx0)/(double)gnx;
  g->dy = (gy1-gy0)/(double)gny;
  g->dz = (gz1-gz0)/(double)gnz;
  g->dV = ((gx1-gx0)/(double)gnx)*
          ((gy1-gy0)/(double)gny)*
          ((gz1-gz0)/(double)gnz);

  g->rdx =  (double)gnx/(gx1-gx0);
  g->rdy =  (double)gny/(gy1-gy0);
  g->rdz =  (double)gnz/(gz1-gz0);
  g->r8V = ((double)gnx/(gx1-gx0))*
           ((double)gny/(gy1-gy0))*
           ((double)gnz/(gz1-gz0))*0.125;

  // We assume that calling the cost function is cheap, so we will call it
  // three times for every grid cell instead of storing the full 3d profile of
  // per-cell costs in ram
  float* xcost = new float[gnx+2];
  float* ycost = new float[gny+2];
  float* zcost = new float[gnz+2];

  // Prefix sum of cost, integrated over y and z. (Should we take max instead of sum?)
  xcost[0] = 0.;
  for(int i = 0; i < gnx; i++) {
    float sum = 0.;
    for(int j = 0; j < gny; j++) {
      for(int k = 0; k < gnz; k++) {
        sum += c.cost(i,j,k);
      }
    }
    xcost[i+1] = xcost[i] + sum;
  }
  xcost[gnx+1] = 1./0.;

  ycost[0] = 0.;
  for(int j = 0; j < gny; j++) {
    float sum = 0.;
    for(int i = 0; i < gny; i++) {
      for(int k = 0; k < gnz; k++) {
        sum += c.cost(i,j,k);
      }
    }
    ycost[j+1] = ycost[j] + sum;
  }
  ycost[gny+1] = 1./0.;

  zcost[0] = 0.;
  for(int k = 0; k < gnz; k++) {
    float sum = 0.;
    for(int i = 0; i < gnz; i++) {
      for(int j = 0; j < gny; j++) {
        sum += c.cost(i,j,k);
      }
    }
    zcost[k+1] = zcost[k] + sum;
  }
  zcost[gnz+1] = 1./0.;

  // Scan for fair boundaries in x
  int nx0 = 0;
  int i = 0;
  while(xcost[i]/xcost[gnx] <= px/float(gpx)) {
    nx0 = i;
    i++;
  }
  int nx1 = 0;
  i = 0;
  while(xcost[i]/xcost[gnx] <= (px+1)/float(gpx)) {
    nx1 = i;
    i++;
  }

  // Scan for fair boundaries in y
  int ny0 = 0;
  int j = 0;
  while(ycost[j]/ycost[gny] <= py/float(gpy)) {
    ny0 = j;
    j++;
  }
  int ny1 = 0;
  j = 0;
  while(ycost[j]/ycost[gny] <= (py+1)/float(gpy)) {
    ny1 = j;
    j++;
  }

  // Scan for fair boundaries in z
  int nz0 = 0;
  int k = 0;
  while(zcost[k]/zcost[gnz] <= pz/float(gpz)) {
    nz0 = k;
    k++;
  }
  int nz1 = 0;
  k = 0;
  while(zcost[k]/zcost[gnz] <= (pz+1)/float(gpz)) {
    nz1 = k;
    k++;
  }

  delete[] xcost;
  delete[] ycost;
  delete[] zcost;

  g->nx0 = nx0;
  g->ny0 = ny0;
  g->nz0 = nz0;

  f = (double)nx0 / (double)gnx; g->x0 = gx0*(1-f) + gx1*f;
  f = (double)ny0 / (double)gny; g->y0 = gy0*(1-f) + gy1*f;
  f = (double)nz0 / (double)gnz; g->z0 = gz0*(1-f) + gz1*f;

  f = (double)nx1 / (double)gnx; g->x1 = gx0*(1-f) + gx1*f;
  f = (double)ny1 / (double)gny; g->y1 = gy0*(1-f) + gy1*f;
  f = (double)nz1 / (double)gnz; g->z1 = gz0*(1-f) + gz1*f;

  // Size the local grid
  size_grid(g,nx1-nx0,ny1-ny0,nz1-nz0);

  // Join the grid to neighbors
  INDEX_TO_RANK(px-1,py,  pz,  rank); join_grid(g,BOUNDARY((-1), 0, 0),rank);
  INDEX_TO_RANK(px,  py-1,pz,  rank); join_grid(g,BOUNDARY( 0,(-1), 0),rank);
  INDEX_TO_RANK(px,  py,  pz-1,rank); join_grid(g,BOUNDARY( 0, 0,(-1)),rank);
  INDEX_TO_RANK(px+1,py,  pz,  rank); join_grid(g,BOUNDARY( 1, 0, 0),rank);
  INDEX_TO_RANK(px,  py+1,pz,  rank); join_grid(g,BOUNDARY( 0, 1, 0),rank);
  INDEX_TO_RANK(px,  py,  pz+1,rank); join_grid(g,BOUNDARY( 0, 0, 1),rank);
}

void
partition_nonuniform_absorbing_box( grid_t * g,
                                    double gx0, double gy0, double gz0,
                                    double gx1, double gy1, double gz1,
                                    int gnx, int gny, int gnz,
                                    int gpx, int gpy, int gpz,
                                    int pbc,
                                    Cost& c ) {
  int px, py, pz;

  partition_nonuniform_periodic_box( g,
                                     gx0, gy0, gz0,
                                     gx1, gy1, gz1,
                                     gnx, gny, gnz,
                                     gpx, gpy, gpz,
                                     c );

  // Override periodic boundary conditions

  RANK_TO_INDEX( world_rank, px,py,pz );

  if( px==0 && gnx>1 ) {
    set_fbc(g,BOUNDARY((-1),0,0),absorb_fields);
    set_pbc(g,BOUNDARY((-1),0,0),pbc);
  }

  if( px==gpx-1 && gnx>1 ) {
    set_fbc(g,BOUNDARY( 1,0,0),absorb_fields);
    set_pbc(g,BOUNDARY( 1,0,0),pbc);
  }

  if( py==0 && gny>1 ) {
    set_fbc(g,BOUNDARY(0,(-1),0),absorb_fields);
    set_pbc(g,BOUNDARY(0,(-1),0),pbc);
  }

  if( py==gpy-1 && gny>1 ) {
    set_fbc(g,BOUNDARY(0, 1,0),absorb_fields);
    set_pbc(g,BOUNDARY(0, 1,0),pbc);
  }

  if( pz==0 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0,(-1)),absorb_fields);
    set_pbc(g,BOUNDARY(0,0,(-1)),pbc);
  }

  if( pz==gpz-1 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0, 1),absorb_fields);
    set_pbc(g,BOUNDARY(0,0, 1),pbc);
  }
}

void
partition_nonuniform_metal_box( grid_t * g,
                                double gx0, double gy0, double gz0,
                                double gx1, double gy1, double gz1,
                                int gnx, int gny, int gnz,
                                int gpx, int gpy, int gpz,
                                Cost& c) {
  int px, py, pz;

  partition_nonuniform_periodic_box( g,
                                     gx0, gy0, gz0,
                                     gx1, gy1, gz1,
                                     gnx, gny, gnz,
                                     gpx, gpy, gpz,
                                     c );

  // Override periodic boundary conditions

  RANK_TO_INDEX( world_rank, px,py,pz );

  if( px==0 && gnx>1 ) {
    set_fbc(g,BOUNDARY((-1),0,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY((-1),0,0),reflect_particles);
  }

  if( px==gpx-1 && gnx>1 ) {
    set_fbc(g,BOUNDARY(1,0,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(1,0,0),reflect_particles);
  }

  if( py==0 && gny>1 ) {
    set_fbc(g,BOUNDARY(0,(-1),0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,(-1),0),reflect_particles);
  }

  if( py==gpy-1 && gny>1 ) {
    set_fbc(g,BOUNDARY(0,1,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,1,0),reflect_particles);
  }

  if( pz==0 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0,(-1)),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,0,(-1)),reflect_particles);
  }

  if( pz==gpz-1 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0,1),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,0,1),reflect_particles);
  }
}

