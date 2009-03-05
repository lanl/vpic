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

#define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {  \
  int _ix, _iy, _iz;                                    \
  _ix  = (rank);  /* ix = ix + gpx*( iy + gpy*iz ) */   \
  _iy  = _ix/gpx; /* iy = iy + gpy*iz */                \
  _ix -= _iy*gpx; /* ix = ix */                         \
  _iz  = _iy/gpy; /* iz = iz */                         \
  _iy -= _iz*gpy; /* iy = iy */                         \
  (ix) = _ix;                                           \
  (iy) = _iy;                                           \
  (iz) = _iz;                                           \
} END_PRIMITIVE

#define INDEX_TO_RANK(ix,iy,iz,rank) BEGIN_PRIMITIVE {  \
  int _ix, _iy, _iz;                                    \
   /* Wrap processor index periodically */              \
   _ix = (ix) % gpx; while(_ix<0) _ix += gpx;           \
   _iy = (iy) % gpy; while(_iy<0) _iy += gpy;           \
   _iz = (iz) % gpz; while(_iz<0) _iz += gpz;           \
   /* Compute the rank */                               \
  (rank) = _ix + gpx*( _iy + gpy*_iz );                 \
} END_PRIMITIVE

void
partition_periodic_box( grid_t * g,
                        double gx0, double gy0, double gz0,
                        double gx1, double gy1, double gz1,
                        int gnx, int gny, int gnz,
                        int gpx, int gpy, int gpz ) {
  double f;
  int rank, nproc, px, py, pz; 

  // Make sure the grid can be setup
  if( g==NULL                      ) ERROR(("Bad grid"));
  if( g->mp==NULL                  ) ERROR(("Bad mp"));
  if( gpx<1 || gpy<1 || gpz<1 ||
      gpx*gpy*gpz!=mp_nproc(g->mp) ) ERROR(("Bad topology"));
  if( gnx<1 || gny<1 || gnz<1      ) ERROR(("Bad res"));
  if( gnx%gpx!=0 ||
      gny%gpy!=0 ||
      gnz%gpz!=0                   ) ERROR(("Incompatible res"));
  rank = mp_rank(g->mp);
  nproc = mp_nproc(g->mp);

  // Setup basic variables
  RANK_TO_INDEX(rank,px,py,pz);

  // keep track of topology for file I/O
  g->gpx = gpx;
  g->gpy = gpy;
  g->gpz = gpz;

  g->dx = (gx1-gx0)/(double)gnx;
  g->dy = (gy1-gy0)/(double)gny;
  g->dz = (gz1-gz0)/(double)gnz;

  g->rdx = (double)gnx/(gx1-gx0);
  g->rdy = (double)gny/(gy1-gy0);
  g->rdz = (double)gnz/(gz1-gz0);

  f = (double) px   /(double)gpx; g->x0 = gx0*(1-f) + gx1*f;
  f = (double) py   /(double)gpy; g->y0 = gy0*(1-f) + gy1*f;
  f = (double) pz   /(double)gpz; g->z0 = gz0*(1-f) + gz1*f;

  f = (double)(px+1)/(double)gpx; g->x1 = gx0*(1-f) + gx1*f;
  f = (double)(py+1)/(double)gpy; g->y1 = gy0*(1-f) + gy1*f;
  f = (double)(pz+1)/(double)gpz; g->z1 = gz0*(1-f) + gz1*f;

  // Size the local grid
  size_grid(g,gnx/gpx,gny/gpy,gnz/gpz);

  // Join the grid to neighbors
  INDEX_TO_RANK(px-1,py,  pz,  rank); join_grid(g,BOUNDARY(-1, 0, 0),rank);
  INDEX_TO_RANK(px,  py-1,pz,  rank); join_grid(g,BOUNDARY( 0,-1, 0),rank);
  INDEX_TO_RANK(px,  py,  pz-1,rank); join_grid(g,BOUNDARY( 0, 0,-1),rank);
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
  int rank, nproc, px, py, pz; 

  partition_periodic_box( g,
                          gx0, gy0, gz0,
                          gx1, gy1, gz1,
                          gnx, gny, gnz,
                          gpx, gpy, gpz );

  // Override periodic boundary conditions

  rank = mp_rank(g->mp);
  nproc = mp_nproc(g->mp);
  RANK_TO_INDEX(rank,px,py,pz);

  // keep track of topology for file I/O
  g->gpx = gpx;
  g->gpy = gpy;
  g->gpz = gpz;

  if( px==0 && gnx>1 ) { 
    set_fbc(g,BOUNDARY(-1,0,0),absorb_fields);
    set_pbc(g,BOUNDARY(-1,0,0),pbc);
  } 
  if( px==gpx-1 && gnx>1 ) {
    set_fbc(g,BOUNDARY( 1,0,0),absorb_fields);
    set_pbc(g,BOUNDARY( 1,0,0),pbc);
  }

  if( py==0 && gny>1 ) { 
    set_fbc(g,BOUNDARY(0,-1,0),absorb_fields);
    set_pbc(g,BOUNDARY(0,-1,0),pbc);
  } 

  if( py==gpy-1 && gny>1 ) {
    set_fbc(g,BOUNDARY(0, 1,0),absorb_fields);
    set_pbc(g,BOUNDARY(0, 1,0),pbc);
  }

  if( pz==0 && gnz>1 ) { 
    set_fbc(g,BOUNDARY(0,0,-1),absorb_fields);
    set_pbc(g,BOUNDARY(0,0,-1),pbc);
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
  int rank, nproc, px, py, pz; 

  partition_periodic_box( g,
                          gx0, gy0, gz0,
                          gx1, gy1, gz1,
                          gnx, gny, gnz,
                          gpx, gpy, gpz );

  // Override periodic boundary conditions

  rank = mp_rank(g->mp);
  nproc = mp_nproc(g->mp);
  RANK_TO_INDEX(rank,px,py,pz);

  // keep track of topology for file I/O
  g->gpx = gpx;
  g->gpy = gpy;
  g->gpz = gpz;

  if( px==0 && gnx>1 ) {
    set_fbc(g,BOUNDARY(-1,0,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(-1,0,0),reflect_particles);
  }

  if( px==gpx-1 && gnx>1 ) {
    set_fbc(g,BOUNDARY(1,0,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(1,0,0),reflect_particles);
  }

  if( py==0 && gny>1 ) {
    set_fbc(g,BOUNDARY(0,-1,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,-1,0),reflect_particles);
  }

  if( py==gpy-1 && gny>1 ) {
    set_fbc(g,BOUNDARY(0,1,0),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,1,0),reflect_particles);
  }

  if( pz==0 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0,-1),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,0,-1),reflect_particles);
  }

  if( pz==gpz-1 && gnz>1 ) {
    set_fbc(g,BOUNDARY(0,0,1),anti_symmetric_fields);
    set_pbc(g,BOUNDARY(0,0,1),reflect_particles);
  }
}
