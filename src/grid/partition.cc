/*
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <iostream> // TODO: delete

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



/**
 * @brief Gives an array one bigger than the input, over the sum of numbers
 * before that index
 *
 * @param data Data to sum over
 *
 * @return Sum
 */
std::vector<int> running_sum(std::vector<int> data)
{
    std::vector<int> sum;
    sum.push_back( 0 );

    // Turn these into running sums
    for (int i = 1; i <= data.size(); i++)
    {
        sum[i] = sum[i-1] + data[i-1];
        std::cout << "Sum at " << i << " = " << sum[i] << std::endl;
    }

    return sum;
}


void
partition_periodic_box( grid_t * g,
                        double gx0, double gy0, double gz0,
                        double gx1, double gy1, double gz1,
                        int gnx, int gny, int gnz,
                        int gpx, int gpy, int gpz,
                        std::vector<int> nx_array,
                        std::vector<int> ny_array,
                        std::vector<int> nz_array
                )

{
  double f;
  int rank, px, py, pz;

  std::cout << "Global grid x "
      << gx0 << ".." << gx1 << ", "
      << gy0 << ".." << gy1 << ", "
      << gz0 << ".." << gz1 << ", "
      << std::endl;
  // Make sure the grid can be setup

  if( !g ) ERROR(( "NULL grid" ));

  if( gpx<1 || gpy<1 || gpz<1 || gpx*gpy*gpz!=world_size )
  {
    ERROR(( "Bad domain decompostion (%ix%ix%i)", gpx, gpy, gpz ));
  }

  //if( gnx<1 || gny<1 || gnz<1 || gnx%gpx!=0 || gny%gpy!=0 || gnz%gpz!=0 )
    //ERROR(( "Bad resolution (%ix%ix%i) for domain decomposition",
            //gnx, gny, gnz, gpx, gpy, gpz ));
  if (nx_array.size() < 1 || ny_array.size() < 1 || nz_array.size() < 1 ||
          nx_array.size() != gpx  ||
          ny_array.size() != gpy  ||
          nz_array.size() != gpz
     )
  {
    ERROR(( "Missmatched topology sizes (%ix%ix%i) for decomposition (%ix%ix%i)",
            nx_array.size(), ny_array.size(), nz_array.size(), gpx, gpy, gpz ));
  }

  std::vector<int> nx_sum = running_sum(nx_array);
  std::vector<int> ny_sum = running_sum(ny_array);
  std::vector<int> nz_sum = running_sum(nz_array);

  // Setup basic variables
  RANK_TO_INDEX( world_rank, px,py,pz );

  // TODO: from here it could call the base function?
  int lnx = nx_array[ px ];
  int lny = ny_array[ py ];
  int lnz = nz_array[ pz ];

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

  // Find what percent we are into the domain
  f = (double) nx_sum[px] / (double)gnx;
  g->x0 = gx0*(1-f) + gx1*f;
  f = (double) ny_sum[py]   /(double)gny;
  g->y0 = gy0*(1-f) + gy1*f;
  f = (double) nz_sum[pz]   /(double)gnz;
  g->z0 = gz0*(1-f) + gz1*f;

  f = (double) nx_sum[px+1] /(double)gnx;
  g->x1 = gx0*(1-f) + gx1*f;
  f = (double) ny_sum[py+1] /(double)gny;
  g->y1 = gy0*(1-f) + gy1*f;

  f = (double) nz_sum[pz+1] /(double)gnz;
  g->z1 = gz0*(1-f) + gz1*f;
  std::cout << "zf " << f << " gnz " << gnz << " px " << px << "nzpx1 " << nz_sum[px+1] << std::endl;


  std::cout << world_rank << " goes from " << g->x0 << " to " << g->x1 << std::endl;
  std::cout << world_rank << " goes from " << g->y0 << " to " << g->y1 << std::endl;
  std::cout << world_rank << " goes from " << g->z0 << " to " << g->z1 << std::endl;

  // Size the local grid
  //size_grid(g,gnx/gpx,gny/gpy,gnz/gpz);
  size_grid(g, lnx, lny, lnz);

  // Join the grid to neighbors
  INDEX_TO_RANK(px-1,py,  pz,  rank); join_grid(g,BOUNDARY((-1), 0, 0),rank);
  INDEX_TO_RANK(px,  py-1,pz,  rank); join_grid(g,BOUNDARY( 0,(-1), 0),rank);
  INDEX_TO_RANK(px,  py,  pz-1,rank); join_grid(g,BOUNDARY( 0, 0,(-1)),rank);
  INDEX_TO_RANK(px+1,py,  pz,  rank); join_grid(g,BOUNDARY( 1, 0, 0),rank);
  INDEX_TO_RANK(px,  py+1,pz,  rank); join_grid(g,BOUNDARY( 0, 1, 0),rank);
  INDEX_TO_RANK(px,  py,  pz+1,rank); join_grid(g,BOUNDARY( 0, 0, 1),rank);

}

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

  f = (double) px   /(double)gpx; g->x0 = gx0*(1-f) + gx1*f;
  f = (double) py   /(double)gpy; g->y0 = gy0*(1-f) + gy1*f;
  f = (double) pz   /(double)gpz; g->z0 = gz0*(1-f) + gz1*f;

  f = (double)(px+1)/(double)gpx; g->x1 = gx0*(1-f) + gx1*f;
  f = (double)(py+1)/(double)gpy; g->y1 = gy0*(1-f) + gy1*f;
  f = (double)(pz+1)/(double)gpz; g->z1 = gz0*(1-f) + gz1*f;

  // Size the local grid
  size_grid(g,gnx/gpx,gny/gpy,gnz/gpz);

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
