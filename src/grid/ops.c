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

#define LOCAL_CELL_ID(x,y,z)  VOXEL(x,y,z, lnx,lny,lnz)
#define REMOTE_CELL_ID(x,y,z) VOXEL(x,y,z, rnx,rny,rnz)

// Everybody must size their local grid in parallel

void
size_grid( grid_t * g,
           int lnx, int lny, int lnz ) {
  int64_t x,y,z;
  int i, j, k;
  int64_t ii, jj, kk;

  if( !g || lnx<1 || lny<1 || lnz<1 ) ERROR(( "Bad args" ));

  // Setup phase 2 data structures
  g->sx =  1;
  g->sy = (lnx+2)*g->sx;
  g->sz = (lny+2)*g->sy;
  g->nv = (lnz+2)*g->sz;
  g->nx = lnx; g->ny = lny; g->nz = lnz;
  for( k=-1; k<=1; k++ )
    for( j=-1; j<=1; j++ )
      for( i=-1; i<=1; i++ )
        g->bc[ BOUNDARY(i,j,k) ] = pec_fields;
  g->bc[ BOUNDARY(0,0,0) ] = world_rank;

  // Setup phase 3 data structures.  This is an ugly kludge to
  // interface phase 2 and phase 3 data structures
  FREE_ALIGNED( g->range );
  MALLOC_ALIGNED( g->range, world_size+1, 16 );
  ii = g->nv; // nv is not 64-bits
  mp_allgather_i64( &ii, g->range, 1 );
  jj = 0;
  g->range[world_size] = 0;
  for( i=0; i<=world_size; i++ ) {
    kk = g->range[i];
    g->range[i] = jj;
    jj += kk;
  }
  g->rangel = g->range[world_rank];
  g->rangeh = g->range[world_rank+1]-1;

  FREE_ALIGNED( g->neighbor );
  MALLOC_ALIGNED( g->neighbor, 6*g->nv, 128 );

  for( z=0; z<=lnz+1; z++ )
    for( y=0; y<=lny+1; y++ )
      for( x=0; x<=lnx+1; x++ ) {
        i = 6*LOCAL_CELL_ID(x,y,z);
        g->neighbor[i+0] = g->rangel + LOCAL_CELL_ID(x-1,y,z);
        g->neighbor[i+1] = g->rangel + LOCAL_CELL_ID(x,y-1,z);
        g->neighbor[i+2] = g->rangel + LOCAL_CELL_ID(x,y,z-1);
        g->neighbor[i+3] = g->rangel + LOCAL_CELL_ID(x+1,y,z);
        g->neighbor[i+4] = g->rangel + LOCAL_CELL_ID(x,y+1,z);
        g->neighbor[i+5] = g->rangel + LOCAL_CELL_ID(x,y,z+1);
        // Set boundary faces appropriately
        if( x==1   ) g->neighbor[i+0] = reflect_particles;
        if( y==1   ) g->neighbor[i+1] = reflect_particles;
        if( z==1   ) g->neighbor[i+2] = reflect_particles;
        if( x==lnx ) g->neighbor[i+3] = reflect_particles;
        if( y==lny ) g->neighbor[i+4] = reflect_particles;
        if( z==lnz ) g->neighbor[i+5] = reflect_particles;
        // Set ghost cells appropriately
        if( x==0 || x==lnx+1 ||
            y==0 || y==lny+1 ||
            z==0 || z==lnz+1 ) {
          g->neighbor[i+0] = reflect_particles;
          g->neighbor[i+1] = reflect_particles;
          g->neighbor[i+2] = reflect_particles;
          g->neighbor[i+3] = reflect_particles;
          g->neighbor[i+4] = reflect_particles;
          g->neighbor[i+5] = reflect_particles;
        }
      }

# if 0
  // Setup the space filling curve
  // FIXME: THIS IS A CRUDE HACK UNTIL A GOOD SFC CAN BE WRITTEN
  // CURRENT SFC IS GOOD FOR THREADING WITH POWER-OF-TWO NUMBER OF THREADS.
  // UP TO AND INCLUDING 8 THREADS.

  FREE_ALIGNED( g->sfc );
  MALLOC_ALIGNED( g->sfc, g->nv, 128 );

  do {
    int off;
    int ox, oy, oz;
    int nx, ny, nz;
    int nx1 = (lnx+2)/2,   ny1 = (lny+2)/2,   nz1 = (lnz+2)/2;
    int nx0 = (lnx+2)-nx1, ny0 = (lny+2)-ny1, nz0 = (lnz+2)-nz1;

    for( z=0; z<=lnz+1; z++ )
      for( y=0; y<=lny+1; y++ )
        for( x=0; x<=lnx+1; x++ ) {
          i = LOCAL_CELL_ID(x,y,z);
          off = 0;
          ox=x; nx=nx0; if(ox>=nx) ox-=nx, off+=nx,                 nx=nx1;
          oy=y; ny=ny0; if(oy>=ny) oy-=ny, off+=ny*(lnx+2),         ny=ny1;
          oz=z; nz=nz0; if(oz>=nz) oz-=nz, off+=nz*(lnx+2)*(lny+2), nz=nz1;
          g->sfc[i] = off + ox + nx*( oy + ny*oz );
        }
  } while(0);
# endif
}

void
join_grid( grid_t * g,
           int boundary,
           int rank ) {
  int lx, ly, lz, lnx, lny, lnz, rx, ry, rz, rnx, rny, rnz, rnc;

  if( !g || boundary<0 || boundary>=27 || boundary==BOUNDARY(0,0,0) ||
      rank<0 || rank>=world_size ) ERROR(( "Bad args" ));

  // Join phase 2 data structures
  g->bc[boundary] = rank;

  // Join phase 3 data structures
  lnx = g->nx;
  lny = g->ny;
  lnz = g->nz;
  rnc = g->range[rank+1] - g->range[rank]; // Note: rnc <~ 2^31 / 6

# define GLUE_FACE(tag,i,j,k,X,Y,Z) BEGIN_PRIMITIVE {           \
    if( boundary==BOUNDARY(i,j,k) ) {                           \
      if( rnc%((ln##Y+2)*(ln##Z+2))!=0 )                        \
        ERROR(("Remote face is incompatible"));                 \
      rn##X = (rnc/((ln##Y+2)*(ln##Z+2)))-2;                    \
      rn##Y = ln##Y;                                            \
      rn##Z = ln##Z;                                            \
      for( l##Z=1; l##Z<=ln##Z; l##Z++ ) {                      \
        for( l##Y=1; l##Y<=ln##Y; l##Y++ ) {                    \
          l##X = (i+j+k)<0 ? 1     : ln##X;                     \
          r##X = (i+j+k)<0 ? rn##X : 1;                         \
          r##Y = l##Y;                                          \
          r##Z = l##Z;                                          \
          g->neighbor[ 6*LOCAL_CELL_ID(lx,ly,lz) + tag ] =      \
            g->range[rank] + REMOTE_CELL_ID(rx,ry,rz);          \
        }                                                       \
      }                                                         \
      return;                                                   \
    }                                                           \
  } END_PRIMITIVE

  GLUE_FACE(0,-1, 0, 0,x,y,z);
  GLUE_FACE(1, 0,-1, 0,y,z,x);
  GLUE_FACE(2, 0, 0,-1,z,x,y);
  GLUE_FACE(3, 1, 0, 0,x,y,z);
  GLUE_FACE(4, 0, 1, 0,y,z,x);
  GLUE_FACE(5, 0, 0, 1,z,x,y);

# undef GLUE_FACE
}

void
set_fbc( grid_t * g,
         int boundary,
         int fbc ) {

  if( !g || boundary<0 || boundary>=27 || boundary==BOUNDARY(0,0,0) ||
      ( fbc!=anti_symmetric_fields && fbc!=symmetric_fields &&
        fbc!=pmc_fields            && fbc!=absorb_fields    ) )
    ERROR(( "Bad args" ));

  g->bc[boundary] = fbc;
}

void
set_pbc( grid_t * g,
         int boundary,
         int pbc ) {
  int lx, ly, lz, lnx, lny, lnz;

  if( !g || boundary<0 || boundary>=27 || boundary==BOUNDARY(0,0,0) || pbc>=0 )
    ERROR(( "Bad args" ));

  lnx = g->nx;
  lny = g->ny;
  lnz = g->nz;

# define SET_PBC(tag,i,j,k,X,Y,Z) BEGIN_PRIMITIVE {             \
    if( boundary==BOUNDARY(i,j,k) ) {                           \
      l##X = (i+j+k)<0 ? 1 : ln##X;                             \
      for( l##Z=1; l##Z<=ln##Z; l##Z++ )                        \
        for( l##Y=1; l##Y<=ln##Y; l##Y++ )                      \
          g->neighbor[ 6*LOCAL_CELL_ID(lx,ly,lz) + tag ] = pbc; \
      return;                                                   \
    }                                                           \
  } END_PRIMITIVE

  SET_PBC(0,-1, 0, 0,x,y,z);
  SET_PBC(1, 0,-1, 0,y,z,x);
  SET_PBC(2, 0, 0,-1,z,x,y);
  SET_PBC(3, 1, 0, 0,x,y,z);
  SET_PBC(4, 0, 1, 0,y,z,x);
  SET_PBC(5, 0, 0, 1,z,x,y);

# undef SET_PBC
}

