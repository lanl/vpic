/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "sf_interface.h"

#define hydro(x,y,z) hydro[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

// Generic looping
#define XYZ_LOOP(xl,xh,yl,yh,zl,zh) \
  for( z=zl; z<=zh; z++ )	    \
    for( y=yl; y<=yh; y++ )	    \
      for( x=xl; x<=xh; x++ )
	      
// x_NODE_LOOP => Loop over all non-ghost nodes at plane x
#define x_NODE_LOOP(x) XYZ_LOOP(x,x,1,ny+1,1,nz+1)
#define y_NODE_LOOP(y) XYZ_LOOP(1,nx+1,y,y,1,nz+1)
#define z_NODE_LOOP(z) XYZ_LOOP(1,nx+1,1,ny+1,z,z)

// Note: synchronize_hydro assumes that hydro has not been adjusted at
// the local domain boundary to account for partial cells

void
synchronize_hydro( hydro_t      * ALIGNED(128) hydro,
                   const grid_t *              g ) {
  int size, face, x, y, z, nx, ny, nz;
  float *p, lw, rw;
  hydro_t *h;

  if( hydro==NULL ) ERROR(("Bad hydro"));
  if( g==NULL     ) ERROR(("Bad grid"));

  local_adjust_hydro( hydro, g );

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;

# define BEGIN_RECV(i,j,k,X,Y,Z) \
  begin_recv_port(i,j,k,( 1 + 14*(n##Y+1)*(n##Z+1) )*sizeof(float),g)

# define BEGIN_SEND(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {      \
    size = ( 1 + 14*(n##Y+1)*(n##Z+1) )*sizeof(float);  \
    p = (float *)size_send_port( i, j, k, size, g );    \
    if( p!=NULL ) {                                     \
      (*(p++)) = g->d##X;                               \
      face = (i+j+k)<0 ? 1 : n##X+1;                    \
      X##_NODE_LOOP(face) {                             \
        h = &hydro(x,y,z);                              \
        (*(p++)) = h->jx;                               \
        (*(p++)) = h->jy;                               \
        (*(p++)) = h->jz;                               \
        (*(p++)) = h->rho;                              \
        (*(p++)) = h->px;                               \
        (*(p++)) = h->py;                               \
        (*(p++)) = h->pz;                               \
        (*(p++)) = h->ke;                               \
        (*(p++)) = h->txx;                              \
        (*(p++)) = h->tyy;                              \
        (*(p++)) = h->tzz;                              \
        (*(p++)) = h->tyz;                              \
        (*(p++)) = h->tzx;                              \
        (*(p++)) = h->txy;                              \
      }                                                 \
      begin_send_port( i, j, k, size, g );              \
    }                                                   \
  } END_PRIMITIVE

# define END_RECV(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {                \
    p = (float *)end_recv_port(i,j,k,g);                        \
    if( p!=NULL ) {                                             \
      rw = (*(p++));                 /* Remote g->d##X */       \
      lw = rw + g->d##X;                                        \
      rw /= lw;                                                 \
      lw = g->d##X/lw;                                          \
      lw += lw;                                                 \
      rw += rw;                                                 \
      face = (i+j+k)<0 ? n##X+1 : 1; /* Twice weighted sum */   \
      X##_NODE_LOOP(face) {                                     \
        h = &hydro(x,y,z);                                      \
        h->jx  = lw*h->jx  + rw*(*(p++));                       \
        h->jy  = lw*h->jy  + rw*(*(p++));                       \
        h->jz  = lw*h->jz  + rw*(*(p++));                       \
        h->rho = lw*h->rho + rw*(*(p++));                       \
        h->px  = lw*h->px  + rw*(*(p++));                       \
        h->py  = lw*h->py  + rw*(*(p++));                       \
        h->pz  = lw*h->pz  + rw*(*(p++));                       \
        h->ke  = lw*h->ke  + rw*(*(p++));                       \
        h->txx = lw*h->txx + rw*(*(p++));                       \
        h->tyy = lw*h->tyy + rw*(*(p++));                       \
        h->tzz = lw*h->tzz + rw*(*(p++));                       \
        h->tyz = lw*h->tyz + rw*(*(p++));                       \
        h->tzx = lw*h->tzx + rw*(*(p++));                       \
        h->txy = lw*h->txy + rw*(*(p++));                       \
      }                                                         \
    }                                                           \
  } END_PRIMITIVE

# define END_SEND(i,j,k,X,Y,Z) end_send_port( i, j, k, g )

  // Exchange x-faces
  BEGIN_SEND(-1, 0, 0,x,y,z);
  BEGIN_SEND( 1, 0, 0,x,y,z);
  BEGIN_RECV(-1, 0, 0,x,y,z);
  BEGIN_RECV( 1, 0, 0,x,y,z);
  END_RECV(-1, 0, 0,x,y,z);
  END_RECV( 1, 0, 0,x,y,z);
  END_SEND(-1, 0, 0,x,y,z);
  END_SEND( 1, 0, 0,x,y,z);

  // Exchange y-faces
  BEGIN_SEND( 0,-1, 0,y,z,x);
  BEGIN_SEND( 0, 1, 0,y,z,x);
  BEGIN_RECV( 0,-1, 0,y,z,x);
  BEGIN_RECV( 0, 1, 0,y,z,x);
  END_RECV( 0,-1, 0,y,z,x);
  END_RECV( 0, 1, 0,y,z,x);
  END_SEND( 0,-1, 0,y,z,x);
  END_SEND( 0, 1, 0,y,z,x);

  // Exchange z-faces
  BEGIN_SEND( 0, 0,-1,z,x,y);
  BEGIN_SEND( 0, 0, 1,z,x,y);
  BEGIN_RECV( 0, 0,-1,z,x,y);
  BEGIN_RECV( 0, 0, 1,z,x,y);
  END_RECV( 0, 0,-1,z,x,y);
  END_RECV( 0, 0, 1,z,x,y);
  END_SEND( 0, 0,-1,z,x,y);
  END_SEND( 0, 0, 1,z,x,y);

# undef BEGIN_RECV
# undef BEGIN_SEND
# undef END_RECV
# undef END_SEND
}

// Because hydro fields are purely diagnostic, correct the hydro along
// local boundaries for particle cell accumulations

void
local_adjust_hydro( hydro_t      * ALIGNED(128) hydro,
                    const grid_t *              g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  int bc, face, x, y, z;
  hydro_t *h;

# define ADJUST_HYDRO(i,j,k,X,Y,Z)              \
  do {                                          \
    bc = g->bc[BOUNDARY(i,j,k)];                \
    if( bc<0 || bc>nproc) {                     \
      face = (i+j+k)<0 ? 1 : n##X+1;            \
      X##_NODE_LOOP(face) {                     \
        h = &hydro(x,y,z);                      \
        h->jx  *= 2;                            \
        h->jy  *= 2;                            \
        h->jz  *= 2;                            \
        h->rho *= 2;                            \
        h->px  *= 2;                            \
        h->py  *= 2;                            \
        h->pz  *= 2;                            \
        h->ke  *= 2;                            \
        h->txx *= 2;                            \
        h->tyy *= 2;                            \
        h->tzz *= 2;                            \
        h->tyz *= 2;                            \
        h->tzx *= 2;                            \
        h->txy *= 2;                            \
      }                                         \
    }                                           \
  } while(0)
  
  ADJUST_HYDRO(-1, 0, 0,x,y,z);
  ADJUST_HYDRO( 0,-1, 0,y,z,x);
  ADJUST_HYDRO( 0, 0,-1,z,x,y);
  ADJUST_HYDRO( 1, 0, 0,x,y,z);
  ADJUST_HYDRO( 0, 1, 0,y,z,x);
  ADJUST_HYDRO( 0, 0, 1,z,x,y);
}
