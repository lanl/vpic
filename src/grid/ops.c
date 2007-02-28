/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <grid.h>

#define LOCAL_CELL_ID(x,y,z)  INDEX_FORTRAN_3(x,y,z,0,lnx+1,0,lny+1,0,lnz+1)
#define REMOTE_CELL_ID(x,y,z) INDEX_FORTRAN_3(x,y,z,0,rnx+1,0,rny+1,0,rnz+1)

/* Everybody must size their local grid in parallel */

void size_grid( grid_t * RESTRICT g, int lnx, int lny, int lnz ) {
  int rank, nproc, i, j, k, x, y, z, lnc;
  INT64_TYPE ii, jj, kk; 

  if( g==NULL )                 { ERROR(("Bad grid"));            return; }
  if( g->mp==NULL             ) { ERROR(("Bad mp"));              return; }
  if( lnx<1 || lny<1 || lnz<1 ) { ERROR(("Bad local grid size")); return; }

  rank = mp_rank(g->mp);
  nproc = mp_nproc(g->mp);

  /* Setup phase 2 data structures */
  g->nx = lnx;
  g->ny = lny;
  g->nz = lnz;
  for( k=-1; k<=1; k++ )
    for( j=-1; j<=1; j++ )
      for( i=-1; i<=1; i++ ) 
        g->bc[ BOUNDARY(i,j,k) ] = pec_fields;
  g->bc[ BOUNDARY(0,0,0) ] = rank;

  /* Setup phase 3 data structures */
  /* This is an ugly kludge to interface phase 2 and phase 3 data structures */
  lnc = (lnx+2)*(lny+2)*(lnz+2);
  if( g->range!=NULL ) free_aligned( g->range );
  g->range = (INT64_TYPE * ALIGNED)
    malloc_aligned( (nproc+1)*sizeof(INT64_TYPE), preferred_alignment );
  if( g->range==NULL ) {
    ERROR(("Could not allocate range array"));
    return;
  }
  ii = lnc; /* lnc is not 64-bits */
  mp_allgather_i64(&ii,g->range,1,g->mp);
  jj = 0;
  g->range[nproc] = 0;
  for( i=0; i<=nproc; i++ ) {
    kk = g->range[i];
    g->range[i] = jj;
    jj += kk;
  }
  g->rangel = g->range[rank];
  g->rangeh = g->range[rank+1]-1;

  if( g->neighbor!=NULL ) free_aligned( g->neighbor );
  g->neighbor = (INT64_TYPE * ALIGNED)
    malloc_aligned( 6*lnc*sizeof(INT64_TYPE), preferred_alignment );

  if( g->neighbor==NULL ) {
    ERROR(("Could not allocate neighbor array"));
    return;
  }
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
        /* Set boundary faces appropriately */
        if( x==1   ) g->neighbor[i+0] = reflect_particles;
        if( y==1   ) g->neighbor[i+1] = reflect_particles;
        if( z==1   ) g->neighbor[i+2] = reflect_particles;
        if( x==lnx ) g->neighbor[i+3] = reflect_particles;
        if( y==lny ) g->neighbor[i+4] = reflect_particles;
        if( z==lnz ) g->neighbor[i+5] = reflect_particles;
        /* Set ghost cells appropriately */
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
}

void join_grid( grid_t * RESTRICT g, int boundary, int rank ) {
  int lx, ly, lz, lnx, lny, lnz, rx, ry, rz, rnx, rny, rnz, rnc;

  if( g==NULL                   ) { ERROR(("Bad grid"));     return; }
  if( g->mp==NULL               ) { ERROR(("Bad mp"));       return; }
  if( boundary<0 ||
      boundary>=27 ||
      boundary==BOUNDARY(0,0,0) ) { ERROR(("Bad boundary")); return; }
  if( rank<0 ||
      rank>=mp_nproc(g->mp)     ) { ERROR(("Bad rank"));     return; }

  /* Join phase 2 data structures */
  g->bc[boundary] = rank;

  /* Join phase 3 data structures */
  lnx = g->nx;
  lny = g->ny;
  lnz = g->nz;
  rnc = g->range[rank+1] - g->range[rank]; /* NOTE: rnc <~ 2^31 / 6 */


# define GLUE_FACE(tag,i,j,k,X,Y,Z) BEGIN_PRIMITIVE {           \
    if( boundary==BOUNDARY(i,j,k) ) {                           \
      if( rnc%((ln##Y+2)*(ln##Z+2))!=0 ) {                      \
        ERROR(("Remote face is incompatible"));                 \
        return;                                                 \
      }                                                         \
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
}

void set_fbc( grid_t *g, int boundary, int fbc ) {
  if( g==NULL                       ) { ERROR(("Bad grid"));      return; }
  if( g->mp==NULL                   ) { ERROR(("Bad mp"));        return; }
  if( boundary<0 ||
      boundary>=27 ||
      boundary==BOUNDARY(0,0,0)     ) { ERROR(("Bad boundary"));  return; }
  if( fbc>=0 &&
      fbc<mp_nproc(g->mp)           ) { ERROR(("Use join_grid")); return; }
  if( fbc!=anti_symmetric_fields &&
      fbc!=symmetric_fields      &&
      fbc!=pmc_fields            &&
      fbc!=absorb_fields            ) { ERROR(("Bad field bc"));  return; }
  g->bc[boundary] = fbc;
}

void set_pbc( grid_t *g, int boundary, int pbc ) {
  int lx, ly, lz, lnx, lny, lnz;

  if( g==NULL                       ) { ERROR(("Bad grid"));        return; }
  if( g->mp==NULL                   ) { ERROR(("Bad mp"));          return; }
  if( boundary<0 ||
      boundary>=27 ||
      boundary==BOUNDARY(0,0,0)     ) { ERROR(("Bad boundary"));    return; }
  if( pbc>=0 && pbc<mp_nproc(g->mp) ) { ERROR(("Use join_grid"));   return; }

/* FIXME:  Fix this error condition to work with custom boundary conditions */
#if 0
  if( pbc!=absorb_particles &&
      pbc!=reflect_particles        ) { ERROR(("Bad particle bc")); return; }
#endif

  if( pbc!=absorb_particles && pbc!=reflect_particles && 
      (-pbc-3<0 || -pbc-3>=g->nb) ) { 
    ERROR(("Bad particle bc")); 
    return; 
  }

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
}


