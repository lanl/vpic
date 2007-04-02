/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

/******************************************************************************
 * local.c sets local boundary conditions. Functions are divided into two
 * categories:
 *   local_ghosts_xxx where xxx = tang_b, norm_e, div_b
 *   - Sets ghosts values of the fields just interior to a local boundary
 *     condition
 *   local_adjust_xxx where xxx = norm_b, tang_e, rhof, rhob, div_e_err
 *   - Directly enforces local boundary conditions on fields
 *****************************************************************************/
#include <field.h>

#define f(x,y,z)         f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
#define hydro(x,y,z) hydro[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

#define XYZ_LOOP(xl,xh,yl,yh,zl,zh)		\
  for( z=zl; z<=zh; z++ )			\
    for( y=yl; y<=yh; y++ )			\
      for( x=xl; x<=xh; x++ )

#define yz_EDGE_LOOP(x) XYZ_LOOP(x,x,1,ny,1,nz+1)
#define zx_EDGE_LOOP(y) XYZ_LOOP(1,nx+1,y,y,1,nz)
#define xy_EDGE_LOOP(z) XYZ_LOOP(1,nx,1,ny+1,z,z)

#define zy_EDGE_LOOP(x) XYZ_LOOP(x,x,1,ny+1,1,nz)
#define xz_EDGE_LOOP(y) XYZ_LOOP(1,nx,y,y,1,nz+1)
#define yx_EDGE_LOOP(z) XYZ_LOOP(1,nx+1,1,ny,z,z)

#define x_NODE_LOOP(x) XYZ_LOOP(x,x,1,ny+1,1,nz+1)
#define y_NODE_LOOP(y) XYZ_LOOP(1,nx+1,y,y,1,nz+1)
#define z_NODE_LOOP(z) XYZ_LOOP(1,nx+1,1,ny+1,z,z)

#define x_FACE_LOOP(x) XYZ_LOOP(x,x,1,ny,1,nz)
#define y_FACE_LOOP(y) XYZ_LOOP(1,nx,y,y,1,nz)
#define z_FACE_LOOP(z) XYZ_LOOP(1,nx,1,ny,z,z)

/*****************************************************************************
 * Local ghosts
 *****************************************************************************/

void local_ghost_tang_b( field_t * ALIGNED f,
			 const grid_t * g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  const float cdt_dx = g->cvac*g->dt/g->dx;
  const float cdt_dy = g->cvac*g->dt/g->dy;
  const float cdt_dz = g->cvac*g->dt/g->dz;
  int bc, face, ghost, x, y, z;
  float decay, drive, higend, t1, t2;
  field_t *fg, *fh;

  /* Absorbing boundary condition is 2nd order accurate implementation of a
     1st order Higend ABC with 15 degree annihilation cone except for 1d
     simulations where the 2nd order accurate implementation of a 1st order
     Mur boundary condition is used. */
  higend = ( nx>1 || ny>1 || nz>1 ) ? 1.03527618 : 1.;

# define APPLY_LOCAL_TANG_B(i,j,k,X,Y,Z)                                 \
  do {                                                                   \
    bc = g->bc[BOUNDARY(i,j,k)];                                         \
    if( bc<0 || bc>nproc) {                                              \
      ghost = (i+j+k)<0 ? 0 : n##X+1;                                    \
      face  = (i+j+k)<0 ? 1 : n##X+1;                                    \
      switch(bc) {                                                       \
      case anti_symmetric_fields:                                        \
	Z##Y##_EDGE_LOOP(ghost) f(x,y,z).cb##Y= f(x-i,y-j,z-k).cb##Y;    \
	Y##Z##_EDGE_LOOP(ghost) f(x,y,z).cb##Z= f(x-i,y-j,z-k).cb##Z;    \
	break;                                                           \
      case symmetric_fields: case pmc_fields:                            \
	Z##Y##_EDGE_LOOP(ghost) f(x,y,z).cb##Y=-f(x-i,y-j,z-k).cb##Y;    \
	Y##Z##_EDGE_LOOP(ghost) f(x,y,z).cb##Z=-f(x-i,y-j,z-k).cb##Z;    \
	break;                                                           \
      case absorb_fields:                                                \
        drive = cdt_d##X*higend;                                         \
        decay = (1-drive)/(1+drive);                                     \
        drive = 2*drive/(1+drive);                                       \
	Z##Y##_EDGE_LOOP(ghost) {                                        \
          fg = &f(x,y,z);                                                \
          fh = &f(x-i,y-j,z-k);                                          \
          X = face;                                                      \
          t1 = cdt_d##X*( f(x-i,y-j,z-k).e##Z - f(x,y,z).e##Z );         \
          t1 = (i+j+k)<0 ? t1 : -t1;                                     \
          X = ghost;                                                     \
          Z++; t2 = f(x-i,y-j,z-k).e##X;                                 \
          Z--; t2 = cdt_d##Z*( t2 - fh->e##X );                          \
          fg->cb##Y = decay*fg->cb##Y + drive*fh->cb##Y - t1 + t2;       \
        }                                                                \
	Y##Z##_EDGE_LOOP(ghost) {                                        \
          fg = &f(x,y,z);                                                \
          fh = &f(x-i,y-j,z-k);                                          \
          X = face;                                                      \
          t1 = cdt_d##X*( f(x-i,y-j,z-k).e##Y - f(x,y,z).e##Y );         \
          t1 = (i+j+k)<0 ? t1 : -t1;                                     \
          X = ghost;                                                     \
          Y++; t2 = f(x-i,y-j,z-k).e##X;                                 \
          Y--; t2 = cdt_d##Y*( t2 - fh->e##X );                          \
          fg->cb##Z = decay*fg->cb##Z + drive*fh->cb##Z + t1 - t2;       \
        }                                                                \
	break;                                                           \
      default:                                                           \
	ERROR(("Bad boundary condition encountered."));                  \
	break;                                                           \
      }                                                                  \
    }                                                                    \
  } while(0)

  APPLY_LOCAL_TANG_B(-1, 0, 0,x,y,z);
  APPLY_LOCAL_TANG_B( 0,-1, 0,y,z,x);
  APPLY_LOCAL_TANG_B( 0, 0,-1,z,x,y);
  APPLY_LOCAL_TANG_B( 1, 0, 0,x,y,z);
  APPLY_LOCAL_TANG_B( 0, 1, 0,y,z,x);
  APPLY_LOCAL_TANG_B( 0, 0, 1,z,x,y);
}

/* Note: local_adjust_div_e zeros the error on the boundaries for absorbing
   boundary conditions. Thus, ghost norm e value is irrevelant */

void local_ghost_norm_e( field_t * ALIGNED f,
                         const grid_t * g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  int bc, face, x, y, z;
  field_t * ALIGNED f0, * ALIGNED f1, * ALIGNED f2;

# define APPLY_LOCAL_NORM_E(i,j,k,X,Y,Z)                        \
  do {                                                          \
    bc = g->bc[BOUNDARY(i,j,k)];                                \
    if( bc<0 || bc>nproc) {                                     \
      face = (i+j+k)<0 ? 0 : n##X+1;                            \
      switch(bc) {                                              \
      case anti_symmetric_fields:                               \
	X##_NODE_LOOP(face) {                                   \
          f0 = &f(x,y,z);                                       \
          f1 = &f(x-i,y-j,z-k);                                 \
          f0->e##X   = f1->e##X;                                \
          f0->tca##X = f1->tca##X;                              \
        }                                                       \
	break;                                                  \
      case symmetric_fields: case pmc_fields:                   \
	X##_NODE_LOOP(face) {                                   \
          f0 = &f(x,y,z);                                       \
          f1 = &f(x-i,y-j,z-k);                                 \
          f0->e##X   = -f1->e##X;                               \
          f0->tca##X = -f1->tca##X;                             \
        }                                                       \
	break;                                                  \
      case absorb_fields:                                       \
	X##_NODE_LOOP(face) {                                   \
          f0 = &f(x,y,z);                                       \
          f1 = &f(x-i,y-j,z-k);                                 \
          f2 = &f(x-i*2,y-j*2,z-k*2);                           \
          f0->e##X   = 2*f1->e##X   - f2->e##X;                 \
          f0->tca##X = 2*f1->tca##X - f2->tca##X;               \
        }                                                       \
	break;                                                  \
      default:                                                  \
	ERROR(("Bad boundary condition encountered."));         \
	break;                                                  \
      }                                                         \
    }                                                           \
  } while(0)

  APPLY_LOCAL_NORM_E(-1, 0, 0,x,y,z);
  APPLY_LOCAL_NORM_E( 0,-1, 0,y,z,x);
  APPLY_LOCAL_NORM_E( 0, 0,-1,z,x,y);
  APPLY_LOCAL_NORM_E( 1, 0, 0,x,y,z);
  APPLY_LOCAL_NORM_E( 0, 1, 0,y,z,x);
  APPLY_LOCAL_NORM_E( 0, 0, 1,z,x,y);
}

void local_ghost_div_b( field_t * ALIGNED f,
                        const grid_t * g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  int bc, face, x, y, z;

# define APPLY_LOCAL_DIV_B(i,j,k,X,Y,Z)					    \
  do {									    \
    bc = g->bc[BOUNDARY(i,j,k)];					    \
    if( bc<0 || bc>nproc) {						    \
      face = (i+j+k)<0 ? 0 : n##X+1;					    \
      switch(bc) {							    \
      case anti_symmetric_fields:					    \
	X##_FACE_LOOP(face) f(x,y,z).div_b_err =  f(x-i,y-j,z-k).div_b_err; \
	break;								    \
      case symmetric_fields: case pmc_fields:				    \
	X##_FACE_LOOP(face) f(x,y,z).div_b_err = -f(x-i,y-j,z-k).div_b_err; \
	break;								    \
      case absorb_fields:						    \
	X##_FACE_LOOP(face) f(x,y,z).div_b_err = 0;			    \
	break;								    \
      default:								    \
	ERROR(("Bad boundary condition encountered."));	                    \
	break;								    \
      }									    \
    }									    \
  } while(0)
  
  APPLY_LOCAL_DIV_B(-1, 0, 0,x,y,z);
  APPLY_LOCAL_DIV_B( 0,-1, 0,y,z,x);
  APPLY_LOCAL_DIV_B( 0, 0,-1,z,x,y);
  APPLY_LOCAL_DIV_B( 1, 0, 0,x,y,z);
  APPLY_LOCAL_DIV_B( 0, 1, 0,y,z,x);
  APPLY_LOCAL_DIV_B( 0, 0, 1,z,x,y);
}

/*****************************************************************************
 * Local adjusts
 *****************************************************************************/

/* FIXME: Specialty edge loops should be added to zero e_tang on local edges
   exclusively to handle concave domain geometries */

void local_adjust_tang_e( field_t * ALIGNED f,
                          const grid_t * g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  int bc, face, x, y, z;
  field_t *fs;

# define ADJUST_TANG_E(i,j,k,X,Y,Z)                                     \
  do {                                                                  \
    bc = g->bc[BOUNDARY(i,j,k)];                                        \
    if( bc<0 || bc>nproc) {                                             \
      face = (i+j+k)<0 ? 1 : n##X+1;                                    \
      switch(bc) {                                                      \
      case anti_symmetric_fields:                                       \
	Y##Z##_EDGE_LOOP(face) {                                        \
          fs = &f(x,y,z);                                               \
          fs->e##Y = 0;                                                 \
          fs->tca##Y = 0;                                               \
        }                                                               \
	Z##Y##_EDGE_LOOP(face) {                                        \
          fs = &f(x,y,z);                                               \
          fs->e##Z = 0;                                                 \
          fs->tca##Z = 0;                                               \
        }                                                               \
	break;                                                          \
      case symmetric_fields: case pmc_fields: case absorb_fields:       \
        break;                                                          \
      default:                                                          \
	ERROR(("Bad boundary condition encountered."));                 \
	break;                                                          \
      }                                                                 \
    }                                                                   \
  } while(0)

  ADJUST_TANG_E(-1, 0, 0,x,y,z);
  ADJUST_TANG_E( 0,-1, 0,y,z,x);
  ADJUST_TANG_E( 0, 0,-1,z,x,y);
  ADJUST_TANG_E( 1, 0, 0,x,y,z);
  ADJUST_TANG_E( 0, 1, 0,y,z,x);
  ADJUST_TANG_E( 0, 0, 1,z,x,y);
}

void local_adjust_norm_b( field_t * ALIGNED f,
                          const grid_t * g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  int bc, face, x, y, z;

# define ADJUST_NORM_B(i,j,k,X,Y,Z)                                     \
  do {                                                                  \
    bc = g->bc[BOUNDARY(i,j,k)];                                        \
    if( bc<0 || bc>nproc) {                                             \
      face = (i+j+k)<0 ? 1 : n##X+1;                                    \
      switch(bc) {                                                      \
      case anti_symmetric_fields: case pmc_fields: case absorb_fields:  \
	break;                                                          \
      case symmetric_fields:                                            \
	X##_FACE_LOOP(face) f(x,y,z).cb##X = 0;                         \
	break;                                                          \
      default:                                                          \
	ERROR(("Bad boundary condition encountered."));                 \
	break;                                                          \
      }                                                                 \
    }                                                                   \
  } while(0)

  ADJUST_NORM_B(-1, 0, 0,x,y,z);
  ADJUST_NORM_B( 0,-1, 0,y,z,x);
  ADJUST_NORM_B( 0, 0,-1,z,x,y);
  ADJUST_NORM_B( 1, 0, 0,x,y,z);
  ADJUST_NORM_B( 0, 1, 0,y,z,x);
  ADJUST_NORM_B( 0, 0, 1,z,x,y);
}

void local_adjust_div_e( field_t * ALIGNED f,
                         const grid_t * g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  int bc, face, x, y, z;

# define ADJUST_DIV_E_ERR(i,j,k,X,Y,Z)			         \
  do {							         \
    bc = g->bc[BOUNDARY(i,j,k)];				 \
    if( bc<0 || bc>nproc) {					 \
      face = (i+j+k)<0 ? 1 : n##X+1;				 \
      switch(bc) {						 \
      case anti_symmetric_fields: case absorb_fields:		 \
        X##_NODE_LOOP(face) f(x,y,z).div_e_err = 0;              \
        break;                                                   \
      case symmetric_fields: case pmc_fields:			 \
        break;                                                   \
      default:							 \
	ERROR(("Bad boundary condition encountered."));          \
	break;							 \
      }								 \
    }								 \
  } while(0)

  ADJUST_DIV_E_ERR(-1, 0, 0,x,y,z);
  ADJUST_DIV_E_ERR( 0,-1, 0,y,z,x);
  ADJUST_DIV_E_ERR( 0, 0,-1,z,x,y);
  ADJUST_DIV_E_ERR( 1, 0, 0,x,y,z);
  ADJUST_DIV_E_ERR( 0, 1, 0,y,z,x);
  ADJUST_DIV_E_ERR( 0, 0, 1,z,x,y);
}

/* anti_symmetric => Opposite sign image charges (zero jf_tang)
   symmetric      => Same sign image charges (double jf_tang) 
   absorbing      => No image charges, half cell accumulation (double jf_tang)
   (rhob/jf_norm account for particles that hit boundary and reflect/stick) */

void local_adjust_jf( field_t * ALIGNED f,
			const grid_t * g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  int bc, face, x, y, z;

# define ADJUST_JF(i,j,k,X,Y,Z)                                         \
  do {                                                                  \
    bc = g->bc[BOUNDARY(i,j,k)];                                        \
    if( bc<0 || bc>nproc) {                                             \
      face = (i+j+k)<0 ? 1 : n##X+1;                                    \
      switch(bc) {                                                      \
      case anti_symmetric_fields:                                       \
	Y##Z##_EDGE_LOOP(face) f(x,y,z).jf##Y = 0;                      \
        Z##Y##_EDGE_LOOP(face) f(x,y,z).jf##Z = 0;                      \
	break;                                                          \
      case symmetric_fields: case pmc_fields: case absorb_fields:       \
	Y##Z##_EDGE_LOOP(face) f(x,y,z).jf##Y *= 2.;                    \
        Z##Y##_EDGE_LOOP(face) f(x,y,z).jf##Z *= 2.;                    \
	break;                                                          \
      default:                                                          \
	ERROR(("Bad boundary condition encountered."));                 \
	break;                                                          \
      }                                                                 \
    }                                                                   \
  } while(0)
  
  ADJUST_JF(-1, 0, 0,x,y,z);
  ADJUST_JF( 0,-1, 0,y,z,x);
  ADJUST_JF( 0, 0,-1,z,x,y);
  ADJUST_JF( 1, 0, 0,x,y,z);
  ADJUST_JF( 0, 1, 0,y,z,x);
  ADJUST_JF( 0, 0, 1,z,x,y);
}

/* anti_symmetric => Opposite sign image charges (zero rhof/rhob)
   symmetric      => Same sign image charges (double rhof)
                  => (double rhof, rhob is already correct)
   absorbing      => No image charges, half cell accumulation (double rhof)
   (rhob/jf_norm account for particles that hit the boundary) */

void local_adjust_rhof( field_t * ALIGNED f,
                        const grid_t * g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  int bc, face, x, y, z;

# define ADJUST_RHOF(i,j,k,X,Y,Z)                                       \
  do {                                                                  \
    bc = g->bc[BOUNDARY(i,j,k)];                                        \
    if( bc<0 || bc>nproc) {                                             \
      face = (i+j+k)<0 ? 1 : n##X+1;                                    \
      switch(bc) {                                                      \
      case anti_symmetric_fields:                                       \
	X##_NODE_LOOP(face) f(x,y,z).rhof = 0;                          \
	break;                                                          \
      case symmetric_fields: case pmc_fields: case absorb_fields:       \
	X##_NODE_LOOP(face) f(x,y,z).rhof *= 2;                         \
        break;                                                          \
      default:                                                          \
	ERROR(("Bad boundary condition encountered."));                 \
	break;                                                          \
      }                                                                 \
    }                                                                   \
  } while(0)
  
  ADJUST_RHOF(-1, 0, 0,x,y,z);
  ADJUST_RHOF( 0,-1, 0,y,z,x);
  ADJUST_RHOF( 0, 0,-1,z,x,y);
  ADJUST_RHOF( 1, 0, 0,x,y,z);
  ADJUST_RHOF( 0, 1, 0,y,z,x);
  ADJUST_RHOF( 0, 0, 1,z,x,y);
}

/* anti_symmetric => Opposite sign image charges (zero rhob)
   symmetric      => Same sign image charges (rhob already correct)
   absorbing      => No image charges, half cell accumulation (rhob already
                     correct) */

void local_adjust_rhob( field_t * ALIGNED f,
                        const grid_t * g ) {
  const int nx = g->nx, ny = g->ny, nz = g->nz, nproc = mp_nproc(g->mp);
  int bc, face, x, y, z;

# define ADJUST_RHOB(i,j,k,X,Y,Z)                                       \
  do {                                                                  \
    bc = g->bc[BOUNDARY(i,j,k)];                                        \
    if( bc<0 || bc>nproc) {                                             \
      face = (i+j+k)<0 ? 1 : n##X+1;                                    \
      switch(bc) {                                                      \
      case anti_symmetric_fields:                                       \
	X##_NODE_LOOP(face) f(x,y,z).rhob = 0;                          \
	break;                                                          \
      case symmetric_fields: case pmc_fields: case absorb_fields:       \
        break;                                                          \
      default:                                                          \
	ERROR(("Bad boundary condition encountered."));                 \
	break;                                                          \
      }                                                                 \
    }                                                                   \
  } while(0)
  
  ADJUST_RHOB(-1, 0, 0,x,y,z);
  ADJUST_RHOB( 0,-1, 0,y,z,x);
  ADJUST_RHOB( 0, 0,-1,z,x,y);
  ADJUST_RHOB( 1, 0, 0,x,y,z);
  ADJUST_RHOB( 0, 1, 0,y,z,x);
  ADJUST_RHOB( 0, 0, 1,z,x,y);
}

/* Because hydro fields are purely diagnostic, correct the hydro along local
   boundaries for particle cell accumulations */

void local_adjust_hydro( hydro_t * ALIGNED hydro,
                         const grid_t * g ) {
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
        h->rho *= 2;                            \
        h->jx  *= 2;                            \
        h->jy  *= 2;                            \
        h->jz  *= 2;                            \
        h->ke  *= 2;                            \
        h->px  *= 2;                            \
        h->py  *= 2;                            \
        h->pz  *= 2;                            \
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
