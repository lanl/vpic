/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#ifndef _grid_h_
#define _grid_h_

#include "../util/util.h"

// Define a "pointer to boundary handler function" type.

// FIXME: THESE PREDECLARATIONS ARE UGLY AND MAKE BABY JESUS CRY

struct particle;
struct particle_mover;
struct field;
struct accumulator;
struct grid;
struct species;
struct particle_injector;

typedef int
(*boundary_handler_t)( void                     * params,
                       struct particle          * r,
                       struct particle_mover    * pm,       
                       struct field             * f,
                       struct accumulator       * a,
                       const struct grid        * g,
                       struct species           * s, 
                       struct particle_injector * ppi, 
                       mt_rng_t                 * rng, 
                       int                        face );

enum boundary_handler_enums {
  INVALID_BOUNDARY       = 0xBADF00D,
  MAX_BOUNDARY_DATA_SIZE = 1024 // Sized to hold data characterizing a boundary model.
};

// FIXME: MODEL_PARAMETER DOES NOT HAVE GUARANTEED NICE ALIGNMENT

typedef struct boundary {
 boundary_handler_t handler;
 char params[MAX_BOUNDARY_DATA_SIZE]; 
} boundary_t;

#define BOUNDARY(i,j,k) (13+(i)+3*(j)+9*(k)) /* FORTRAN -1:1,-1:1,-1:1 */

enum grid_enums {

  // Phase 2 boundary conditions
  anti_symmetric_fields = -1, // E_tang = 0
  pec_fields            = -1,
  metal_fields          = -1,
  symmetric_fields      = -2, // B_tang = 0, B_norm = 0
  pmc_fields            = -3, // B_tang = 0, B_norm floats
  absorb_fields         = -4, // Gamma = 0

  // Phase 3 boundary conditions
  reflect_particles = -1, // Cell boundary should reflect particles
  absorb_particles  = -2  // Cell boundary should absorb particles

  // Symmetry in the field boundary conditions refers to image charge
  // sign
  //
  // Anti-symmetric -> Image charges are opposite signed (ideal metal)
  //                   Boundary rho/j are accumulated over partial voxel+image
  // Symmetric      -> Image charges are same signed (symmetry plane or pmc)
  //                   Boundary rho/j are accumulated over partial voxel+image
  // Absorbing      -> No image charges
  //                   Boundary rho/j are accumulated over partial voxel only
  //
  // rho     -> Anti-symmetric      | rho     -> Symmetric
  // jf_tang -> Anti-symmetric      | jf_tang -> Symmetric
  // E_tang  -> Anti-symmetric      | E_tang  -> Symmetric
  // B_norm  -> Anti-symmetric + DC | B_norm  -> Symmetric      (see note)
  // B_tang  -> Symmetric           | B_tang  -> Anti-symmetric
  // E_norm  -> Symmetric           | E_norm  -> Anti-symmetric (see note)
  // div B   -> Symmetric           | div B   -> Anti-symmetric
  // 
  // Note: B_norm is tricky. For a symmetry plane, B_norm on the
  // boundary must be zero as there are no magnetic charges (a
  // non-zero B_norm would imply an infinitesimal layer of magnetic
  // charge). However, if a symmetric boundary is interpreted as a
  // perfect magnetic conductor, B_norm could be present due to
  // magnetic conduction surface charges. Even though there are no
  // bulk volumetric magnetic charges to induce a surface magnetic
  // charge, I think that radiation/waveguide modes/etc could (the
  // total surface magnetic charge in the simulation would be zero
  // though). As a result, symmetric and pmc boundary conditions are
  // treated separately. Symmetric and pmc boundaries are identical
  // except the symmetric boundaries explicitly zero boundary
  // B_norm. Note: anti-symmetric and pec boundary conditions would
  // have the same issue if norm E was located directly on the
  // boundary. However, it is not so this problem does not arise.
  //
  // Note: Absorbing boundary conditions make no effort to clean
  // divergence errors on them. They assume that the ghost div b is
  // zero and force the surface div e on them to be zero. This means
  // ghost norm e can be set to any value on absorbing boundaries.

};

typedef struct grid {
  mp_handle mp;           // Communications handle
  float dt, cvac, eps0;   // System of units
  float damp;             // Radiation damping parameter
                          // FIXME: DOESN'T GO HERE

  // Phase 2 grid data structures 
  float x0, y0, z0;       // Min corner local domain (must be coherent)
  float x1, y1, z1;       // Max corner local domain (must be coherent)
  float dx, dy, dz;       // Cell dimensions (CONVENIENCE ... USE
                          // x0,x1 WHEN DECIDING WHICH NODE TO USE!)
  float rdx, rdy, rdz;    // Inverse voxel dimensions (CONVENIENCE)
  int   nx, ny, nz;       // Local voxel mesh resolution.  Voxels are
                          // indexed FORTRAN style 0:nx+1,0:ny+1,0:nz+1
                          // with voxels 1:nx,1:ny,1:nz being non-ghost
                          // voxels.
  int   bc[27];           // (-1:1,-1:1,-1:1) FORTRAN indexed array of
                          // boundary conditions to apply at domain edge
                          // 0 ... nproc-1 ... comm boundary condition
                          // <0 ... locally applied boundary condition

  // Phase 3 grid data structures

  // NOTE: VOXEL INDEXING LIMITS NUMBER OF VOXELS TO 2^31 (INCLUDING
  // GHOSTS) PER NODE.  NEIGHBOR INDEXING FURTHER LIMITS TO
  // (2^31)/6.  THE EMITTER COMPONENT ID INDEXING FURTHER LIMITS TO
  // 2^27 PER NODE.  THE LIMIT IS 2^64 OVER ALL NODES THOUGH.

  int64_t * ALIGNED(16) range;
                          // (0:nproc) indexed array giving range of
                          // global indexes of voxel owned by each
                          // processor.  Replicated on each processor.
                          // (range[rank]:range[rank+1]-1) are global
                          // voxels owned by processor "rank".  Note:
                          // range[rank+1]-range[rank] <~ 2^31 / 6

  int64_t * ALIGNED(128) neighbor;
                          // (0:5,0:local_num_voxel-1) FORTRAN indexed
                          // array neighbor(0:5,lidx) are the global
                          // indexes of neighboring voxels of the
                          // voxel with local index "lidx".  Negative
                          // if neighbor is a boundary condition.

  int64_t rangel, rangeh; // Redundant for move_p performance reasons:
                          //   rangel = range[rank]
                          //   rangeh = range[rank+1]-1.
                          // Note: rangeh-rangel <~ 2^31 / 6

  // Enable user-defined boundary handlers

  int nb;                 // Number of custom boundary conditions
  boundary_t * boundary;  // Head of array of boundary_t.
                          // FIXME: DOESN'T GO HERE!

} grid_t;

// Given a voxel mesh coordinates (on 0:nx+1,0:ny+1,0:nz+1) and
// voxel mesh resolution (nx,ny,nz), return the index of that voxel.

#define VOXEL(x,y,z, nx,ny,nz) ((x) + ((nx)+2)*((y) + ((ny)+2)*(z)))

// Advance the voxel mesh index (v) and corresponding voxel mesh
// coordinates (x,y,z) in a region with min- and max-corners of
// (xl,yl,zl) and (xh,yh,zh) of a (nx,ny,nz) resolution voxel mesh in
// FORTRAN ordering.  Results will not be valid (v,x,y,z) are not in
// the region or if (v,x,y,z) is the last voxel in that region.
//
// This macro is not robust.  Macro arguments should be safe against
// multiple evaluation.  Further, this macro is not semantically a
// single statement.  (It is meant for use in high performance stencil
// inner loops.)
//
// This is written with seeming extraneously if tests in order to get
// the compiler to generate branceless conditional move and add 
// instructions (none of the branches below are actual branches in
// assembly).

#define NEXT_VOXEL(v,x,y,z, xl,xh, yl,yh, zl,zh, nx,ny,nz) \
  (v)++;                                                   \
  (x)++;                                                   \
  if( (x)>(xh) ) (v) +=  (nx)-(xh)+(xl)+1;                 \
  if( (x)>(xh) ) (y)++;                                    \
  if( (x)>(xh) ) (x) = (xl);                               \
  if( (y)>(yh) ) (v) += ((ny)-(yh)+(yl)+1)*((nx)+2);       \
  if( (y)>(yh) ) (z)++;                                    \
  if( (y)>(yh) ) (y) = (yl)

BEGIN_C_DECLS

// In boundary_handler.c

// Note that boundaries _copy_ the state given by ip into their own
// storage.  Thus, to change the state of a given boundary, one has to
// look up the boundary and change it there.  Example usage:
//
//   reflux_boundary_params_t bc;
//
//   bc.species[0] = ion;
//   bc.ux[0]      = sqrt(k*T/(m*c^2))
//   ...
//   reflux_boundary = add_boundary( g, reflux_handler, &bc );
// FIXME: WHY MAKE USER PASS &bc ... SILLY
#define add_boundary(g,bh,ip) IUO_add_boundary((g),(bh),(ip),sizeof(*(ip)))

int
IUO_add_boundary( grid_t *g,
		  boundary_handler_t bh,
		  const void * initial_params,
		  int sz );

// In grid_structors.c

grid_t *
new_grid( void );

void
delete_grid( grid_t * g );

// In ops.c

void
size_grid( grid_t * g, int lnx, int lny, int lnz );

void
join_grid( grid_t * g, int bound, int rank );

void
set_fbc( grid_t *g, int bound, int fbc );

void
set_pbc( grid_t *g, int bound, int pbc );

// In partition.c

// g->{n,d}{x,y,z} is _coherent_ on all nodes in the domain after
// these calls as are g->{x,y,z}{0,1}.  Due to the vagaries of
// floating point, though g->nx*g->dx may not be the exactly the same
// as g->x1-g->x0 though.  Thus matters when doing things like
// robustly converting global position coordinates to/from local index
// + offset position coordinates.
//
// The robust procedure to convert _from_ a global coordinate to a
// local coordinate is:
//
// (1) Test if this node has ownership of the point using
// g->{x,y,z}{0,1}.  Points with x==g->x1 exactly boundaries should be
// considered part of the local domain only if the corresponding
// x-boundary condition is local.  Similarly for y- and z-.
//
// (2) If this node has ownership of the point, compute the relative
// voxel and offset of the x-coordinate via
// g->nx*((x-g->x0)/(g->x1-g->x0)), _NOT_ (x-g->x0)/g->dx and _NOT_
// (x-g->x0)*(1/g->dx)!  Similarly for y and z.
//
// (3) Break the voxel and offsets into integer and fractional parts.
// Particles exactly on the far wall should have their fractional
// particles set to 1 and their integer parts subtracted by 1.  Double
// the fractional part and subtract by one to get the voxel centered
// offset.  Convert the local voxel coordinates into a local voxel index
// using VOXEL above.
//
// Reverse this protocol to robustly convert from voxel+offset to
// global coordinates.  Due to the vagaries of floating point, the
// inverse process may not be exact.

void
partition_periodic_box( grid_t *g,
			double gx0, double gy0, double gz0,
			double gx1, double gy1, double gz1,
                        int gnx, int gny, int gnz,
                        int gpx, int gpy, int gpz );

void
partition_absorbing_box( grid_t *g,
                         double gx0, double gy0, double gz0,
                         double gx1, double gy1, double gz1,
                         int gnx, int gny, int gnz,
                         int gpx, int gpy, int gpz,
                         int pbc );

void
partition_metal_box( grid_t *g,
                     double gx0, double gy0, double gz0,
                     double gx1, double gy1, double gz1,
                     int gnx, int gny, int gnz,
                     int gpx, int gpy, int gpz );

// In grid_comm.c

// FIXME: SHOULD TAKE A RAW PORT INDEX INSTEAD OF A PORT COORDS

// Start receiving a message from the node.
// Only one message recv may be pending at a time on a given port.

void
begin_recv_port( int i,    // x port coord ([-1,0,1])
                 int j,    // y port coord ([-1,0,1])
                 int k,    // z port coord ([-1,0,1])
                 int size, // Expected size in bytes
                 const grid_t * g );

// Returns pointer to the buffer that begin send will use for the next
// send on the given port.  The buffer is guaranteed to have enough
// room for size bytes.  This is only valid to call if no sends on
// that port are pending.

void * ALIGNED(16)
size_send_port( int i,    // x port coord ([-1,0,1])
                int j,    // y port coord ([-1,0,1])
                int k,    // z port coord ([-1,0,1])
                int size, // Needed send size in bytes
                const grid_t * g );

// Begin sending size bytes of the buffer out the given port.  Only
// one message send may be pending at a time on a given port.  (FIXME:
// WHAT HAPPENS IF SIZE_SEND_PORT size < begin_send_port
// size??)

void
begin_send_port( int i,    // x port coord ([-1,0,1])
                 int j,    // y port coord ([-1,0,1])
                 int k,    // z port coord ([-1,0,1])
                 int size, // Number of bytes to send (in bytes)
                 const grid_t * g );

// Complete the pending recv on the given port.  Only valid to call if
// there is a pending recv.  Returns pointer to a buffer containing
// the received data.  (FIXME: WHAT HAPPENS IF EXPECTED RECV SIZE
// GIVEN IN BEGIN_RECV DOES NOT MATCH END_RECV??)

void * ALIGNED(16)
end_recv_port( int i, // x port coord ([-1,0,1])
               int j, // y port coord ([-1,0,1])
               int k, // z port coord ([-1,0,1])
               const grid_t * g );

// Complete the pending send on the given port.  Only valid to call if
// there is a pending send on the port.  Note that this guarantees
// that send port is available to the caller for additional use, not
// necessarily that the message has arrived at the destination of the
// port.

void
end_send_port( int i, // x port coord ([-1,0,1])
               int j, // y port coord ([-1,0,1])
               int k, // z port coord ([-1,0,1])
               const grid_t * g );

// In distribute_voxels.c

// Given a block of voxels to be processed, determine the number of
// the first voxel (v,x,y,z) a particular job assigned to a pipeline
// should process and return the number of voxels to process.
//
// It is assumed that the pipelines will process voxels in FORTRAN
// ordering (e.g. inner loop increments x-index).
//
// jobs are indexed from 0 to n_job-1.  jobs are _always_ have the
// number of voxels an integer multiple of the bundle size.  If job 
// is set to n_job, this function will determine the parameters of
// the final incomplete bundle.
//
// FIXME: To acccomodate PPE and SPE builds, the distribute_voxels
// function is deprecated.  Use the macro instead.

int
distribute_voxels( int x0, int x1,     // range of x-indices (inclusive)
                   int y0, int y1,     // range of y-indices (inclusive)
                   int z0, int z1,     // range of z-indices (inclusive)
                   int bundle,         // number of voxels in a bundle
                   int job, int n_job, // job ... on [0,n_job-1]
                   int * _x, int * _y, int * _z );

#define DISTRIBUTE_VOXELS( x0,x1, y0,y1, z0,z1, b, p,P, x,y,z,nv ) do { \
    int _x0=(x0), _y0=(y0), _z0=(z0), _b=(b), _p=(p), _P=(P);           \
    int _nx = (x1)-_x0+1, _ny = (y1)-_y0+1, _nv = _nx*_ny*((z1)-_z0+1); \
    double _t = (double)( _nv/_b ) / (double)_P;                        \
    int          _x=_b*(int)( _t*(double)(_p  ) + 0.5 ), _y, _z;        \
    if( _p<_P ) _nv=_b*(int)( _t*(double)(_p+1) + 0.5 );                \
    _nv -= _x;     /* x = (x-x0) + nx*( (y-y0) + ny*(z-z0) ) */         \
    _y   = _x/_nx; /* y =               (y-y0) + ny*(z-z0)   */         \
    _z   = _y/_ny; /* z =                           (z-z0)   */         \
    _x  -= _y*_nx; /* x = (x-x0)                             */         \
    _y  -= _z*_ny; /* y =               (y-y0)               */         \
    (x)  = _x+_x0;                                                      \
    (y)  = _y+_y0;                                                      \
    (z)  = _z+_z0;                                                      \
    (nv) = _nv;                                                         \
  } while(0)

END_C_DECLS

#endif
