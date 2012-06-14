/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Revised and extended from earlier V4PIC versions
 *
 */

#include <iostream> // For std::cerr and friends
#include "vpic/vpic.hxx"

//-----------------------------------------------------------------------------

#define begin_globals struct user_global_t
#define global ((struct user_global_t *)user_global)

#define begin_initialization                                      \
void                                                              \
vpic_simulation::user_initialization( int num_cmdline_arguments,  \
                                      char ** cmdline_argument )

#define begin_diagnostics \
void                      \
vpic_simulation::user_diagnostics( void ) 

#define begin_particle_injection \
void                             \
vpic_simulation::user_particle_injection( void )

#define begin_current_injection \
void                            \
vpic_simulation::user_current_injection( void )

#define begin_field_injection \
void                          \
vpic_simulation::user_field_injection( void )

#define begin_particle_collisions \
void                              \
vpic_simulation::user_particle_collisions( void )

#define repeat( count ) for( int64_t _remain=(int64_t)(count); _remain; _remain-- )

#define _SIM_LOG_PREFIX \
  __FILE__ "(" EXPAND_AND_STRINGIFY(__LINE__) ")[" << rank() << "]: "
#define sim_log_local(x) std::cerr << _SIM_LOG_PREFIX << x << std::endl
#define sim_log(x) do {                                \
    if( rank()==0 ) {                                  \
      std::cerr << _SIM_LOG_PREFIX << x << std::endl;  \
      std::cerr.flush();                               \
    }                                                  \
  } while(0)

//-----------------------------------------------------------------------------

// These macros provide support for setting materials, boundary
// conditions and field values inside and on the surface of regions.
//
// Most macros work by providing a logical expression in terms of
// double precision coordinates x,y,z that are non-zero if a point is
// inside the intended region and 0 if not. The field macros also take
// several other equations to set field values. For example:
//
// set_region_field( x>0 && sqrt(x*x+y*y+z*z)<1,  // A half-sphere region 
//                   sin(k*x), 0, 0,              // electric field
//                   0, sin(k*x), bz );           // magnetic field
//
// There are two types of regions, point regions and regular regions.
//
// A material value or field component is inside a point region if its
// location is inside the region. A boundary condition face is inside
// a point region if all corners of the face are inside the
// region. Otherwise, a face is on the partially inside a point region
// if some corner of the face is inside the region.
//
// A regular region has two parts: an interior and a
// surface. set_region_bc further divides the surface into an interior
// surface and an exterior surface.  The mapping of the region to the
// grid is dictated by the solely by the location of cell centers.
//
// Interior cells are cells whose centers are inside the
// region. Exterior cells are cells whose centers are outside the
// region.
// 
// Surface faces are faces for which one associated cell-center is
// inside the region. Interior faces are faces where both associated
// cell-centers are inside the region. Interior surface faces are
// faces whose associated cell-center is inside the region but
// neighbor cell-center is outside the region. The exterior surface
// faces are faces whose associated cell-center is outside the region
// but neighbor cell-center is inside the region.
//
// Surface edges are edges for which up to 3 associated cell-centers
// are inside the region. Interior edges are edges where all
// associated cell-centers are inside the region.
//
// Surface nodes are nodes for which up to 7 one associated
// cell-centers are inside the region. Interior nodes are nodes where
// all associated cell-centers are inside the region.

// Define a region that fills the whole simulation

#define everywhere (x == x) && (y == y) && (z == z)
//#define everywhere 1

// Define a macro to allow different parts of a region to be selected.
// Note: get_particle_bc_id returns  0 if NULL is passed to it
// Note: get_material_bc_id returns -1 if NULL is passed to it

#define leave_unchanged NULL

// FIXME: THESE GLOBAL POSITION CALCULATIONS NEED TO BE MADE MORE RIGOROUS

#define set_point_region_material( rgn, rmat ) do {                 \
    const material_id _rmat = get_material_id( (rmat) );            \
    if( _rmat==-1 ) break;                                          \
    const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;    \
    const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;    \
    const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;    \
    for( int _k=0; _k<_nz+2; _k++ ) { const double _zn = _z0 + _dz*(_k-1), _zc = _z0 + _dz*(_k-0.5); \
    for( int _j=0; _j<_ny+2; _j++ ) { const double _yn = _y0 + _dy*(_j-1), _yc = _y0 + _dy*(_j-0.5); field_t * _f = &field(0,_j,_k); \
    for( int _i=0; _i<_nx+2; _i++ ) { const double _xn = _x0 + _dx*(_i-1), _xc = _x0 + _dx*(_i-0.5); double x, y, z; \
          x = _xn; y = _yn; z = _zn; if( (rgn) ) _f->nmat  = _rmat; \
          x = _xc;                   if( (rgn) ) _f->ematx = _rmat; \
                   y = _yc;          if( (rgn) ) _f->fmatz = _rmat; \
                            z = _zc; if( (rgn) ) _f->cmat  = _rmat; \
                   y = _yn;          if( (rgn) ) _f->fmaty = _rmat; \
          x = _xn;                   if( (rgn) ) _f->ematz = _rmat; \
                   y = _yc;          if( (rgn) ) _f->fmatx = _rmat; \
                            z = _zn; if( (rgn) ) _f->ematy = _rmat; \
          _f++;                                                     \
    }}}                                                             \
  } while(0)

#define set_point_region_bc( rgn, ipbc, epbc ) do {            		 \
    const int64_t _ipbc = get_particle_bc_id( (particle_bc_t *)(ipbc) ); \
    const int64_t _epbc = get_particle_bc_id( (particle_bc_t *)(epbc) ); \
    if( !_ipbc && !_epbc ) break;                                        \
    const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;         \
    const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;         \
    const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;         \
    int64_t * _n0 = grid->neighbor;                                      \
    for( int _k=1; _k<_nz+1; _k++ ) { const double _zn = _z0 + _dz*(_k-1), _zh = _z0 + _dz*_k; \
    for( int _j=1; _j<_ny+1; _j++ ) { const double _yn = _y0 + _dy*(_j-1), _yh = _y0 + _dy*_j; int64_t * _n = _n0 + 6*voxel(1,_j,_k); \
    for( int _i=1; _i<_nx+1; _i++ ) { const double _xn = _x0 + _dx*(_i-1), _xh = _x0 + _dx*_i; double x, y, z; \
          int _r000, _r100, _r010, _r110, _r001, _r101, _r011, _r111;    \
          x = _xn; y = _yn; z = _zn; _r000 = (rgn);                      \
          x = _xh;                   _r100 = (rgn);                      \
          x = _xn; y = _yh;          _r010 = (rgn);                      \
          x = _xh;                   _r110 = (rgn);                      \
          x = _xn; y = _yn; z = _zh; _r001 = (rgn);                      \
          x = _xh;                   _r101 = (rgn);                      \
          x = _xn; y = _yh;          _r011 = (rgn);                      \
          x = _xh;                   _r111 = (rgn);                      \
          if( _epbc ) {                                                  \
            if( _r000 || _r010 || _r001 || _r011 ) _n[0] = _epbc;        \
            if( _r000 || _r001 || _r100 || _r101 ) _n[1] = _epbc;        \
            if( _r000 || _r100 || _r010 || _r110 ) _n[2] = _epbc;        \
            if( _r100 || _r110 || _r101 || _r111 ) _n[3] = _epbc;        \
            if( _r010 || _r011 || _r110 || _r111 ) _n[4] = _epbc;        \
            if( _r001 || _r101 || _r011 || _r111 ) _n[5] = _epbc;        \
          }                                                              \
          if( _ipbc ) {                                                  \
            if( _r000 && _r010 && _r001 && _r011 ) _n[0] = _ipbc         \
            if( _r000 && _r001 && _r100 && _r101 ) _n[1] = _ipbc         \
            if( _r000 && _r100 && _r010 && _r110 ) _n[2] = _ipbc         \
            if( _r100 && _r110 && _r101 && _r111 ) _n[3] = _ipbc         \
            if( _r010 && _r011 && _r110 && _r111 ) _n[4] = _ipbc         \
            if( _r001 && _r101 && _r011 && _r111 ) _n[5] = _ipbc         \
          }                                                              \
          _n += 6;                                                       \
    }}}                                                                  \
  } while(0)

// The equations are strictly evaluated inside the region
#define set_point_region_field( rgn,                                     \
                                eqn_ex, eqn_ey, eqn_ez,                  \
                                eqn_bx, eqn_by, eqn_bz ) do {            \
    const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;         \
    const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;         \
    const double _c  = grid->cvac;                                       \
    const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;         \
    for( int _k=0; _k<_nz+2; _k++ ) { const double _zn = _z0 + _dz*(_k-1), _zc = _z0 + _dz*(_k-0.5); \
    for( int _j=0; _j<_ny+2; _j++ ) { const double _yn = _y0 + _dy*(_j-1), _yc = _y0 + _dy*(_j-0.5); field_t * _f = &field(0,_j,_k); \
    for( int _i=0; _i<_nx+2; _i++ ) { const double _xn = _x0 + _dx*(_i-1), _xc = _x0 + _dx*(_i-0.5); double x, y, z; \
          x = _xn; y = _yn; z = _zn; /* No node fields */                \
          x = _xc;                   if( (rgn) ) _f->ex  =    (eqn_ex);  \
                   y = _yc;          if( (rgn) ) _f->cbz = _c*(eqn_bz);  \
                            z = _zc; /* No cell fields */                \
                   y = _yn;          if( (rgn) ) _f->cby = _c*(eqn_by);  \
          x = _xn;                   if( (rgn) ) _f->ez  =    (eqn_ez);  \
                   y = _yc;          if( (rgn) ) _f->cbx = _c*(eqn_bx);  \
                            z = _zn; if( (rgn) ) _f->ey  =    (eqn_ey);  \
          _f++;                                                          \
    }}}                                                                  \
  } while(0)

#define set_region_material( rgn, vmat, smat ) do {                    \
    const material_id _vmat = get_material_id( (vmat) );               \
    const material_id _smat = get_material_id( (smat) );               \
    if( _vmat==-1 && _smat==-1 ) break;                                \
    const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;       \
    const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;       \
    const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;       \
    for( int _k=0; _k<_nz+2; _k++ ) { const double _zl = _z0 + _dz*(_k-1.5), _zc = _z0 + _dz*(_k-0.5); \
    for( int _j=0; _j<_ny+2; _j++ ) { const double _yl = _y0 + _dy*(_j-1.5), _yc = _y0 + _dy*(_j-0.5); field_t *_f = &field(0,_j,_k); \
    for( int _i=0; _i<_nx+2; _i++ ) { const double _xl = _x0 + _dx*(_i-1.5), _xc = _x0 + _dx*(_i-0.5); double x, y, z; \
          int _rccc, _rlcc, _rclc, _rllc, _rccl, _rlcl, _rcll, _rlll;  \
          x = _xc; y = _yc; z = _zc; _rccc = (rgn);		       \
          x = _xl;                   _rlcc = (rgn);		       \
          x = _xc; y = _yl;          _rclc = (rgn);		       \
          x = _xl;                   _rllc = (rgn);		       \
          x = _xc; y = _yc; z = _zl; _rccl = (rgn);		       \
          x = _xl;                   _rlcl = (rgn);		       \
          x = _xc; y = _yl;          _rcll = (rgn);		       \
          x = _xl;                   _rlll = (rgn);		       \
          if( _smat!=-1 ) {                                            \
            if( _rccc || _rclc || _rccl || _rcll )  _f->ematx = _smat; \
            if( _rccc || _rccl || _rlcc || _rlcl )  _f->ematy = _smat; \
            if( _rccc || _rlcc || _rclc || _rllc )  _f->ematz = _smat; \
            if( _rccc || _rlcc )                    _f->fmatx = _smat; \
            if( _rccc || _rclc )                    _f->fmaty = _smat; \
            if( _rccc || _rccl )                    _f->fmatz = _smat; \
            if( _rccc || _rlcc || _rclc || _rllc ||		       \
                _rccl || _rlcl || _rcll || _rlll )  _f->nmat  = _smat; \
          }                                                            \
          if( _vmat!=-1 ) {                                            \
            if( _rccc && _rclc && _rccl && _rcll )  _f->ematx = _vmat; \
            if( _rccc && _rccl && _rlcc && _rlcl )  _f->ematy = _vmat; \
            if( _rccc && _rlcc && _rclc && _rllc )  _f->ematz = _vmat; \
            if( _rccc && _rlcc )                    _f->fmatx = _vmat; \
            if( _rccc && _rclc )                    _f->fmaty = _vmat; \
            if( _rccc && _rccl )                    _f->fmatz = _vmat; \
            if( _rccc && _rlcc && _rclc && _rllc &&		       \
                _rccl && _rlcl && _rcll && _rlll )  _f->nmat  = _vmat; \
            if( _rccc )                             _f->cmat  = _vmat; \
          }							       \
          _f++;                                                        \
    }}}                                                                \
  } while(0)

#define set_region_bc( rgn, vpbc, ipbc, epbc ) do {                      \
    const int64_t _vpbc = get_particle_bc_id( (particle_bc_t *)(vpbc) ); \
    const int64_t _ipbc = get_particle_bc_id( (particle_bc_t *)(ipbc) ); \
    const int64_t _epbc = get_particle_bc_id( (particle_bc_t *)(epbc) ); \
    if( !_vpbc && !_ipbc && !_epbc ) break;                              \
    const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;         \
    const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;         \
    const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;         \
    int64_t * _n0 = grid->neighbor;                                      \
    for( int _k=1; _k<_nz+1; _k++ ) { const double _zl = _z0 + _dz*(_k-1.5), _zc = _z0 + _dz*(_k-0.5), _zh = _z0 + _dz*(_k+0.5); \
    for( int _j=1; _j<_ny+1; _j++ ) { const double _yl = _y0 + _dy*(_j-1.5), _yc = _y0 + _dy*(_j-0.5), _yh = _y0 + _dy*(_j+0.5); int64_t * _n = _n0 + 6*voxel(1,_j,_k); \
    for( int _i=1; _i<_nx+1; _i++ ) { const double _xl = _x0 + _dx*(_i-1.5), _xc = _x0 + _dx*(_i-0.5), _xh = _x0 + _dx*(_i+0.5); double x, y, z; \
          int _rc, _r0, _r1, _r2, _r3, _r4, _r5;                         \
          x = _xc; y = _yc; z = _zc; _rc = (rgn);                        \
          x = _xl; y = _yc; z = _zc; _r0 = (rgn);                        \
          x = _xc; y = _yl; z = _zc; _r1 = (rgn);                        \
          x = _xc; y = _yc; z = _zl; _r2 = (rgn);                        \
          x = _xh; y = _yc; z = _zc; _r3 = (rgn);                        \
          x = _xc; y = _yh; z = _zc; _r4 = (rgn);                        \
          x = _xc; y = _yc; z = _zh; _r5 = (rgn);                        \
          if( _vpbc ) {                                                  \
            if( _rc && _r0  ) _n[0] = _vpbc;                             \
            if( _rc && _r1  ) _n[1] = _vpbc;                             \
            if( _rc && _r2  ) _n[2] = _vpbc;                             \
            if( _rc && _r3  ) _n[3] = _vpbc;                             \
            if( _rc && _r4  ) _n[4] = _vpbc;                             \
            if( _rc && _r5  ) _n[5] = _vpbc;                             \
          }                                                              \
          if( _ipbc ) {                                                  \
            if( _rc && !_r0 ) _n[0] = _ipbc;                             \
            if( _rc && !_r1 ) _n[1] = _ipbc;                             \
            if( _rc && !_r2 ) _n[2] = _ipbc;                             \
            if( _rc && !_r3 ) _n[3] = _ipbc;                             \
            if( _rc && !_r4 ) _n[4] = _ipbc;                             \
            if( _rc && !_r5 ) _n[5] = _ipbc;                             \
          }                                                              \
          if( _epbc ) {                                                  \
            if( !_rc && _r0 ) _n[0] = _epbc;                             \
            if( !_rc && _r1 ) _n[1] = _epbc;                             \
            if( !_rc && _r2 ) _n[2] = _epbc;                             \
            if( !_rc && _r3 ) _n[3] = _epbc;                             \
            if( !_rc && _r4 ) _n[4] = _epbc;                             \
            if( !_rc && _r5 ) _n[5] = _epbc;                             \
          }                                                              \
          _n += 6;                                                       \
    }}}                                                                  \
  } while(0)

// rgn is a logical equation that specifies the interior of the volume
// emitter.  This mechanism is only efficient for volumeteric emission
// processes that occupy a small portion of the simulation volume.
// For volumetric emission processes that occupy the entire simulation
// volume, recommend using the begin_particle_injection { }; input
// deck segment.

#define define_volume_emitter( e, rgn ) do {                      \
    /* Count the number of cells in the emitter */                \
    const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;  \
    const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;  \
    const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;  \
    int _nc = 0;                                                  \
    for( int _k=1; _k<_nz+1; _k++ ) { const double _zc = _z0 + _dz*(_k-0.5); \
    for( int _j=1; _j<_ny+1; _j++ ) { const double _yc = _y0 + _dy*(_j-0.5); \
    for( int _i=1; _i<_nx+1; _i++ ) { const double _xc = _x0 + _dx*(_i-0.5); double x, y, z; \
          x = _xc; y = _yc; z = _zc; if( (rgn) ) _nc++;           \
    }}}                                                           \
    /* Define the emitter */                                      \
    int32_t * _c = size_emitter( define_emitter( (e) ), _nc );    \
    _nc = 0;                                                      \
    for( int _k=1; _k<_nz+1; _k++ ) { const double _zc = _z0 + _dz*(_k-0.5); \
    for( int _j=1; _j<_ny+1; _j++ ) { const double _yc = _y0 + _dy*(_j-0.5); \
    for( int _i=1; _i<_nx+1; _i++ ) { const double _xc = _x0 + _dx*(_i-0.5); double x, y, z; \
          x = _xc; y = _yc; z = _zc;                              \
          if( (rgn) ) _c[_nc++] =                                 \
            COMPONENT_ID( voxel(_i,_j,_k), BOUNDARY(0,0,0) );     \
    }}}                                                           \
  } while(0)

// rgn is a logical equation.
// rgn = true for interior of region
// rgn = false for exterior of region
// A surface emitter emits into the exterior of the region.

#define define_surface_emitter( e, rgn ) do {                    \
    /* Count the number of faces in emitter surface */           \
    const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0; \
    const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz; \
    const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz; \
    int _nf = 0;                                                 \
    for( int _k=1; _k<_nz+1; _k++ ) { const double _zl = _z0 + _dz*(_k-1.5), _zc = _z0 + _dz*(_k-0.5), _zh = _z0 + _dz*(_k+0.5); \
    for( int _j=1; _j<_ny+1; _j++ ) { const double _yl = _y0 + _dy*(_j-1.5), _yc = _y0 + _dy*(_j-0.5), _yh = _y0 + _dy*(_j+0.5); \
    for( int _i=1; _i<_nx+1; _i++ ) { const double _xl = _x0 + _dx*(_i-1.5), _xc = _x0 + _dx*(_i-0.5), _xh = _x0 + _dx*(_i+0.5); double x, y, z; \
          int _rc, _r0, _r1, _r2, _r3, _r4, _r5;                 \
          x = _xc; y = _yc; z = _zc; _rc = (rgn);                \
          x = _xl; y = _yc; z = _zc; _r0 = (rgn);                \
          x = _xc; y = _yl; z = _zc; _r1 = (rgn);                \
          x = _xc; y = _yc; z = _zl; _r2 = (rgn);                \
          x = _xh; y = _yc; z = _zc; _r3 = (rgn);                \
          x = _xc; y = _yh; z = _zc; _r4 = (rgn);                \
          x = _xc; y = _yc; z = _zh; _r5 = (rgn);                \
          if( !_rc && _r0 ) _nf++;                               \
          if( !_rc && _r1 ) _nf++;                               \
          if( !_rc && _r2 ) _nf++;                               \
          if( !_rc && _r3 ) _nf++;                               \
          if( !_rc && _r4 ) _nf++;                               \
          if( !_rc && _r5 ) _nf++;                               \
    }}}                                                          \
    /* Define the emitter */                                     \
    int32_t * _c = size_emitter( define_emitter( (e) ), _nf );   \
    _nf = 0;                                                     \
    for( int _k=1; _k<_nz+1; _k++ ) { const double _zl = _z0 + _dz*(_k-1.5), _zc = _z0 + _dz*(_k-0.5), _zh = _z0 + _dz*(_k+0.5); \
    for( int _j=1; _j<_ny+1; _j++ ) { const double _yl = _y0 + _dy*(_j-1.5), _yc = _y0 + _dy*(_j-0.5), _yh = _y0 + _dy*(_j+0.5); \
    for( int _i=1; _i<_nx+1; _i++ ) { const double _xl = _x0 + _dx*(_i-1.5), _xc = _x0 + _dx*(_i-0.5), _xh = _x0 + _dx*(_i+0.5); double x, y, z; \
          int _rc, _r0, _r1, _r2, _r3, _r4, _r5;                 \
          x = _xc; y = _yc; z = _zc; _rc = (rgn);                \
          x = _xl; y = _yc; z = _zc; _r0 = (rgn);                \
          x = _xc; y = _yl; z = _zc; _r1 = (rgn);                \
          x = _xc; y = _yc; z = _zl; _r2 = (rgn);                \
          x = _xh; y = _yc; z = _zc; _r3 = (rgn);                \
          x = _xc; y = _yh; z = _zc; _r4 = (rgn);                \
          x = _xc; y = _yc; z = _zh; _r5 = (rgn);                \
          if( !_rc && _r0 ) _c[_nf++] = COMPONENT_ID( voxel(_i,_j,_k), BOUNDARY(-1, 0, 0) ); \
          if( !_rc && _r1 ) _c[_nf++] = COMPONENT_ID( voxel(_i,_j,_k), BOUNDARY( 0,-1, 0) ); \
          if( !_rc && _r2 ) _c[_nf++] = COMPONENT_ID( voxel(_i,_j,_k), BOUNDARY( 0, 0,-1) ); \
          if( !_rc && _r3 ) _c[_nf++] = COMPONENT_ID( voxel(_i,_j,_k), BOUNDARY( 1, 0, 0) ); \
          if( !_rc && _r4 ) _c[_nf++] = COMPONENT_ID( voxel(_i,_j,_k), BOUNDARY( 0, 1, 0) ); \
          if( !_rc && _r5 ) _c[_nf++] = COMPONENT_ID( voxel(_i,_j,_k), BOUNDARY( 0, 0, 1) ); \
    }}}                                                          \
  } while(0)

// The equations are only evaluated inside the mesh-mapped region
// (This is not strictly inside the region)
#define set_region_field( rgn,                                        \
                          eqn_ex, eqn_ey, eqn_ez,                     \
                          eqn_bx, eqn_by, eqn_bz ) do {               \
    const double _x0 = grid->x0, _y0 = grid->y0, _z0 = grid->z0;      \
    const double _dx = grid->dx, _dy = grid->dy, _dz = grid->dz;      \
    const double _c  = grid->cvac;                                    \
    const int    _nx = grid->nx, _ny = grid->ny, _nz = grid->nz;      \
    for( int _k=0; _k<_nz+2; _k++ ) { const double _zl = _z0 + _dz*(_k-1.5), _ze = _z0 + _dz*_k, _zc = _z0 + _dz*(_k-0.5); \
    for( int _j=0; _j<_ny+2; _j++ ) { const double _yl = _y0 + _dy*(_j-1.5), _ye = _y0 + _dy*_j, _yc = _y0 + _dy*(_j-0.5); field_t *_f = &field(0,_j,_k); \
    for( int _i=0; _i<_nx+2; _i++ ) { const double _xl = _x0 + _dx*(_i-1.5), _xe = _x0 + _dx*_i, _xc = _x0 + _dx*(_i-0.5); double x, y, z; \
          int _rccc, _rlcc, _rclc, _rllc, _rccl, _rlcl, _rcll;        \
          x = _xc; y = _yc; z = _zc; _rccc = (rgn);                   \
          x = _xl;                   _rlcc = (rgn);                   \
          x = _xc; y = _yl;          _rclc = (rgn);                   \
          x = _xl;                   _rllc = (rgn);                   \
          x = _xc; y = _yc; z = _zl; _rccl = (rgn);                   \
          x = _xl;                   _rlcl = (rgn);                   \
          x = _xc; y = _yl;          _rcll = (rgn);                   \
          x = _xc; y = _ye; z = _ze; if( _rccc || _rclc || _rccl || _rcll ) _f->ex  =    (eqn_ex); \
          x = _xe; y = _yc; z = _ze; if( _rccc || _rccl || _rlcc || _rlcl ) _f->ey  =    (eqn_ey); \
          x = _xe; y = _ye; z = _zc; if( _rccc || _rlcc || _rclc || _rllc ) _f->ez  =    (eqn_ez); \
          x = _xe; y = _yc; z = _zc; if( _rccc || _rlcc )                   _f->cbx = _c*(eqn_bx); \
          x = _xc; y = _ye; z = _zc; if( _rccc || _rclc )                   _f->cby = _c*(eqn_by); \
          x = _xc; y = _yc; z = _ze; if( _rccc || _rccl )                   _f->cbz = _c*(eqn_bz); \
          _f++;                                                       \
    }}}                                                               \
  } while(0)

// In main.cxx

void
checkpt( const char * fbase,
         int tag );

//-----------------------------------------------------------------------------

// Include the users input deck
#include EXPAND_AND_STRINGIFY(INPUT_DECK)
