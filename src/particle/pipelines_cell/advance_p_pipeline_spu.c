#ifdef CELL_SPU_BUILD

#define IN_particle_pipeline
#include <particle_pipelines.h>
#include <v4c_spu.h>
#include <spu_mfcio.h>

#ifdef IN_HARNESS
#include <profile.h>
#endif

// DMA tag usage:
//  0: 2 - Particle buffer 0 (read, write, mover write)
//  3: 5 - Particle buffer 1 (read, write, mover write)
//  6: 8 - Particle buffer 2 (read, write, mover write)
//  9:16 - Interpolator cache
// 17:25 - Accumulator cache
// 26:29 - Neighbor cache
// 30:30 - Pipeline input arguments
// 31:31 - Pipeline return values

///////////////////////////////////////////////////////////////////////////////
// Pipeline Input / Output arguments

DECLARE_ALIGNED_ARRAY( advance_p_pipeline_args_t, 128, args, 1 );
DECLARE_ALIGNED_ARRAY( particle_mover_seg_t,      128, seg,  1 );
int pipeline_rank;
int n_pipeline;

///////////////////////////////////////////////////////////////////////////////
// Pipeline buffers

// NP_BLOCK_TARGET is such that a 16Kb of particle data is located in
// a particle block (the maximum amount of data that can be moved in a
// single DMA transfer).  Triple buffering requires 3 such particle
// blocks.  So local particle storage takes up 48Kb.  Each local
// particle has a mover associated with it such that a SPU pipeline
// cannot run out movers during particle processing.  This costs 24Kb.

#define NP_BLOCK_TARGET 512
DECLARE_ALIGNED_ARRAY( particle_t,       128, local_p,  3*NP_BLOCK_TARGET );
DECLARE_ALIGNED_ARRAY( particle_mover_t, 128, local_pm, 3*NP_BLOCK_TARGET );

///////////////////////////////////////////////////////////////////////////////
// External memory caches

// #define CACHE_STATS

// Since DMA transfers seem optimized for 128-bytes at aligned
// addresses, set up the cache lines to use 128-byte aligned lines
// 128-bytes in size.

// Interpolator cache: 512 cached interpolators (64Kb). Roughly four
// second nearest neighborhoods (125 voxels) around a particle can exist
// within the cache.

#undef CACHE_NAME
#undef CACHED_TYPE
#undef CACHE_TYPE
#undef CACHELINE_LOG2SIZE
#undef CACHE_LOG2NWAY
#undef CACHE_LOG2NSETS
#undef CACHE_SET_TAGID
#undef CACHE_READ_X4

#define CACHE_NAME           interpolator_cache
#define CACHED_TYPE          interpolator_t
#define CACHE_TYPE           0             /* r/o */
#define CACHELINE_LOG2SIZE   7             /* 1 per line - 128 byte lines */
#define CACHE_LOG2NWAY       2             /* 4 way */
#define CACHE_LOG2NSETS      7             /* 128 lines per way */
#define CACHE_SET_TAGID(set) (9+(set)&0x7) /* tags 9:16 */
#include <cache-api.h>

#define PTR_INTERPOLATOR(v) cache_rw( interpolator_cache,               \
                                      args->f0 +                        \
                                      (v)*sizeof(interpolator_t) )

// Accumulator cache: 512 cached accumulators (32Kb).  Roughly four
// second nearest neighborhood of accumulators can exist within the
// cache.

#undef CACHE_NAME
#undef CACHED_TYPE
#undef CACHE_TYPE
#undef CACHELINE_LOG2SIZE
#undef CACHE_LOG2NWAY
#undef CACHE_LOG2NSETS
#undef CACHE_SET_TAGID
#undef CACHE_READ_X4

#define CACHE_NAME           accumulator_cache
#define CACHED_TYPE          accumulator_t
#define CACHE_TYPE           1              /* r/w */
#define CACHELINE_LOG2SIZE   7              /* 2 per line - 128 byte lines */
#define CACHE_LOG2NWAY       2              /* 4 way */
#define CACHE_LOG2NSETS      6              /* 64 lines per way */
#define CACHE_SET_TAGID(set) (17+(set)&0x7) /* tags 17:25 */
#include <cache-api.h>

#define PTR_ACCUMULATOR(v) cache_rw( accumulator_cache,                 \
                                     args->a0 +                         \
                                     (v)*sizeof(accumulator_t) )

// Neighbor cache: 2048 cached cell adjacencies (16K).  Roughly three
// second nearest neighborhoods of voxel adjacencies can exist within the
// cache.

#undef CACHE_NAME
#undef CACHED_TYPE
#undef CACHE_TYPE
#undef CACHELINE_LOG2SIZE
#undef CACHE_LOG2NWAY
#undef CACHE_LOG2NSETS
#undef CACHE_SET_TAGID
#undef CACHE_READ_X4

#define CACHE_NAME           neighbor_cache
#define CACHED_TYPE          int64_t
#define CACHE_TYPE           0              /* r/o */
#define CACHELINE_LOG2SIZE   7              /* 16 per line - 128 byte lines */
#define CACHE_LOG2NWAY       2              /* 4 way */
#define CACHE_LOG2NSETS      5              /* 32 lines per way */
#define CACHE_SET_TAGID(set) (25+(set)&0x3) /* tags 25:29 */
#include <cache-api.h>

#define NEIGHBOR(v,face) cache_rd( neighbor_cache,                      \
                                   args->neighbor +                     \
                                   6*sizeof(int64_t)*(v) +              \
                                   (face)*sizeof(int64_t) )

///////////////////////////////////////////////////////////////////////////////
// Computational kernel

// move_p moves the particle m->p by m->dispx, m->dispy, m->dispz
// depositing particle current as it goes. If the particle was moved
// sucessfully (particle mover is no longer in use), this returns 0.
// If the particle interacted with something this routine could not
// handle, this routine returns 1 (particle mover is still in use).
// On a successful move, the particle position is updated and dispx,
// dispy and dispz are zerod.  On a partial move, the particle
// position is updated to the point where the particle interacted and
// dispx, dispy, dispz contain the remaining particle displacement.
// The displacements are the physical displacments normalized current
// cell size.
//
// Because move_p is internal use only and frequently called, it does
// not check its input arguments.  Callers are responsible for
// insuring valid arguments.

static int
move_p_spu( particle_t       * ALIGNED(32) p,
            particle_mover_t * ALIGNED(16) pm ) {
  float s_midx, s_midy, s_midz;
  float s_dispx, s_dispy, s_dispz;
  float s_dir[3];
  float v0, v1, v2, v3, v4, v5;
  int type;
  int64_t neighbor;
  float * ALIGNED(16) a;

  // FIXME: THIS CAN BE PARTIALLY HORIZONTAL SIMD ACCELERATED

  for(;;) {
    s_midx = p->dx;
    s_midy = p->dy;
    s_midz = p->dz;

    s_dispx = pm->dispx;
    s_dispy = pm->dispy;
    s_dispz = pm->dispz;

    s_dir[0] = (s_dispx>0) ? 1 : -1;
    s_dir[1] = (s_dispy>0) ? 1 : -1;
    s_dir[2] = (s_dispz>0) ? 1 : -1;
    
    // Compute the twice the fractional distance to each potential
    // streak/cell face intersection.

    v0 = (s_dispx==0) ? 3.4e38 : (s_dir[0]-s_midx)/s_dispx;
    v1 = (s_dispy==0) ? 3.4e38 : (s_dir[1]-s_midy)/s_dispy;
    v2 = (s_dispz==0) ? 3.4e38 : (s_dir[2]-s_midz)/s_dispz;

    // Determine the fractional length and type of current streak. The
    // streak ends on either the first face intersected by the
    // particle track or at the end of the particle track.
    // 
    //   type 0,1 or 2 ... streak ends on a x,y or z-face respectively
    //   type 3        ... streak ends at end of the particle track

    /**/      v3=2,  type=3;
    if(v0<v3) v3=v0, type=0;
    if(v1<v3) v3=v1, type=1;
    if(v2<v3) v3=v2, type=2;
    v3 *= 0.5;

    // Compute the midpoint and the normalized displacement of the streak

    s_dispx *= v3;
    s_dispy *= v3;
    s_dispz *= v3;

    s_midx += s_dispx;
    s_midy += s_dispy;
    s_midz += s_dispz;

    // Accumulate the streak.  Note: accumulator values are 4 times
    // the total physical charge that passed through the appropriate
    // current quadrant in a time-step

    v5 = p->q*s_dispx*s_dispy*s_dispz*(1./3.);
    a = (float * ALIGNED(64))PTR_ACCUMULATOR(p->i); // No need to lock

#   define ACCUMULATE_J(X,Y,Z,offset)                                   \
    v4  = p->q*s_disp##X; /* v2 = q ux                            */    \
    v1  = v4*s_mid##Y;    /* v1 = q ux dy                         */    \
    v0  = v4-v1;          /* v0 = q ux (1-dy)                     */    \
    v1 += v4;             /* v1 = q ux (1+dy)                     */    \
    v4  = 1+s_mid##Z;     /* v4 = 1+dz                            */    \
    v2  = v0*v4;          /* v2 = q ux (1-dy)(1+dz)               */    \
    v3  = v1*v4;          /* v3 = q ux (1+dy)(1+dz)               */    \
    v4  = 1-s_mid##Z;     /* v4 = 1-dz                            */    \
    v0 *= v4;             /* v0 = q ux (1-dy)(1-dz)               */    \
    v1 *= v4;             /* v1 = q ux (1+dy)(1-dz)               */    \
    v0 += v5;             /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */    \
    v1 -= v5;             /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */    \
    v2 -= v5;             /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */    \
    v3 += v5;             /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */    \
    a[offset+0] += v0;                                                  \
    a[offset+1] += v1;                                                  \
    a[offset+2] += v2;                                                  \
    a[offset+3] += v3

    ACCUMULATE_J( x,y,z, 0 );
    ACCUMULATE_J( y,z,x, 4 );
    ACCUMULATE_J( z,x,y, 8 );

#   undef ACCUMULATE_J

    // Compute the remaining particle displacment

    pm->dispx -= s_dispx;
    pm->dispy -= s_dispy;
    pm->dispz -= s_dispz;

    // Compute the new particle offset

    p->dx += s_dispx+s_dispx;
    p->dy += s_dispy+s_dispy;
    p->dz += s_dispz+s_dispz;

    // If an end streak, return success (should be ~50% of the time)

    if( type==3 ) return 0;

    // Determine if the cell crossed into a local cell or if it hit a
    // boundary.  Convert the coordinate system accordingly.  Note:
    // Crossing into a local cell should happen the other ~50% of
    // time; hitting a structure and parallel domain boundary should
    // usually be a rare event.  Note: the entry / exit coordinate for
    // the particle is guaranteed to be +/-1 _exactly_ for the
    // particle.

    // FIXME: WOULD SWITCHING ON TYPE BE FASTER ON THE SPUS?

    v0 = s_dir[type];
    neighbor = NEIGHBOR( p->i, ((v0>0)?3:0) + type );
    if( neighbor<args->rangel || neighbor>args->rangeh ) { // Hit a boundary
      (&(p->dx))[type] = v0;                               // Put on boundary
      if( neighbor!=reflect_particles ) return 1;          // Cannot handle it
      (&(p->ux))[type] = -(&(p->ux))[type];
      (&(pm->dispx))[type] = -(&(pm->dispx))[type];
    } else {
      p->i = neighbor - args->rangel; // Compute local index of neighbor
      /**/                            // Note: neighbor-args->rangel < 2^31 / 6
      (&(p->dx))[type] = -v0;         // Convert coordinate system
    }
  }
  return 0; // Never get here ... avoid compiler warning
}

// FIXME: Using restricted pointers makes this worse(!) on gcc??
// FIXME: Branch hints makes this worse(!) on gcc?? (No change on xlc.)

static int                             // Return number of movers used
advance_p_pipeline_spu( particle_t       * ALIGNED(128) p, // Particle array
                        particle_mover_t * ALIGNED(16)  pm, // Mover array
                        int idx,       // Index of first particle
                        int nq ) {     // Number of particle quads
  USING_V4C;

  const vec_float4 qdt_2mc        = VEC_FLOAT4( args->qdt_2mc );
  const vec_float4 cdt_dx         = VEC_FLOAT4( args->cdt_dx  );
  const vec_float4 cdt_dy         = VEC_FLOAT4( args->cdt_dy  );
  const vec_float4 cdt_dz         = VEC_FLOAT4( args->cdt_dz  );
  const vec_float4 one            = VEC_FLOAT4(  1.           );
  const vec_float4 one_third      = VEC_FLOAT4(  1./3.        );
  const vec_float4 one_half       = VEC_FLOAT4(  1./2.        );
  const vec_float4 two_fifteenths = VEC_FLOAT4(  2./15.       );
  const vec_float4 neg_one        = VEC_FLOAT4( -1.           );

  vec_float4 p0r, p1r, p2r, p3r;       // transposed position
  vec_float4 p0u, p1u, p2u, p3u;       // transposed momentum

  vec_float4 dx, dy, dz; vec_int4 i;   // particle position
  vec_float4 ux, uy, uz, q;            // particle momentum and charge

  vec_float4 ex0,  dexdy,  dexdz, d2exdydz; // ex interp coeff
  vec_float4 ey0,  deydz,  deydx, d2eydzdx; // ey interp coeff
  vec_float4 ez0,  dezdx,  dezdy, d2ezdxdy; // ez interp coeff
  vec_float4 cbx0, dcbxdx, cby0,  dcbydy; // bx and by interp coeff
  vec_float4 cbz0, dcbzdz, v12,   v13; // bz interp coeff and scratch

  vec_float4 hax, hay, haz;            // half particle electric field accel
  vec_float4 cbx, cby, cbz;            // particle magnetic field

  vec_float4 ux0, uy0, uz0;            // pre Boris rotation momentum
  vec_float4 v14, cbs, ths, v15, v16, v17; // Boris rotation scalars
  vec_float4 wx0, wy0, wz0;            // intermediate Boris rotation momentum
  vec_float4 uxh, uyh, uzh;            // new particle momentum

  vec_float4 rgamma;                   // new inverse gamma
  vec_float4 ddx, ddy, ddz;            // cell-normalized particle displacment
  vec_float4 dxh, dyh, dzh;            // half-step position rel to start cell
  vec_float4 dx1, dy1, dz1;            // new position rel to start cell
  vec_uint4  outbnd;                   // Boolean (true if particle left cell)

  vec_float4 qa, ccc;
  vec_float4 a0x, a1x, a2x, a3x, a4x;
  vec_float4 a0y, a1y, a2y, a3y, a4y;
  vec_float4 a0z, a1z, a2z, a3z, a4z;

  int i0, i1, i2, i3;

  const interpolator_t * ALIGNED(128) fi;
  accumulator_t * ALIGNED(64) ja;
  int nm = 0;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=4, idx+=4 ) {

    // Load the particle quad positions

    LOAD_4x1( &p[0].dx, p0r );
    LOAD_4x1( &p[1].dx, p1r );
    LOAD_4x1( &p[2].dx, p2r );
    LOAD_4x1( &p[3].dx, p3r );
    TRANSPOSE( p0r, p1r, p2r, p3r );
    dx = p0r;
    dy = p1r;
    dz = p2r;
    i  = (vec_int4)p3r;

    // Interpolate fields

    i0 = EXTRACT( i, 0 );
    fi = PTR_INTERPOLATOR(i0);
    LOAD_4x1( &fi->ex,  ex0  );
    LOAD_4x1( &fi->ey,  ey0  );
    LOAD_4x1( &fi->ez,  ez0  );
    LOAD_4x1( &fi->cbx, cbx0 );
    LOAD_4x1( &fi->cbz, cbz0 );

    i1 = EXTRACT( i, 1 );
    fi = PTR_INTERPOLATOR(i1);
    LOAD_4x1( &fi->ex,  dexdy  );
    LOAD_4x1( &fi->ey,  deydz  );
    LOAD_4x1( &fi->ez,  dezdx  );
    LOAD_4x1( &fi->cbx, dcbxdx );
    LOAD_4x1( &fi->cbz, dcbzdz );

    i2 = EXTRACT( i, 2 );
    fi = PTR_INTERPOLATOR(i2);
    LOAD_4x1( &fi->ex,  dexdz );
    LOAD_4x1( &fi->ey,  deydx );
    LOAD_4x1( &fi->ez,  dezdy );
    LOAD_4x1( &fi->cbx, cby0  );
    LOAD_4x1( &fi->cbz, v12   );

    i3 = EXTRACT( i, 3 );
    fi = PTR_INTERPOLATOR(i3);
    LOAD_4x1( &fi->ex,  d2exdydz );
    LOAD_4x1( &fi->ey,  d2eydzdx );
    LOAD_4x1( &fi->ez,  d2ezdxdy );
    LOAD_4x1( &fi->cbx, dcbydy   );
    LOAD_4x1( &fi->cbz, v13      );

    TRANSPOSE( ex0, dexdy, dexdz, d2exdydz );
    hax = MUL( qdt_2mc, FMA( FMA( d2exdydz, dy, dexdz ), dz,
                             FMA( dexdy,    dy, ex0   ) ) );

    TRANSPOSE( ey0, deydz, deydx, d2eydzdx );
    hay = MUL( qdt_2mc, FMA( FMA( d2eydzdx, dz, deydx ), dx,
                             FMA( deydz,    dz, ey0 ) ) );

    TRANSPOSE( ez0, dezdx, dezdy, d2ezdxdy );
    haz = MUL( qdt_2mc, FMA( FMA( d2ezdxdy, dx, dezdy ), dy,
                             FMA( dezdx,    dx, ez0 ) ) );

    TRANSPOSE( cbx0, dcbxdx, cby0, dcbydy );
    cbx = FMA( dcbxdx, dx, cbx0 );
    cby = FMA( dcbydy, dy, cby0 );

    HALF_TRANSPOSE( cbz0, dcbzdz, v12, v13 );
    cbz = FMA( dcbzdz, dz, cbz0 );

    // Update momentum.  Note: Could eliminate a dependency in v14 calc
    // if willing to play fast and loose with numerics (saves about a spu
    // clock per particle).

    LOAD_4x1( &p[0].ux, p0u );
    LOAD_4x1( &p[1].ux, p1u );
    LOAD_4x1( &p[2].ux, p2u );
    LOAD_4x1( &p[3].ux, p3u );
    TRANSPOSE( p0u, p1u, p2u, p3u );
    ux  = p0u;
    uy  = p1u;
    uz  = p2u;
    q   = p3u;
    ux0 = ADD( ux, hax );
    uy0 = ADD( uy, hay );
    uz0 = ADD( uz, haz );
    v14 = MUL( qdt_2mc, RSQRT( ADD( one, FMA( ux0,ux0, FMA( uy0,uy0, MUL(uz0,uz0) ) ) ) ) );
    cbs = FMA( cbx,cbx, FMA( cby,cby, MUL(cbz,cbz) ) ); // |cB|^2
    ths = MUL( MUL( v14,v14 ), cbs );  // |theta|^2 = |wc dt/2|^2
    v15 = MUL( v14, FMA( FMA( two_fifteenths, ths, one_third ), ths, one ) );
    v16 = MUL( v15, RCP( FMA( MUL(v15,v15), cbs, one ) ) );
    v17 = ADD( v16, v16 );
    wx0 =      FMA( FMS( uy0,cbz, MUL(uz0,cby) ), v15, ux0 );
    wy0 =      FMA( FMS( uz0,cbx, MUL(ux0,cbz) ), v15, uy0 );
    wz0 =      FMA( FMS( ux0,cby, MUL(uy0,cbx) ), v15, uz0 );
    uxh = ADD( FMA( FMS( wy0,cbz, MUL(wz0,cby) ), v17, ux0 ), hax );
    uyh = ADD( FMA( FMS( wz0,cbx, MUL(wx0,cbz) ), v17, uy0 ), hay );
    uzh = ADD( FMA( FMS( wx0,cby, MUL(wy0,cbx) ), v17, uz0 ), haz );
    p0u = uxh;
    p1u = uyh;
    p2u = uzh;
    /* p3u is unchanged */
    TRANSPOSE( p0u, p1u, p2u, p3u );
    STORE_4x1( p0u, &p[0].ux );
    STORE_4x1( p1u, &p[1].ux );
    STORE_4x1( p2u, &p[2].ux );
    STORE_4x1( p3u, &p[3].ux );
    
    // Update the position of inbnd particles

    rgamma = RSQRT( ADD( one, FMA( uxh,uxh, FMA( uyh,uyh, MUL(uzh,uzh) ) ) ) );
    ddx    = MUL( MUL( uxh, cdt_dx ), rgamma );
    ddy    = MUL( MUL( uyh, cdt_dy ), rgamma );
    ddz    = MUL( MUL( uzh, cdt_dz ), rgamma ); // voxel normalized movement
    dxh    = ADD( dx,  ddx );
    dyh    = ADD( dy,  ddy );
    dzh    = ADD( dz,  ddz );          // Half step position
    dx1    = ADD( dxh, ddx );
    dy1    = ADD( dyh, ddy );
    dz1    = ADD( dzh, ddz );          // New particle position
    outbnd = OR( OR( OR( CMPLT(dx1,neg_one), CMPGT(dx1,one) ),
                     OR( CMPLT(dy1,neg_one), CMPGT(dy1,one) ) ),
                     OR( CMPLT(dz1,neg_one), CMPGT(dz1,one) ) );
    p0r    = MERGE( outbnd, dx, dx1 ); // Do not update outbnd particles
    p1r    = MERGE( outbnd, dy, dy1 );
    p2r    = MERGE( outbnd, dz, dz1 );
    /* p3r is unchanged */
    TRANSPOSE( p0r, p1r, p2r, p3r );
    STORE_4x1( p0r, &p[0].dx );
    STORE_4x1( p1r, &p[1].dx );
    STORE_4x1( p2r, &p[2].dx );
    STORE_4x1( p3r, &p[3].dx );
   
    // Accumulate current of inbnd particles
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step

    qa  = CZERO( outbnd, q );          // Do not accumulate outbnd particles
    /**/                               // The half step position is not the
    /**/                               // streak midpoint for them!
    ccc = MUL( MUL( one_third, qa ), MUL( ddx, MUL( ddy, ddz ) ) );
    /**/                               // Charge conservation correction

#   define ACCUMULATE_J(X,Y,Z)                                               \
    a4##X = MUL(qa,   dd##X);   /* a4 = qa ddx                            */ \
    a1##X = MUL(a4##X,d##Y##h); /* a1 = qa ddx dyh                        */ \
    a0##X = SUB(a4##X,a1##X);   /* a0 = qa ddx (1-dyh)                    */ \
    a1##X = ADD(a1##X,a4##X);   /* a1 = qa ddx (1+dyh)                    */ \
    a4##X = ADD(one,  d##Z##h); /* a4 = 1+dzh                             */ \
    a2##X = MUL(a0##X,a4##X);   /* a2 = qa ddx (1-dyh)(1+dzh)             */ \
    a3##X = MUL(a1##X,a4##X);   /* a3 = qa ddx (1+dyh)(1+dzh)             */ \
    a4##X = SUB(one,  d##Z##h); /* a4 = 1-dzh                             */ \
    a0##X = MUL(a0##X,a4##X);   /* a0 = qa ddx (1-dyh)(1-dzh)             */ \
    a1##X = MUL(a1##X,a4##X);   /* a1 = qa ddx (1+dyh)(1-dzh)             */ \
    a0##X = ADD(a0##X,ccc);     /* a0 = qa ddx [(1-dyh)(1-dzh)+ddy*ddz/3] */ \
    a1##X = SUB(a1##X,ccc);     /* a1 = qa ddx [(1+dyh)(1-dzh)-ddy*ddz/3] */ \
    a2##X = SUB(a2##X,ccc);     /* a2 = qa ddx [(1-dyh)(1+dzh)-ddy*ddz/3] */ \
    a3##X = ADD(a3##X,ccc);     /* a3 = qa ddx [(1+dyh)(1+dzh)+ddy*ddz/3] */ \
    TRANSPOSE(a0##X,a1##X,a2##X,a3##X)

    ACCUMULATE_J( x,y,z );
    ACCUMULATE_J( y,z,x );
    ACCUMULATE_J( z,x,y );

#   undef ACCUMULATE_J

    ja = PTR_ACCUMULATOR(i0);
    INCREMENT_4x1( ja->jx, a0x );
    INCREMENT_4x1( ja->jy, a0y );
    INCREMENT_4x1( ja->jz, a0z );

    ja = PTR_ACCUMULATOR(i1);
    INCREMENT_4x1( ja->jx, a1x );
    INCREMENT_4x1( ja->jy, a1y );
    INCREMENT_4x1( ja->jz, a1z );

    ja = PTR_ACCUMULATOR(i2);
    INCREMENT_4x1( ja->jx, a2x );
    INCREMENT_4x1( ja->jy, a2y );
    INCREMENT_4x1( ja->jz, a2z );

    ja = PTR_ACCUMULATOR(i3);
    INCREMENT_4x1( ja->jx, a3x );
    INCREMENT_4x1( ja->jy, a3y );
    INCREMENT_4x1( ja->jz, a3z );

    // Update position and accumulate outbnd

#   define MOVE_OUTBND( N )                   \
    if( EXTRACT( outbnd, N ) ) {              \
      pm->dispx = EXTRACT( ddx, N );          \
      pm->dispy = EXTRACT( ddy, N );          \
      pm->dispz = EXTRACT( ddz, N );          \
      pm->i     = idx + N;                    \
      if( move_p_spu( p+N, pm ) ) pm++, nm++; \
    }

    MOVE_OUTBND(0);
    MOVE_OUTBND(1);
    MOVE_OUTBND(2);
    MOVE_OUTBND(3);

#   undef MOVE_OUTBND
    
  }

  return nm;
}

///////////////////////////////////////////////////////////////////////////////
// main (workload distribution and data buffering)

// FIXME: util functionality is not compiled for the spu
// FIXME: UNIQUE EXTENSIONS FOR SPU OBJS AND EXECS!

int
main( uint64_t spu_id,
      uint64_t argp,
      uint64_t envp ) {
  particle_t       * ALIGNED(128) p_block[3];
  int idx[3];
  int np_block[3];

  particle_mover_t * ALIGNED(128) pm_block[3];
  int nm_block[3];

  int buffer;
  const int next[3] = { 1, 2, 0 };
  const int prev[3] = { 2, 0, 1 };

  int np, next_idx, itmp;

# ifdef IN_HARNESS
  prof_clear();
  prof_start();
# endif

  // Get the pipeline arguments from the dispatcher

  mfc_get( args,
           argp,
           sizeof(*args),
           30, 0, 0 );
  mfc_write_tag_mask( (1<<30) );
  mfc_read_tag_status_all();

  pipeline_rank = envp & 0xffffffff;
  n_pipeline    = envp >> 32; // Note: pipeline_rank<n_pipeline

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, next_idx, np );

  // Determine which movers are reserved for this pipeline
  // Movers (16 bytes) are reserved for pipelines in multiples of 8
  // such that the set of particle movers reserved for a pipeline is
  // 128-bit aligned and a multiple of 128-bits in size. 

  args->max_nm -= args->np&15; // Insure host gets enough
  if( args->max_nm<0 ) args->max_nm = 0;
  DISTRIBUTE( args->max_nm, 8, pipeline_rank, n_pipeline, itmp, seg->max_nm );
  seg->pm        = args->pm + itmp*sizeof(particle_mover_t);
  seg->nm        = 0;
  seg->n_ignored = 0;

  // Determine which accumulator array is reserved for this pipeline

  args->a0 += sizeof(accumulator_t)*(1+pipeline_rank)*
    POW2_CEIL((args->nx+2)*(args->ny+2)*(args->nz+2),2);
  
  // Process the particles assigned to this pipeline with triple buffering

# define BEGIN_GET_PBLOCK(buffer) do {                          \
                                                                \
    /* Determine the start and size of the block */             \
    idx[buffer]      = next_idx;                                \
    np_block[buffer] = NP_BLOCK_TARGET;                         \
    if( np_block[buffer]>np ) np_block[buffer] = np;            \
    next_idx += np_block[buffer];                               \
    np       -= np_block[buffer];                               \
                                                                \
    /* If we have a block, start reading it into the buffer */  \
    if( np_block[buffer] )                                      \
      mfc_get( p_block[buffer],                                 \
               args->p0 + idx[buffer]*sizeof(particle_t),       \
               np_block[buffer]*sizeof(particle_t),             \
               3*(buffer)+0, 0, 0 );                            \
  } while(0)
  
# define END_GET_PBLOCK(buffer) do {                                    \
    /* If we have a block, stop reading it into the buffer */           \
    if( np_block[buffer] ) {                                            \
      mfc_write_tag_mask( (np_block[buffer]?1:0)<<(3*(buffer)+0) );     \
      mfc_read_tag_status_all();                                        \
    }                                                                   \
  } while(0)
  
# define PROCESS_PBLOCK(buffer)                                         \
  nm_block[buffer] = advance_p_pipeline_spu( p_block[buffer],           \
                                             pm_block[buffer],          \
                                             idx[buffer],               \
                                             np_block[buffer]>>2 )

  // FIXME: mfc list for this??
# define BEGIN_PUT_PBLOCK(buffer) do {                                  \
                                                                        \
    /* If we have a block, begin writing the buffer into the block */   \
    if( np_block[buffer] )                                              \
      mfc_put( p_block[buffer],                                         \
               args->p0 + idx[buffer]*sizeof(particle_t),               \
               np_block[buffer]*sizeof(particle_t),                     \
               3*(buffer)+1, 0, 0 );                                    \
                                                                        \
    /* Begin writing the movers corresponding to this block */          \
    /* Ignore movers that would overflow the mover segment */           \
    itmp = seg->nm + nm_block[buffer] - seg->max_nm;                    \
    if( itmp > 0 ) {                                                    \
      seg->n_ignored += itmp;                                           \
      nm_block[buffer] = seg->max_nm - seg->nm;                         \
    }                                                                   \
    if( nm_block[buffer] ) {                                            \
      mfc_put( pm_block[buffer],                                        \
               seg->pm + seg->nm*sizeof(particle_mover_t),              \
               nm_block[buffer]*sizeof(particle_mover_t),               \
               3*(buffer)+2, 0, 0 );                                    \
      seg->nm += nm_block[buffer];                                      \
    }                                                                   \
                                                                        \
  } while(0)

# define END_PUT_PBLOCK(buffer) do {                                    \
    if( np_block[buffer] || nm_block[buffer] ) {                        \
      mfc_write_tag_mask( (np_block[buffer]?1:0)<<(3*(buffer)+1) |      \
                          (nm_block[buffer]?1:0)<<(3*(buffer)+2) );     \
      mfc_read_tag_status_all();                                        \
    }                                                                   \
    np_block[buffer] = 0;                                               \
    nm_block[buffer] = 0;                                               \
  } while(0)

  p_block[0]  = local_p;                      np_block[0] = 0;
  p_block[1]  = local_p  + NP_BLOCK_TARGET;   np_block[1] = 0;
  p_block[2]  = local_p  + NP_BLOCK_TARGET*2; np_block[2] = 0;

  pm_block[0] = local_pm;                     nm_block[0] = 0;
  pm_block[1] = local_pm + NP_BLOCK_TARGET;   nm_block[1] = 0;
  pm_block[2] = local_pm + NP_BLOCK_TARGET*2; nm_block[2] = 0;

  BEGIN_GET_PBLOCK(0);
  for( buffer=0; np_block[buffer]>0; buffer=next[buffer] ) {
    BEGIN_GET_PBLOCK( next[buffer] );
    END_GET_PBLOCK( buffer );
    PROCESS_PBLOCK( buffer );
    BEGIN_PUT_PBLOCK( buffer );
    END_PUT_PBLOCK( prev[buffer] );
  }
  END_PUT_PBLOCK( prev[buffer] );

  // Flush the caches. Only accumulator_cache needs to be flushed as
  // the other caches are read-only.

  cache_flush( accumulator_cache );

  // Write the pipeline return values back

  mfc_put( seg,
           args->seg + pipeline_rank*sizeof(particle_mover_seg_t),
           sizeof(particle_mover_seg_t),
           31, 0, 0 );
  mfc_write_tag_mask( (1<<31) );
  mfc_read_tag_status_all();
  
# ifdef IN_HARNESS
  prof_stop();
# endif

  return 0;
}

#endif // CELL_SPU_BUILD
