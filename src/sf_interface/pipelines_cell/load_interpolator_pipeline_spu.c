#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include "../sf_interface_private.h"

#include <spu_mfcio.h>

typedef struct local_field {
  vec_float4 e, b;
} local_field_t;

typedef struct local_interpolator {
  vec_float4 ex, ey, ez, bxy, bz;
  vec_llong2 n01, n23, n45;
} local_interpolator_t;

// Stencil cache setup

#define CACHED_TYPE        local_field_t
#define CACHE_NAME         stencil
#define CACHELINE_LOG2SIZE 5  /* 1 local field per line */
#define CACHE_LOG2NSETS    10 /* 1024 line per way */
#define CACHE_LOG2NWAY     2  /* 4-way (one way per non-unit stride) */
#define CACHE_TYPE         0  /* Read only */    
#define CACHE_SET_TAGID(set) ((set)&0xf) /* Cache will use channels 0:15 */
#include "../../util/overlays/cache-api.h"

#define CACHE_EA_BASE      args->f
#define CACHE_EA_STRIDE    sizeof(field_t)

#define CACHE_INIT() cache_init( stencil, CACHE_EA_BASE )
#define CACHE_RD(i)  (*cache_wait_rw( stencil, (i)*CACHE_EA_STRIDE ))

// Stencil prefetch setup

#define PREFETCH_LOOKAHEAD 16

#define PREFETCH(i) cache_touch( stencil, (i)*CACHE_EA_STRIDE )

// Output writeback setup

#define WB_BUF       local_fi /* Local writeback buffers */
#define WB_N_BUF     4        /* Number of writeback buffers */
#define WB_EA_BASE   fi       /* Writeback base external address */
#define WB_EA_STRIDE sizeof(interpolator_t) /* Writeback ext byte stride */
#define WB_TAG_BASE  16       /* Use tags WB_TAG:WB_TAG+WB_N_BUF-1 */

#define WB_BEGIN(b,i) mfc_put( &WB_BUF[(b)], WB_EA_BASE+(i)*WB_EA_STRIDE, \
                               sizeof(WB_BUF[0]), WB_TAG_BASE+(b), 0, 0 )
#define WB_END(b)     mfc_write_tag_mask( 1 << (WB_TAG_BASE+(b)) );     \
                      mfc_read_tag_status_all()
#define WB_INC(b)     b++; if( b==WB_N_BUF ) b=0

// Spu fun stuff (FIXME: PROBABLY SHOULD BE MOVED TO A COMMON PLACE)

#define PERM4( a, b, i, j, k, l ) spu_shuffle( (a), (b), ((vec_uchar16)   \
  { 4*(i), 4*(i)+1, 4*(i)+2, 4*(i)+3, 4*(j), 4*(j)+1, 4*(j)+2, 4*(j)+3,   \
    4*(k), 4*(k)+1, 4*(k)+2, 4*(k)+3, 4*(l), 4*(l)+1, 4*(l)+2, 4*(l)+3 }) )

#define PERM2( a, b, i, j ) spu_shuffle( (a), (b), ((vec_uchar16)           \
  { 8*(i), 8*(i)+1, 8*(i)+2, 8*(i)+3, 8*(i)+4, 8*(i)+5, 8*(i)+6, 8*(i)+7,   \
    8*(j), 8*(j)+1, 8*(j)+2, 8*(j)+3, 8*(j)+4, 8*(j)+5, 8*(j)+6, 8*(j)+7 }) )

#define pos +0.
#define neg -0.
#define NEG2( a, s0, s1 ) spu_xor( (a), ((vec_double2){ s0, s1 }) )

#define VEC_DOUBLE2( a,i, b,j ) spu_extend( PERM4((a),(b),(i),(i),(j)+4,(j)+4) )

void
_SPUEAR_load_interpolator_pipeline_spu(
    load_interpolator_pipeline_args_t * args,
    int pipeline_rank,
    int n_pipeline ) {
  MEM_PTR( interpolator_t, 128 ) fi = args->fi;
  MEM_PTR( int64_t,        128 ) nb = args->nb;

  local_interpolator_t * RESTRICT ALIGNED(128) local_fi;

  local_field_t f0;
  local_field_t fx, fy, fz;
  local_field_t fyz, fzx, fxy;
  int x, y, z, v, n_voxel;

  const int nx = args->nx, sx = 1;
  const int ny = args->ny, sy = (nx+2)*sx;
  const int nz = args->nz, sz = (ny+2)*sy;

  vec_float4 fourth = (vec_float4){ 0.25f, 0.25f, 0.25f, 0.25f };
  vec_float4 half   = (vec_float4){ 0.5f,  0.5f,  0.5f,  0.5f  };

  vec_double2 vd0, vd1, vd2, vd3, s30, d30, s12, d12;
  vec_float4  vf0, vf1, vf2, vf3;
  vec_float4  ex, ey, ez, bxy, bz;
  vec_llong2 local_nb[3];

  int X, Y, Z, V, b;

  // Process voxels assigned to this pipeline

  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 1, pipeline_rank, n_pipeline,
                     x,y,z,n_voxel ); v = VOXEL(x,y,z, nx,ny,nz);
# define VOXEL_INC(v,x,y,z) NEXT_VOXEL(v,x,y,z, 1,nx, 1,ny, 1,nz, nx,ny,nz)

  // Setup cache, prefetching and writebacks
  CACHE_INIT();
  V = v; X = x; Y = y; Z = z;
  for( b=0; b<PREFETCH_LOOKAHEAD; b++ )
    if( LIKELY( b<n_voxel ) ) {
      PREFETCH( V           );
      PREFETCH( V + sy      );
      PREFETCH( V + sz      );
      PREFETCH( V + sy + sz );
      VOXEL_INC(V,X,Y,Z);
    }
  SPU_MALLOC( local_fi, WB_N_BUF, 128 );
  b = 0;

  for( ; n_voxel; n_voxel-- ) {
    mfc_get( local_nb, nb+6*sizeof(int64_t)*v, 6*sizeof(int64_t), 31, 0, 0 );

    // Prefetch voxels likely to be used in the future
    if( LIKELY( n_voxel>PREFETCH_LOOKAHEAD ) ) {
      PREFETCH( V           );
      PREFETCH( V + sy      );
      PREFETCH( V + sz      );
      PREFETCH( V + sy + sz );
      VOXEL_INC(V,X,Y,Z);
    }

    // Read stencil inputs
    f0  = CACHE_RD( v           ); fx  = CACHE_RD( v + sx           );
    fy  = CACHE_RD( v + sy      ); fxy = CACHE_RD( v + sx + sy      );
    fz  = CACHE_RD( v      + sz ); fzx = CACHE_RD( v + sx      + sz );
    fyz = CACHE_RD( v + sy + sz );

#   define x 0
#   define y 1
#   define z 2

    // Compute ex and ez interpolation coefficients
    vd0 =  spu_extend(  f0.e            );        //       w0x               w0z
    vd1 = VEC_DOUBLE2(  fy.e,x,  fx.e,z );        //       w1x               w1z
    vd2 = VEC_DOUBLE2(  fz.e,x,  fy.e,z );        //       w2x               w2z
    vd3 = VEC_DOUBLE2( fyz.e,x, fxy.e,z );        //       w3x               w3z
    s30 = spu_add( vd3, vd0 );                    //     w3x+w0x           w3z+w0z
    d30 = spu_sub( vd3, vd0 );                    //     w3x-w0x           w3z-w0z
    s12 = spu_add( vd1, vd2 );                    //     w1x+w2x           w1z+w2z
    d12 = spu_sub( vd1, vd2 );                    //     w1x-w2x           w1z-w2z
    vf0 = spu_roundtf( spu_add( s30, s12 ) );     //   ex        0        ez       0
    vf1 = spu_roundtf( spu_add( d30, d12 ) );     //  dexdy      0       dezdx     0
    vf2 = spu_roundtf( spu_sub( d30, d12 ) );     //  dexdz      0       dezdy     0
    vf3 = spu_roundtf( spu_sub( s30, s12 ) );     // d2exdydz    0      d2ezdxdy   0
    vf0 = spu_or( vf0, PERM4(vf1,vf1, 1,0,3,2) ); //   ex       dexdy     ez      dezdx
    vf1 = spu_or( vf2, PERM4(vf3,vf3, 1,0,3,2) ); //  dexdz    d2exdydz  dezdy   d2ezdxdy
    ex  = spu_mul( fourth, PERM4(vf0,vf1, 0,1,4,5) );
    ez  = spu_mul( fourth, PERM4(vf0,vf1, 2,3,6,7) );

    // Compute ey interpolation coefficients
    vd0 = VEC_DOUBLE2( fzx.e,y,  fz.e,y );        //       w3y               w1y
    vd1 = VEC_DOUBLE2(  f0.e,y,  fx.e,y );        //       w0y               w2y
    vd2 = spu_add( vd0, vd1 );                    //     w3y+w0y           w1y+w2y
    vd3 = spu_sub( vd0, vd1 );                    //     w3y-w0y           w1y-w2y
    vd0 = PERM2(vd2,vd3, 0,2);                    //     w3y+w0y           w3y-w0y
    vd1 = PERM2(vd2,vd3, 1,3);                    //     w1y+w2y           w1y-w2y
    vf0 = spu_roundtf( spu_add( vd0, vd1 ) );     //   ey        0       deydz     0
    vf1 = spu_roundtf( spu_sub( vd0, vd1 ) );     // d2eydzdx    0       deydx     0
    ey  = spu_mul( fourth, PERM4(vf0,vf1, 0,2,6,4) );

    // Compute bx and by interpolation coefficients
    vd0 = VEC_DOUBLE2(  f0.b,x,  f0.b,y );        //       w0x               w0y
    vd1 = VEC_DOUBLE2(  fx.b,x,  fy.b,y );        //       w1x               w1y
    vf0 = spu_roundtf( spu_add( vd1, vd0 ) );     //    cbx      0        cby      0
    vf1 = spu_roundtf( spu_sub( vd1, vd0 ) );     //   dcbdx     0       dcbydy    0
    bxy = spu_mul( half, PERM4(vf0,vf1, 0,4,2,6) );

    // Compute bz interpolation coefficients
    vd0 = VEC_DOUBLE2(  f0.b,z,  fz.b,z );        //       w0z               w1z
    vd1 = NEG2( PERM2(vd0,vd0, 1,0), pos,neg );   //       w1z              -w0z
    vf0 = spu_roundtf( spu_add( vd0, vd1 ) );     //    cbz      0       dcbzdz    0
    bz  = spu_mul( half, PERM4(vf0,vf0, 0,2,1,3) );

#   undef z
#   undef y
#   undef x

    // Write stencil outputs
    WB_END(b);
    local_fi[b].ex  = ex;
    local_fi[b].ey  = ey;
    local_fi[b].ez  = ez;
    local_fi[b].bxy = bxy;
    local_fi[b].bz  = bz;
    mfc_write_tag_mask( 1 << 31 ); mfc_read_tag_status_all();
    local_fi[b].n01 = local_nb[0];
    local_fi[b].n23 = local_nb[1];
    local_fi[b].n45 = local_nb[2];
    WB_BEGIN(b,v); WB_INC(b);

    // Advance to the next voxel
    VOXEL_INC(v,x,y,z);
  }
}
