#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include "../sf_interface_private.h"

#include <spu_mfcio.h>

// DMA tag usage:
// Tags  0:15 - stencil cache
// Tags 16:23 -  input stream
// Tags 24:31 - output stream

// Define the buffered I/O streams

#define IN_TAG(b)   (16+(b))
#define OUT_TAG(b)  (24+(b))
#define N_BUF 8 /* Must be a power of two */

typedef struct ls_neighbor {
  vec_llong2 n01, n23, n45;
} ls_neighbor_t;

typedef struct ls_interp {
  vec_float4 ex, ey, ez, bxy, bz;
  vec_llong2 n01, n23, n45;
} ls_interp_t;

// Define the stencil cache

typedef struct ls_field {
  vec_float4 e, b;
} ls_field_t;

#define CACHED_TYPE        ls_field_t
#define CACHE_NAME         stencil
#define CACHELINE_LOG2SIZE 5  /* 1 local field per line */
#define CACHE_LOG2NSETS    10 /* 1024 line per way */
#define CACHE_LOG2NWAY     2  /* 4-way (one way per non-unit stride) */
#define CACHE_TYPE         0  /* Read only */    
#define CACHE_SET_TAGID(set) ((set)&0xf) /* Cache will use channels 0:15 */
#include "../../util/overlays/cache-api.h"

// Spu fun stuff (FIXME: PROBABLY SHOULD BE MOVED TO A COMMON PLACE)

#define OFFSET_OF(t,f) ((size_t)(&(((t *)(NULL))->f)))

#define PERM4(a,b, i,j,k,l) spu_shuffle( (a), (b), ((vec_uchar16)           \
    { 4*(i), 4*(i)+1, 4*(i)+2, 4*(i)+3, 4*(j), 4*(j)+1, 4*(j)+2, 4*(j)+3,   \
      4*(k), 4*(k)+1, 4*(k)+2, 4*(k)+3, 4*(l), 4*(l)+1, 4*(l)+2, 4*(l)+3 }) )

#define PERM2(a,b, i,j) spu_shuffle( (a), (b), ((vec_uchar16)                 \
    { 8*(i), 8*(i)+1, 8*(i)+2, 8*(i)+3, 8*(i)+4, 8*(i)+5, 8*(i)+6, 8*(i)+7,   \
      8*(j), 8*(j)+1, 8*(j)+2, 8*(j)+3, 8*(j)+4, 8*(j)+5, 8*(j)+6, 8*(j)+7 }) )

#define VEC_DOUBLE2(a,i, b,j) spu_extend( PERM4((a),(b),(i),(i),(j)+4,(j)+4) )

void
_SPUEAR_load_interpolator_pipeline_spu(
    load_interpolator_pipeline_args_t * args,
    int pipeline_rank,
    int n_pipeline ) {
  MEM_PTR( interpolator_t, 128 ) fi = args->fi;
  MEM_PTR( int64_t,        128 ) nb = args->nb;
  cache_init( stencil, args->f + OFFSET_OF(field_t,ex) );

  ls_neighbor_t * RESTRICT ALIGNED(128) nb_in;  SPU_MALLOC(nb_in, N_BUF,128);
  ls_interp_t   * RESTRICT ALIGNED(128) fi_out; SPU_MALLOC(fi_out,N_BUF,128);

  ls_field_t f0;
  ls_field_t fx, fy, fz;
  ls_field_t fyz, fzx, fxy;
  int v, x, y, z, n_voxel;

  int nx = args->nx, sx = 1;
  int ny = args->ny, sy = (nx+2)*sx;
  int nz = args->nz, sz = (ny+2)*sy;

  vec_float4 fourth = (vec_float4){ 0.25f, 0.25f, 0.25f, 0.25f };
  vec_float4 half   = (vec_float4){ 0.5f,  0.5f,  0.5f,  0.5f  };

  vec_double2 vd0, vd1, vd2, vd3, s30, d30, s12, d12;
  vec_float4  vf0, vf1, vf2, vf3;
  vec_float4  ex, ey, ez, bxy, bz;

  int V, X, Y, Z, b;

  // Determine which voxels are assigned here

  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 1, pipeline_rank, n_pipeline,
                     x,y,z,n_voxel ); v = VOXEL(x,y,z, nx,ny,nz);
# define VOXEL_INC(v,x,y,z) NEXT_VOXEL(v,x,y,z, 1,nx, 1,ny, 1,nz, nx,ny,nz)

  // Fill up the voxel pipeline

  V = v; X = x; Y = y; Z = z;
  for( b=0; b<N_BUF; b++ )
    if( LIKELY( n_voxel>b ) ) {
      mfc_get( nb_in + b, nb + V*6*sizeof(int64_t),
               6*sizeof(int64_t), IN_TAG(b), 0, 0 );
      cache_touch( stencil,  V       *sizeof(field_t) );
      cache_touch( stencil, (V+sy   )*sizeof(field_t) );
      cache_touch( stencil, (V   +sz)*sizeof(field_t) );
      cache_touch( stencil, (V+sy+sz)*sizeof(field_t) );
      VOXEL_INC(V,X,Y,Z);
    }

  // For each voxel assigned here

  b = 0;
  for( ; n_voxel; n_voxel-- ) {

    // Load the stencil inputs

    f0  = cache_wait_rd( stencil,  v          *sizeof(field_t), 1 );
    fx  = cache_wait_rd( stencil, (v+sx      )*sizeof(field_t), 0 );
    fy  = cache_wait_rd( stencil, (v   +sy   )*sizeof(field_t), 0 );
    fxy = cache_wait_rd( stencil, (v+sx+sy   )*sizeof(field_t), 0 );
    fz  = cache_wait_rd( stencil, (v      +sz)*sizeof(field_t), 0 );
    fzx = cache_wait_rd( stencil, (v+sx   +sz)*sizeof(field_t), 0 );
    fyz = cache_wait_rd( stencil, (v   +sy+sz)*sizeof(field_t), 0 );

    // Compute the stencil values

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
    vd1 = spu_xor( PERM2(vd0,vd0, 1,0),
                   ((vec_double2){ +0.,-0.}) );   //       w1z              -w0z
    vf0 = spu_roundtf( spu_add( vd0, vd1 ) );     //    cbz      0       dcbzdz    0
    bz  = spu_mul( half, PERM4(vf0,vf0, 0,2,1,3) );

#   undef z
#   undef y
#   undef x

    // Finishing reading in the input voxel data, reduce the stencil
    // values to it and starting writing the output voxel data

    mfc_write_tag_mask( (1<<IN_TAG(b)) | (1<<OUT_TAG(b)) );
    mfc_read_tag_status_all();
    fi_out[b].ex  = ex;
    fi_out[b].ey  = ey;
    fi_out[b].ez  = ez;
    fi_out[b].bxy = bxy;
    fi_out[b].bz  = bz;
    fi_out[b].n01 = nb_in[b].n01;
    fi_out[b].n23 = nb_in[b].n23;
    fi_out[b].n45 = nb_in[b].n45;
    mfc_put( fi_out + b, fi + v*sizeof(interpolator_t),
             sizeof(interpolator_t), OUT_TAG(b), 0, 0 );
    if( LIKELY( n_voxel>N_BUF ) ) {
      mfc_get( nb_in + b, nb + V*6*sizeof(int64_t),
               6*sizeof(int64_t), IN_TAG(b), 0, 0 );
      cache_touch( stencil,  V       *sizeof(field_t) );
      cache_touch( stencil, (V+sy   )*sizeof(field_t) );
      cache_touch( stencil, (V   +sz)*sizeof(field_t) );
      cache_touch( stencil, (V+sy+sz)*sizeof(field_t) );
      VOXEL_INC(V,X,Y,Z);
    }

    // Advance to the next voxel

    VOXEL_INC(v,x,y,z); b = (b+1) & (N_BUF-1);
  }
}
