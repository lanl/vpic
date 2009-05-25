#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include "../sf_interface_private.h"

#include <spu_mfcio.h>

// DMA tag usage:
// Tags  0:15 - stencil cache
// Tags 16:23 -  input stream
// Tags 24:31 - output stream

#define IN_TAG(b)  (16+(b))
#define OUT_TAG(b) (24+(b))
#define N_BUF 8 /* Must be a power of two */

// Define the stencil cache
// FIXME: COULD MAKE CACHE MORE EFFICIENT AND IGNORE THIS PADDING IN
// THE ACCUMULATOR FOR 64-BYTE ALIGNMENT

typedef struct local_accumulator {
  vec_float4 jx;
  vec_float4 jy;
  vec_float4 jz;
  vec_float4 pad;
} local_accumulator_t;

#define CACHED_TYPE        local_accumulator_t
#define CACHE_NAME         stencil
#define CACHELINE_LOG2SIZE 6  /* 1 local accumulator per line */
#define CACHE_LOG2NSETS    9  /* 9 line per way */
#define CACHE_LOG2NWAY     2  /* 4-way (one way per non-unit stride) */
#define CACHE_TYPE         0  /* Read only */    
#define CACHE_SET_TAGID(set) ((set)&0xf) /* Cache will use channels 0:15 */
#include "../../util/overlays/cache-api.h"

#define UNCACHE(v) /* FIXME: IMPLEMENT THIS */

// Spu fun stuff (FIXME: PROBABLY SHOULD BE MOVED TO A COMMON PLACE)

#define PERM4(a,b, i,j,k,l) spu_shuffle((a),(b), ((vec_uchar16)         \
    { 4*(i), 4*(i)+1, 4*(i)+2, 4*(i)+3, 4*(j), 4*(j)+1, 4*(j)+2, 4*(j)+3, \
      4*(k), 4*(k)+1, 4*(k)+2, 4*(k)+3, 4*(l), 4*(l)+1, 4*(l)+2, 4*(l)+3 }))

#define PERM2(a,b, i,j) spu_shuffle( (a), (b), ((vec_uchar16)           \
    { 8*(i), 8*(i)+1, 8*(i)+2, 8*(i)+3, 8*(i)+4, 8*(i)+5, 8*(i)+6, 8*(i)+7, \
      8*(j), 8*(j)+1, 8*(j)+2, 8*(j)+3, 8*(j)+4, 8*(j)+5, 8*(j)+6, 8*(j)+7 }) )

#define VEC_DOUBLE2(a,i, b,j) spu_extend( PERM4((a),(b),(i),(i),(j)+4,(j)+4) )

void
_SPUEAR_unload_accumulator_pipeline_spu(
    unload_accumulator_pipeline_args_t * args,
    int pipeline_rank,
    int n_pipeline ) {
  MEM_PTR( field_t, 128 ) f = args->f;
  cache_init( stencil, args->a );

  vec_float4 * RESTRICT ALIGNED(128) jf_in;  SPU_MALLOC( jf_in,  N_BUF, 128 );
  vec_float4 * RESTRICT ALIGNED(128) jf_out; SPU_MALLOC( jf_out, N_BUF, 128 );

  local_accumulator_t  a0;
  local_accumulator_t  ax,  ay,  az;
  local_accumulator_t ayz, azx, axy;
  int x, y, z, v, n_voxel;

  int nx = args->nx, sx = 1;
  int ny = args->ny, sy = (nx+2)*sx;
  int nz = args->nz, sz = (ny+2)*sy;

  vec_double2 cxz = { (double)args->cx, (double)args->cz };
  vec_double2 cy0 = { (double)args->cy, 0.               };
  vec_double2 jxz, jyw;
  vec_float4 f0;

  int X, Y, Z, V, b;

  // Process voxels assigned to this pipeline
  DISTRIBUTE_VOXELS( 1,nx+1, 1,ny+1, 1,nz+1, 1, pipeline_rank, n_pipeline,
                     x,y,z,n_voxel ); v = VOXEL(x,y,z, nx,ny,nz);
# define VOXEL_INC(v,x,y,z) NEXT_VOXEL(v,x,y,z, 1,nx, 1,ny, 1,nz, nx,ny,nz)

  // Fill up the voxel pipeline

  V = v; X = x; Y = y; Z = z;
  for( b=0; b<N_BUF; b++ )
    if( LIKELY( n_voxel>b ) ) {
      mfc_get( jf_in + b, f + V*sizeof(field_t) + 3*4*sizeof(float),
               4*sizeof(float), IN_TAG(b), 0, 0 );
      cache_touch( stencil, (V-sy-sz)*sizeof(accumulator_t) );
      cache_touch( stencil, (V   -sz)*sizeof(accumulator_t) );
      cache_touch( stencil, (V-sy   )*sizeof(accumulator_t) );
      cache_touch( stencil,  V       *sizeof(accumulator_t) );
      VOXEL_INC(V,X,Y,Z);
    }

  // For each voxel

  b = 0;
  for( ; n_voxel; n_voxel-- ) {

    // Load the stencil inputs
    ayz = *cache_wait_rd( stencil, (v   -sy-sz)*sizeof(accumulator_t) );
    azx = *cache_wait_rd( stencil, (v-sx   -sz)*sizeof(accumulator_t) );
    az  = *cache_wait_rd( stencil, (v      -sz)*sizeof(accumulator_t) );
    axy = *cache_wait_rd( stencil, (v-sx-sy   )*sizeof(accumulator_t) );
    ay  = *cache_wait_rd( stencil, (v   -sy   )*sizeof(accumulator_t) );
    ax  = *cache_wait_rd( stencil, (v-sx      )*sizeof(accumulator_t) );
    a0  = *cache_wait_rd( stencil,  v          *sizeof(accumulator_t) );
    
    // Compute the stencil values
    jxz = spu_mul( cxz, spu_add( spu_add( VEC_DOUBLE2( a0.jx,0,  a0.jz,0),
                                          VEC_DOUBLE2( ay.jx,1,  ax.jz,1)),
                                 spu_add( VEC_DOUBLE2( az.jx,2,  ay.jz,2),
                                          VEC_DOUBLE2(ayz.jx,3, axy.jz,3))));
    /**/                  jyw =  spu_add( VEC_DOUBLE2( a0.jy,0,  ax.jy,2),
                                          VEC_DOUBLE2( az.jy,1, azx.jy,3) );
    jyw = spu_mul( cy0, spu_add( PERM2(jyw,jyw, 0,0), PERM2(jyw,jyw, 1,1) ) );

    // Finishing reading in the input voxel data, reduce the stencil to it
    // and starting writing the output voxel data
    mfc_write_tag_mask( (1<<IN_TAG(b)) | (1<<OUT_TAG(b)) );
    mfc_read_tag_status_all();
    f0 = jf_in[b];
    jf_out[b] = PERM4( spu_roundtf( spu_add( spu_extend( f0        ), jxz ) ),
                       spu_roundtf( spu_add( VEC_DOUBLE2(f0,1, f0,3), jyw ) ),
                       0,4,2,6 );
    mfc_put( jf_out+b, f + v*sizeof(field_t) + 3*4*sizeof(float),
             4*sizeof(float), OUT_TAG(b), 0, 0 );

    if( LIKELY( n_voxel>N_BUF ) ) {
      // Discard items we know aren't useful anymore
      UNCACHE( v      - sy - sz );
      UNCACHE( v - sx      - sz );
      UNCACHE( v - sx - sy      );
      
      // Prefetch items we know we will need
      mfc_get( &f_in[b].jf, f + V*sizeof(field_t) + 3*4*sizeof(float),
               4*sizeof(float), IN_TAG(b), 0, 0 );
      cache_touch( stencil, (V-sy-sz)*sizeof(accumulator_t) );
      cache_touch( stencil, (V   -sz)*sizeof(accumulator_t) );
      cache_touch( stencil, (V-sy   )*sizeof(accumulator_t) );
      cache_touch( stencil,  V       *sizeof(accumulator_t) );
      VOXEL_INC(V,X,Y,Z);
    }

    // Advance to the next voxel
    VOXEL_INC(v,x,y,z); b = (b+1) & (N_BUF-1);
  }
}

#if 0
#include <stdio.h>
#define DELAY() for( volatile int delay[1]={1000000}; delay[0]; delay[0]-- )
#define INSPECT_F(v) do { printf( #v"=%f\n", v ); fflush(stdout); DELAY(); } while(0)
#define INSPECT_I(v) do { printf( #v"=%i\n", v ); fflush(stdout); DELAY(); } while(0)
#define INSPECT_VF(v) do { printf( #v"=%f %f %f %f\n", spu_extract(v,0), spu_extract(v,1), spu_extract(v,2), spu_extract(v,3) ); fflush(stdout); DELAY(); } while(0)
#define INSPECT_VD(v) do { printf( #v"=%f %f\n", spu_extract(v,0), spu_extract(v,1) ); fflush(stdout); DELAY(); } while(0)
#endif
