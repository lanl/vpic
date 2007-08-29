#define IN_particle_pipeline
#define HAS_SPU_PIPELINE
#include <particle_pipelines.h>
#include <v4c_spu.h>
#include <spu_mfcio.h>

#ifdef IN_HARNESS
#include <profile.h>
#endif

#ifndef likely
#define likely(_c)      __builtin_expect((_c), 1)
#define unlikely(_c)    __builtin_expect((_c), 0)
#endif

// DMA tag usage:
//  0: 2 - Particle buffer 0 (read, write, mover write)
//  3: 5 - Particle buffer 1 (read, write, mover write)
//  6: 8 - Particle buffer 2 (read, write, mover write)
//  9:16 - Interpolator cache
// 17:25 - Accumulator cache
// 26:29 - Unused presently
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

// NP_BLOCK=512 is such that a 16Kb of particle data is located in a
// particle block (the maximum amount of data that can be moved in a
// single DMA transfer).  Triple buffering requires 3 such particle
// blocks.  So local particle storage takes up 48Kb.  Each local
// particle has a mover associated with it such that a SPU pipeline
// cannot run out movers during particle processing.  This costs 24Kb.

#define NP_BLOCK 512

DECLARE_ALIGNED_ARRAY( particle_t,       128, p_stream, 3*NP_BLOCK );
DECLARE_ALIGNED_ARRAY( particle_mover_t, 128, m_stream, 3*NP_BLOCK );

///////////////////////////////////////////////////////////////////////////////
// External memory caches - custom protocol

// VOXEL_CACHE_N_LINE gives the number of voxels to cache locally.
// This must be a power of two no less than NP_BLOCK.

#define VOXEL_CACHE_N_LINE 512

// Storage for locally cached interpolators and accumulators.  Since
// particle blocks are no more than N_CACHE_N_LINE in size, there is
// enough room in the cache to accomodate the worst case situation of
// every particle in a particle block being in different voxels.

DECLARE_ALIGNED_ARRAY( interpolator_t, 128, i_cache, VOXEL_CACHE_N_LINE );
DECLARE_ALIGNED_ARRAY( accumulator_t,  128, a_cache, VOXEL_CACHE_N_LINE );

// The cache tags. voxel_cache_addr gives the external address
// (actually a local voxel index) of the interpolator and accumulator
// data that is cached there.  It is 0xffffffff if there is cache line
// is empty.  voxel_cache_next gives the cache line that will be
// replaced after this one.  voxel_cache_prev gives the cache line
// that will be replaced immediately before this one.

DECLARE_ALIGNED_ARRAY( uint32_t, 128, voxel_cache_addr, VOXEL_CACHE_N_LINE );
DECLARE_ALIGNED_ARRAY( uint32_t, 128, voxel_cache_prev, VOXEL_CACHE_N_LINE );
DECLARE_ALIGNED_ARRAY( uint32_t, 128, voxel_cache_next, VOXEL_CACHE_N_LINE );

// The cache hash table.  This provides a hash table that maps keys
// (addresses / voxel index) to values (cache lines).
// VOXEL_CACHE_N_HASH must be a power of two larger than
// VOXEL_CACHE_N_LINE.  VOXEL_CACHE_HASH hashes an external "address"
// (actually a voxel index) into a hash table index.  voxel_cache_key
// and voxel_cache_val give the key-value pairs in the hash table.

#define VOXEL_CACHE_N_HASH 4096

DECLARE_ALIGNED_ARRAY( uint32_t, 128, voxel_cache_key, VOXEL_CACHE_N_HASH );
DECLARE_ALIGNED_ARRAY( uint32_t, 128, voxel_cache_val, VOXEL_CACHE_N_HASH );

#define VOXEL_CACHE_GOLDEN ((uint32_t)2654435761U)
#define VOXEL_CACHE_HASH(addr) \
  (((addr)*VOXEL_CACHE_GOLDEN)&(VOXEL_CACHE_N_HASH-1))

// Cache miscelleanous.  voxel_lru gives the least recently used cache
// line.  {I,A}_CACHE_[FIRST,LAST]_DMA_CHANNEL gives the [FIRST,LAST]
// DMA channels that {interpolator, accumulator} operations should
// use.  There should be a power of two number of DMA channels
// assigned to each type of operation.  a_writeback provides a buffer
// to hold data about pending writeback operations.

uint32_t voxel_cache_lru;

#define I_CACHE_DMA_CHANNEL_FIRST 9
#define I_CACHE_DMA_CHANNEL_LAST  16
#define A_CACHE_DMA_CHANNEL_FIRST 17
#define A_CACHE_DMA_CHANNEL_LAST  25

#define I_CACHE_DMA_CHANNEL(addr) ( I_CACHE_DMA_CHANNEL_FIRST + \
    ( (addr) & ( I_CACHE_DMA_CHANNEL_LAST - I_CACHE_DMA_CHANNEL_FIRST ) ) )

#define A_CACHE_DMA_CHANNEL(addr) ( A_CACHE_DMA_CHANNEL_FIRST + \
    ( (addr) & ( A_CACHE_DMA_CHANNEL_LAST - A_CACHE_DMA_CHANNEL_FIRST ) ) )

DECLARE_ALIGNED_ARRAY( accumulator_t, 128, a_cache_writeback,
                       A_CACHE_DMA_CHANNEL_LAST -
                       A_CACHE_DMA_CHANNEL_FIRST + 1 );

// voxel_cache_init initializes the interpolator/accumulator cache
// data structures.

static __inline void
voxel_cache_init( void ) {
  uint32_t line;

  memset( voxel_cache_addr, 0xff, sizeof(uint32_t)*VOXEL_CACHE_N_LINE );
  memset( voxel_cache_key,  0xff, sizeof(uint32_t)*VOXEL_CACHE_N_HASH );
  memset( voxel_cache_val,  0xff, sizeof(uint32_t)*VOXEL_CACHE_N_HASH );

  for( line=0; line<VOXEL_CACHE_N_LINE; line++ ) {
    voxel_cache_prev[line] = line-1;
    voxel_cache_next[line] = line+1;
  }
  voxel_cache_prev[0                   ] = VOXEL_CACHE_N_LINE - 1;
  voxel_cache_next[VOXEL_CACHE_N_LINE-1] = 0;

  voxel_cache_lru = 0;
}

// voxel_cache_fetch will start fetching the data at "addr" into the
// cache.  The cache line which will store the data is returned.  The
// data is not be guaranteed to be available until voxel_cache_wait is
// called.  WARNING: This routine may silently overwrite cached data
// if used naively.  Read the guarantee that "voxel_cache_wait"
// provides _carefully_.
// 
// The user has the option of providing the hash of the
// voxel_cache_fetch in case it they can use vector instructions to
// produce the hash faster.

#define voxel_cache_fetch( addr ) \
  _voxel_cache_fetch( addr, VOXEL_CACHE_HASH( addr ) )

static __inline uint32_t
_voxel_cache_fetch( uint32_t addr,
                 uint32_t hash_addr ) {
  uint32_t old_addr;
  uint32_t hash_i, hash_j, hash_k;
  uint32_t line, prev, next;

  // Lookup addr in the cache - O(1)

  hash_i = hash_addr;
  while( unlikely( voxel_cache_key[hash_i]!=addr & /* YES! A BIT AND */
                   voxel_cache_key[hash_i]!=0xffffffff ) )
    hash_i = ( hash_i + 1 ) & ( VOXEL_CACHE_N_HASH - 1 );
  line = voxel_cache_val[hash_i];

  if( likely( line!=0xffffffff ) ) { // Cache hit
    
    // Mark this cache line as the most recently used
    
    if( unlikely( voxel_cache_next[line]!=voxel_cache_lru ) ) {
      if( unlikely( line==voxel_cache_lru ) )
        voxel_cache_lru = voxel_cache_next[voxel_cache_lru];
      else {

        // Remove this cache line from the cache line replacement
        // sequence.  Then insert this cache line just before the
        // least recently used line in the replacement sequence.
      
        next = voxel_cache_next[line];
        prev = voxel_cache_prev[line];
        voxel_cache_next[prev] = next;
        voxel_cache_prev[next] = prev;
           
        next = voxel_cache_lru;
        prev = voxel_cache_prev[next];
        voxel_cache_next[line] = next;
        voxel_cache_next[line] = prev;
        voxel_cache_next[prev] = line;
        voxel_cache_prev[next] = line;

      }
    }
    
  } else { // Cache miss ... Unlikely
    
    // We will fill the least recently used line with the data at addr

    line = voxel_cache_lru;

    // Read in the new interpolator for this line

    mfc_get( &i_cache[line], args->f0+sizeof(interpolator_t)*addr,
             sizeof(interpolator_t), I_CACHE_DMA_CHANNEL(addr), 0, 0 );

    // This accumulator DMA operation is _very_ tricky.  The
    // accumulator cache read must have a fence to insure that
    // previous accumulator cache writes have completed before the
    // read. Note that accumulator cache writes to a channel do not
    // need a fence because there will never be two cache writes to
    // the same address without at least one intervening cache read.
    // (Think about it.)  It is desirable to use more than one DMA
    // channel in cache operations but fenced DMA operations only
    // apply to a single DMA channel.  Thus, to make the fence
    // effective such that all writes to a particular address complete
    // before a read of that address, all writes and the reads to that
    // address must use the same DMA channel.  Given the cache is
    // fully associative (there is no relation between the a cache
    // line and the address of the data it holds), cache line
    // information is not useful to determine which DMA channel to use
    // for the writes or reads.  Thus, the address of the transaction
    // must be used to to determine the DMA channel.

    // This brings up a second complexity.  When writing back the old
    // contents of a cache line being replaced, the write is
    // immediately followed by read to that same cache line.  However,
    // because the DMA channel for these transactions is determined by
    // the addresses of these operations, these operations may use
    // different DMA channels.  So that these operations do not
    // collide (i.e. the read of the new data replaces the old data
    // before the old data is written), it is necessary to copy the
    // data into a writeback buffer.

    // This brings up a third complexity.  So that a writeback on does
    // not clobber a previous writeback on a given DMA channel, we
    // need to wait until the channel used for the writeback is free.

    // Flush out the old accumulator as necessary.
   
    old_addr = voxel_cache_addr[line];

    if( likely( old_addr!=0xffffffff ) ) {

      uint32_t channel = A_CACHE_DMA_CHANNEL(old_addr);     

      // Wait for the writeback channel to become free

      mfc_write_tag_mask( 1<< channel );
      mfc_read_tag_status_all();

      // Buffer the data to writeback

      a_cache_writeback[channel] = a_cache[line];

      // Writeback the data

      mfc_put( &a_cache_writeback[channel],
               args->a0+sizeof(accumulator_t)*old_addr,
               sizeof(accumulator_t), channel, 0, 0 );

    }

    // Read in the new accumulator for this line

    mfc_getf( &a_cache[line], args->a0+sizeof(accumulator_t)*addr,
              sizeof(accumulator_t), A_CACHE_DMA_CHANNEL(addr), 0, 0 );

    // Update the least recently used line.  The new least recently
    // used line is the one that follows this line in the cache line
    // replacement sequence.

    voxel_cache_lru = voxel_cache_next[line];

    // Update the cache tag
    
    voxel_cache_addr[line] = addr;

    // If we had to flush previously cached data, remove the map of
    // old_addr to line - O(1).  Note that if we had a write back, we
    // are guaranteed that old_addr _is_ in the hash table.

    if( likely( old_addr!=0xffffffff ) ) {

      hash_i = VOXEL_CACHE_HASH(old_addr);
      while( unlikely( voxel_cache_key[hash_i]!=old_addr ) )
        hash_i = ( hash_i + 1 ) & ( VOXEL_CACHE_N_HASH - 1 );     

      hash_j = hash_i;
      for(;;) {
        hash_j = ( hash_j + 1 ) & ( VOXEL_CACHE_N_HASH - 1 );
        if( likely( voxel_cache_key[hash_j]==0xffffffff ) ) break;
        hash_k = VOXEL_CACHE_HASH( voxel_cache_key[hash_j] );
        if( ( hash_j>hash_i & ( hash_k<=hash_i | hash_k>hash_j ) ) |
            ( hash_j<hash_i & ( hash_k<=hash_i & hash_k>hash_j ) ) ) {
          // YES! BIT OPS!
          voxel_cache_key[hash_i] = voxel_cache_key[hash_j];
          voxel_cache_val[hash_i] = voxel_cache_val[hash_j];
          hash_i = hash_j;
        }
      }

      voxel_cache_key[hash_i] = 0xffffffff;
      voxel_cache_val[hash_i] = 0xffffffff;

    }

    // Map addr to line - O(1).  Note that we are guaranteed that addr
    // _is_ _not_ in the hash table; if addr was in the hash table,
    // this would not have been a cache miss.

    hash_i = hash_addr;
    while( unlikely( voxel_cache_key[hash_i]!=0xffffffff  ) )
      hash_i = ( hash_i + 1 ) & ( VOXEL_CACHE_N_HASH - 1 );
    voxel_cache_key[hash_i] = addr;
    voxel_cache_val[hash_i] = line;
  
  }

  return line;
}

// voxel_cache_wait completes all pending cache operations.  After
// voxel_cache_wait returns, the cache is _guaranteed_ to contain at
// least the previous VOXEL_CACHE_N_LINE fetches of data.  Even more
// specifically, after voxel_cache_wait returns, the cache is
// _guaranteed_ to contain the data for the previous
// VOXEL_CACHE_N_LINE unique addresses fetched.  For example, if you
// have a cache with 128 lines and perform 129 fetches before calling
// voxel_cache_wait, you _are_ guaranteed that the data needed for the
// last 128 fetches is in the cache.  You are _not_ guaranteed that
// the data for the first fetch is in the cache unless the were no
// more than 128 unique addresses in the fetches.

static __inline void
voxel_cache_wait( void ) {
  mfc_write_tag_mask( MASK_BIT_RANGE( uint32_t,
                                      I_CACHE_DMA_CHANNEL_FIRST,
                                      I_CACHE_DMA_CHANNEL_LAST ) |
                      MASK_BIT_RANGE( uint32_t,
                                      A_CACHE_DMA_CHANNEL_FIRST,
                                      A_CACHE_DMA_CHANNEL_LAST ) );
  mfc_read_tag_status_all();
}

// voxel_cache_writeback writes back all the accumulators stored in
// the cache to memory.  FIXME: COULD THIS BE MORE EFFICIENT IN ITS
// USE OF DMA CHANNELS?  AND IS THE PARANOIA NECESSARY?

static __inline int
voxel_cache_writeback( void ) {
  uint32_t addr;
  uint32_t line;
  
  voxel_cache_wait(); // Paranoid?

  for( line=0; line<VOXEL_CACHE_N_LINE; line++ ) {
    addr = voxel_cache_addr[line];
    if( addr==0xffffffff ) continue; // Unlikely??
    mfc_put( &a_cache[line], args->a0+sizeof(accumulator_t)*addr,
             sizeof(accumulator_t), A_CACHE_DMA_CHANNEL(addr), 0, 0 );
  }

  voxel_cache_wait(); // Paranoid?
}

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

static __inline int
move_p_spu( particle_t       * ALIGNED(32) p,
            particle_mover_t * ALIGNED(16) pm ) {
  float s_midx, s_midy, s_midz;
  float s_dispx, s_dispy, s_dispz;
  float s_dir[3];
  float v0, v1, v2, v3, v4, v5;
  int type;
  int64_t neighbor;
  uint32_t line;

  // FIXME: THIS CAN BE PARTIALLY HORIZONTAL SIMD ACCELERATED

  for(;;) {

    line = voxel_cache_fetch( p->i );

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

    voxel_cache_wait();

#   define ACCUMULATE_J(X,Y,Z)                                          \
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
    a_cache[line].j##X[0] += v0;                                        \
    a_cache[line].j##X[1] += v1;                                        \
    a_cache[line].j##X[2] += v2;                                        \
    a_cache[line].j##X[3] += v3

    ACCUMULATE_J(x,y,z);
    ACCUMULATE_J(y,z,x);
    ACCUMULATE_J(z,x,y);

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
    neighbor = i_cache[line].neighbor[ ((v0>0)?3:0) + type ];
    if( unlikely( neighbor<args->rangel |         // YES! BIT OR
                  neighbor>args->rangeh ) ) {     // Hit a boundary
      (&(p->dx))[type] = v0;                      // Put on boundary
      if( neighbor!=reflect_particles ) return 1; // Cannot handle it
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

static __inline int                         // Return number of movers used
advance_p_pipeline_spu( particle_t       * __restrict ALIGNED(128) p,  // Particle array
                        particle_mover_t * __restrict ALIGNED(16)  pm, // Mover array
                        int idx,   // Index of first particle
                        int np ) { // Number of quad particle quads
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

  const vec_uint4  v_i_cache      = VEC_UINT4( i_cache        );
  const vec_uint4  v_a_cache      = VEC_UINT4( a_cache        );

  vec_float4  p0r,  p1r,  p2r,  p3r;                                             vec_float4  p4r,  p5r,  p6r,  p7r;                                             vec_float4  p8r,  p9r, p10r, p11r;                                             vec_float4 p12r, p13r, p14r, p15r;
  vec_float4  p0u,  p1u,  p2u,  p3u;                                             vec_float4  p4u,  p5u,  p6u,  p7u;                                             vec_float4  p8u,  p9u, p10u, p11u;                                             vec_float4 p12u, p13u, p14u, p15u;

  vec_uint4 vp;                                                                  vec_uint4 vp_1;                                                                vec_uint4 vp_2;                                                                vec_uint4   vp_3;

  vec_float4 dx,   dy,   dz;   vec_uint4 i;                                      vec_float4 dx_1, dy_1, dz_1; vec_uint4 i_1;                                    vec_float4 dx_2, dy_2, dz_2; vec_uint4 i_2;                                    vec_float4 dx_3, dy_3, dz_3; vec_uint4 i_3;
  vec_float4 ux,   uy,   uz,   q;                                                vec_float4 ux_1, uy_1, uz_1, q_1;                                              vec_float4 ux_2, uy_2, uz_2, q_2;                                              vec_float4 ux_3, uy_3, uz_3, q_3; 

  vec_float4 ex0,    dexdy,    dexdz,   d2exdydz;                                vec_float4 ex0_1,  dexdy_1,  dexdz_1, d2exdydz_1;                              vec_float4 ex0_2,  dexdy_2,  dexdz_2, d2exdydz_2;                              vec_float4 ex0_3,  dexdy_3,  dexdz_3, d2exdydz_3;
  vec_float4 ey0,    deydz,    deydx,   d2eydzdx;                                vec_float4 ey0_1,  deydz_1,  deydx_1, d2eydzdx_1;                              vec_float4 ey0_2,  deydz_2,  deydx_2, d2eydzdx_2;                              vec_float4 ey0_3,  deydz_3,  deydx_3, d2eydzdx_3;
  vec_float4 ez0,    dezdx,    dezdy,   d2ezdxdy;                                vec_float4 ez0_1,  dezdx_1,  dezdy_1, d2ezdxdy_1;                              vec_float4 ez0_2,  dezdx_2,  dezdy_2, d2ezdxdy_2;                              vec_float4 ez0_3,  dezdx_3,  dezdy_3, d2ezdxdy_3;
  vec_float4 cbx0,   dcbxdx,   cby0,    dcbydy;                                  vec_float4 cbx0_1, dcbxdx_1, cby0_1,  dcbydy_1;                                vec_float4 cbx0_2, dcbxdx_2, cby0_2,  dcbydy_2;                                vec_float4 cbx0_3, dcbxdx_3, cby0_3,  dcbydy_3;
  vec_float4 cbz0,   dcbzdz,   v12,     v13;                                     vec_float4 cbz0_1, dcbzdz_1, v12_1,   v13_1;                                   vec_float4 cbz0_2, dcbzdz_2, v12_2,   v13_2;                                   vec_float4 cbz0_3, dcbzdz_3, v12_3,   v13_3; 

  vec_float4 hax,   hay,   haz;                                                  vec_float4 hax_1, hay_1, haz_1;                                                vec_float4 hax_2, hay_2, haz_2;                                                vec_float4 hax_3, hay_3, haz_3;
  vec_float4 cbx,   cby,   cbz;                                                  vec_float4 cbx_1, cby_1, cbz_1;                                                vec_float4 cbx_2, cby_2, cbz_2;                                                vec_float4 cbx_3, cby_3, cbz_3;
  vec_float4 ux0,   uy0,   uz0;                                                  vec_float4 ux0_1, uy0_1, uz0_1;                                                vec_float4 ux0_2, uy0_2, uz0_2;                                                vec_float4 ux0_3, uy0_3, uz0_3;
  vec_float4 v14,   cbs,   ths;                                                  vec_float4 v14_1, cbs_1, ths_1;                                                vec_float4 v14_2, cbs_2, ths_2;                                                vec_float4 v14_3, cbs_3, ths_3;
  vec_float4 v15,   v16,   v17;                                                  vec_float4 v15_1, v16_1, v17_1;                                                vec_float4 v15_2, v16_2, v17_2;                                                vec_float4 v15_3, v16_3, v17_3;
  vec_float4 wx0,   wy0,   wz0;                                                  vec_float4 wx0_1, wy0_1, wz0_1;                                                vec_float4 wx0_2, wy0_2, wz0_2;                                                vec_float4 wx0_3, wy0_3, wz0_3;
  vec_float4 uxh,   uyh,   uzh;                                                  vec_float4 uxh_1, uyh_1, uzh_1;                                                vec_float4 uxh_2, uyh_2, uzh_2;                                                vec_float4 uxh_3, uyh_3, uzh_3;

  vec_float4 rgamma;                                                             vec_float4 rgamma_1;                                                           vec_float4 rgamma_2;                                                           vec_float4 rgamma_3; 
  vec_float4 ddx,   ddy,   ddz;                                                  vec_float4 ddx_1, ddy_1, ddz_1;                                                vec_float4 ddx_2, ddy_2, ddz_2;                                                vec_float4 ddx_3, ddy_3, ddz_3;
  vec_float4 dxh,   dyh,   dzh;                                                  vec_float4 dxh_1, dyh_1, dzh_1;                                                vec_float4 dxh_2, dyh_2, dzh_2;                                                vec_float4 dxh_3, dyh_3, dzh_3;
  vec_float4 dx1,   dy1,   dz1;                                                  vec_float4 dx1_1, dy1_1, dz1_1;                                                vec_float4 dx1_2, dy1_2, dz1_2;                                                vec_float4 dx1_3, dy1_3, dz1_3;
  vec_uint4  outbnd;                                                             vec_uint4  outbnd_1;                                                           vec_uint4  outbnd_2;                                                           vec_uint4  outbnd_3; 
  uint32_t gather;                                                               uint32_t gather_1;                                                             uint32_t gather_2;                                                             uint32_t gather_3;

  vec_float4 qa,   ccc;                                                          vec_float4 qa_1, ccc_1;                                                        vec_float4 qa_2, ccc_2;                                                        vec_float4 qa_3, ccc_3;
 
  vec_float4 a0x,   a1x,   a2x,   a3x,   a4x;                                    vec_float4 a0x_1, a1x_1, a2x_1, a3x_1, a4x_1;                                  vec_float4 a0x_2, a1x_2, a2x_2, a3x_2, a4x_2;                                  vec_float4 a0x_3, a1x_3, a2x_3, a3x_3, a4x_3;
  vec_float4 a0y,   a1y,   a2y,   a3y,   a4y;                                    vec_float4 a0y_1, a1y_1, a2y_1, a3y_1, a4y_1;                                  vec_float4 a0y_2, a1y_2, a2y_2, a3y_2, a4y_2;                                  vec_float4 a0y_3, a1y_3, a2y_3, a3y_3, a4y_3;
  vec_float4 a0z,   a1z,   a2z,   a3z,   a4z;                                    vec_float4 a0z_1, a1z_1, a2z_1, a3z_1, a4z_1;                                  vec_float4 a0z_2, a1z_2, a2z_2, a3z_2, a4z_2;                                  vec_float4 a0z_3, a1z_3, a2z_3, a3z_3, a4z_3;

  const interpolator_t * __restrict ALIGNED(128) fi;
  accumulator_t        * __restrict ALIGNED(64)  ja;

  int m, n, nm, new_nm;

  DECLARE_ALIGNED_ARRAY( uint32_t, 128, line, NP_BLOCK );

  // Prefetch interpolator and accumulator data for this block

  for( n=0; n<np; n++ ) line[n] = voxel_cache_fetch( p[n].i ) << 6; // FIXME: UGLY!
  voxel_cache_wait();

  // Process the particle quads for this pipeline

  for( n=nm=0; n<np; n+=16 ) {

    // Load the particle quad positions

    LOAD_4x1( &p[n+ 0].dx,  p0r );                                                 LOAD_4x1( &p[n+ 4].dx,  p4r );                                                 LOAD_4x1( &p[n+ 8].dx,  p8r );                                                 LOAD_4x1( &p[n+12].dx, p12r );
    LOAD_4x1( &p[n+ 1].dx,  p1r );                                                 LOAD_4x1( &p[n+ 5].dx,  p5r );                                                 LOAD_4x1( &p[n+ 9].dx,  p9r );                                                 LOAD_4x1( &p[n+13].dx, p13r );
    LOAD_4x1( &p[n +2].dx,  p2r );                                                 LOAD_4x1( &p[n+ 6].dx,  p6r );                                                 LOAD_4x1( &p[n+10].dx, p10r );                                                 LOAD_4x1( &p[n+14].dx, p14r );
    LOAD_4x1( &p[n +3].dx,  p3r );                                                 LOAD_4x1( &p[n+ 7].dx,  p7r );                                                 LOAD_4x1( &p[n+11].dx, p11r );                                                 LOAD_4x1( &p[n+15].dx, p15r );
    TRANSPOSE(  p0r,  p1r,  p2r,  p3r );                                           TRANSPOSE(  p4r,  p5r,  p6r,  p7r );                                           TRANSPOSE(  p8r,  p9r, p10r, p11r );                                           TRANSPOSE( p12r, p13r, p14r, p15r );
    dx   =  p0r;                                                                   dx_1 =  p4r;                                                                   dx_2 =  p8r;                                                                   dx_3 = p12r;
    dy   =  p1r;                                                                   dy_1 =  p5r;                                                                   dy_2 =  p9r;                                                                   dy_3 = p13r;
    dz   =  p2r;                                                                   dz_1 =  p6r;                                                                   dz_2 = p10r;                                                                   dz_3 = p14r;
    LOAD_UINT_4x1( &line[n   ], i   );                                             LOAD_UINT_4x1( &line[n+ 4], i_1 );                                             LOAD_UINT_4x1( &line[n+ 8], i_2 );                                             LOAD_UINT_4x1( &line[n+12], i_3 );

    // Interpolate fields

#   define LOAD_INTERPOLATOR(VP,I,a,b,c,d,e)    \
    fi = (const interpolator_t *)EXTRACT(VP,I); \
    LOAD_4x1( &fi->ex,  a );                    \
    LOAD_4x1( &fi->ey,  b );                    \
    LOAD_4x1( &fi->ez,  c );                    \
    LOAD_4x1( &fi->cbx, d );                    \
    LOAD_4x1( &fi->cbz, e )
     
    vp   = ADD( v_i_cache, ADD(i,  i  ) );                                         vp_1 = ADD( v_i_cache, ADD(i_1,i_1) );                                         vp_2 = ADD( v_i_cache, ADD(i_2,i_2) );                                         vp_3 = ADD( v_i_cache, ADD(i_3,i_3) );

    LOAD_INTERPOLATOR( vp,   0, ex0,        ey0,        ez0,        cbx0,     cbz0     ); 
    LOAD_INTERPOLATOR( vp,   1, dexdy,      deydz,      dezdx,      dcbxdx,   dcbzdz   );
    LOAD_INTERPOLATOR( vp,   2, dexdz,      deydx,      dezdy,      cby0,     v12      );
    LOAD_INTERPOLATOR( vp,   3, d2exdydz,   d2eydzdx,   d2ezdxdy,   dcbydy,   v13      );

    LOAD_INTERPOLATOR( vp_1, 0, ex0_1,      ey0_1,      ez0_1,      cbx0_1,   cbz0_1   );
    LOAD_INTERPOLATOR( vp_1, 1, dexdy_1,    deydz_1,    dezdx_1,    dcbxdx_1, dcbzdz_1 );
    LOAD_INTERPOLATOR( vp_1, 2, dexdz_1,    deydx_1,    dezdy_1,    cby0_1,   v12_1    );
    LOAD_INTERPOLATOR( vp_1, 3, d2exdydz_1, d2eydzdx_1, d2ezdxdy_1, dcbydy_1, v13_1    );

    LOAD_INTERPOLATOR( vp_2, 0, ex0_2,      ey0_2,      ez0_2,      cbx0_2,   cbz0_2   ); 
    LOAD_INTERPOLATOR( vp_2, 1, dexdy_2,    deydz_2,    dezdx_2,    dcbxdx_2, dcbzdz_2 );
    LOAD_INTERPOLATOR( vp_2, 2, dexdz_2,    deydx_2,    dezdy_2,    cby0_2,   v12_2    );
    LOAD_INTERPOLATOR( vp_2, 3, d2exdydz_2, d2eydzdx_2, d2ezdxdy_2, dcbydy_2, v13_2    );

    LOAD_INTERPOLATOR( vp_3, 0, ex0_3,      ey0_3,      ez0_3,      cbx0_3,   cbz0_3   );
    LOAD_INTERPOLATOR( vp_3, 1, dexdy_3,    deydz_3,    dezdx_3,    dcbxdx_3, dcbzdz_3 );
    LOAD_INTERPOLATOR( vp_3, 2, dexdz_3,    deydx_3,    dezdy_3,    cby0_3,   v12_3    );
    LOAD_INTERPOLATOR( vp_3, 3, d2exdydz_3, d2eydzdx_3, d2ezdxdy_3, dcbydy_3, v13_3    );

    TRANSPOSE( ex0,   dexdy,   dexdz,   d2exdydz   );                              TRANSPOSE( ex0_1, dexdy_1, dexdz_1, d2exdydz_1 );                              TRANSPOSE( ex0_2, dexdy_2, dexdz_2, d2exdydz_2 );                              TRANSPOSE( ex0_3, dexdy_3, dexdz_3, d2exdydz_3 );
    hax   = MUL( qdt_2mc, FMA( FMA( d2exdydz,   dy,   dexdz   ), dz,   FMA( dexdy,   dy,   ex0   ) ) );
    hax_1 = MUL( qdt_2mc, FMA( FMA( d2exdydz_1, dy_1, dexdz_1 ), dz_1, FMA( dexdy_1, dy_1, ex0_1 ) ) );
    hax_2 = MUL( qdt_2mc, FMA( FMA( d2exdydz_2, dy_2, dexdz_2 ), dz_2, FMA( dexdy_2, dy_2, ex0_2 ) ) );
    hax_3 = MUL( qdt_2mc, FMA( FMA( d2exdydz_3, dy_3, dexdz_3 ), dz_3, FMA( dexdy_3, dy_3, ex0_3 ) ) );

    TRANSPOSE( ey0,   deydz,   deydx,   d2eydzdx   );                              TRANSPOSE( ey0_1, deydz_1, deydx_1, d2eydzdx_1 );                              TRANSPOSE( ey0_2, deydz_2, deydx_2, d2eydzdx_2 );                              TRANSPOSE( ey0_3, deydz_3, deydx_3, d2eydzdx_3 );
    hay   = MUL( qdt_2mc, FMA( FMA( d2eydzdx,   dz,   deydx   ), dx,   FMA( deydz,   dz,   ey0   ) ) );
    hay_1 = MUL( qdt_2mc, FMA( FMA( d2eydzdx_1, dz_1, deydx_1 ), dx_1, FMA( deydz_1, dz_1, ey0_1 ) ) );
    hay_2 = MUL( qdt_2mc, FMA( FMA( d2eydzdx_2, dz_2, deydx_2 ), dx_2, FMA( deydz_2, dz_2, ey0_2 ) ) );
    hay_3 = MUL( qdt_2mc, FMA( FMA( d2eydzdx_3, dz_3, deydx_3 ), dx_3, FMA( deydz_3, dz_3, ey0_3 ) ) );

    TRANSPOSE( ez0,   dezdx,   dezdy,   d2ezdxdy   );                              TRANSPOSE( ez0_1, dezdx_1, dezdy_1, d2ezdxdy_1 );                              TRANSPOSE( ez0_2, dezdx_2, dezdy_2, d2ezdxdy_2 );                              TRANSPOSE( ez0_3, dezdx_3, dezdy_3, d2ezdxdy_3 );
    haz   = MUL( qdt_2mc, FMA( FMA( d2ezdxdy,   dx,   dezdy   ), dy,   FMA( dezdx,   dx,   ez0   ) ) );
    haz_1 = MUL( qdt_2mc, FMA( FMA( d2ezdxdy_1, dx_1, dezdy_1 ), dy_1, FMA( dezdx_1, dx_1, ez0_1 ) ) );
    haz_2 = MUL( qdt_2mc, FMA( FMA( d2ezdxdy_2, dx_2, dezdy_2 ), dy_2, FMA( dezdx_2, dx_2, ez0_2 ) ) );
    haz_3 = MUL( qdt_2mc, FMA( FMA( d2ezdxdy_3, dx_3, dezdy_3 ), dy_3, FMA( dezdx_3, dx_3, ez0_3 ) ) );

    TRANSPOSE( cbx0,   dcbxdx,   cby0,   dcbydy   );                               TRANSPOSE( cbx0_1, dcbxdx_1, cby0_1, dcbydy_1 );                               TRANSPOSE( cbx0_2, dcbxdx_2, cby0_2, dcbydy_2 );                               TRANSPOSE( cbx0_3, dcbxdx_3, cby0_3, dcbydy_3 );
    cbx   = FMA( dcbxdx,   dx,   cbx0   );                                         cbx_1 = FMA( dcbxdx_1, dx_1, cbx0_1 );                                         cbx_2 = FMA( dcbxdx_2, dx_2, cbx0_2 );                                         cbx_3 = FMA( dcbxdx_3, dx_3, cbx0_3 );
    cby   = FMA( dcbydy,   dy,   cby0   );                                         cby_1 = FMA( dcbydy_1, dy_1, cby0_1 );                                         cby_2 = FMA( dcbydy_2, dy_2, cby0_2 );                                         cby_3 = FMA( dcbydy_3, dy_3, cby0_3 );

    HALF_TRANSPOSE( cbz0,   dcbzdz,   v12,   v13   );                              HALF_TRANSPOSE( cbz0_1, dcbzdz_1, v12_1, v13_1 );                              HALF_TRANSPOSE( cbz0_2, dcbzdz_2, v12_2, v13_2 );                              HALF_TRANSPOSE( cbz0_3, dcbzdz_3, v12_3, v13_3 );
    cbz   = FMA( dcbzdz,   dz,   cbz0   );                                         cbz_1 = FMA( dcbzdz_1, dz_1, cbz0_1 );                                         cbz_2 = FMA( dcbzdz_2, dz_2, cbz0_2 );                                         cbz_3 = FMA( dcbzdz_3, dz_3, cbz0_3 );

    // Update momentum.  Note: Could eliminate a dependency in v14 calc
    // if willing to play fast and loose with numerics (saves about a spu
    // clock per particle).

    LOAD_4x1( &p[n+ 0].ux,  p0u );                                                 LOAD_4x1( &p[n+ 4].ux,  p4u );                                                 LOAD_4x1( &p[n+ 8].ux,  p8u );                                                 LOAD_4x1( &p[n+12].ux, p12u );
    LOAD_4x1( &p[n+ 1].ux,  p1u );                                                 LOAD_4x1( &p[n+ 5].ux,  p5u );                                                 LOAD_4x1( &p[n+ 9].ux,  p9u );                                                 LOAD_4x1( &p[n+13].ux, p13u );
    LOAD_4x1( &p[n+ 2].ux,  p2u );                                                 LOAD_4x1( &p[n+ 6].ux,  p6u );                                                 LOAD_4x1( &p[n+10].ux, p10u );                                                 LOAD_4x1( &p[n+14].ux, p14u );
    LOAD_4x1( &p[n+ 3].ux,  p3u );                                                 LOAD_4x1( &p[n+ 7].ux,  p7u );                                                 LOAD_4x1( &p[n+11].ux, p11u );                                                 LOAD_4x1( &p[n+15].ux, p15u );
    TRANSPOSE(  p0u,  p1u,  p2u,  p3u );                                           TRANSPOSE(  p4u,  p5u,  p6u,  p7u );                                           TRANSPOSE(  p8u,  p9u, p10u, p11u );                                           TRANSPOSE( p12u, p13u, p14u, p15u );
    ux    = p0u;                                                                   ux_1  = p4u;                                                                   ux_2  = p8u;                                                                   ux_3  = p12u;
    uy    = p1u;                                                                   uy_1  = p5u;                                                                   uy_2  = p9u;                                                                   uy_3  = p13u;
    uz    = p2u;                                                                   uz_1  = p6u;                                                                   uz_2  = p10u;                                                                  uz_3  = p14u;
    q     = p3u;                                                                   q_1   = p7u;                                                                   q_2   = p11u;                                                                  q_3   = p15u;
    ux0   = ADD( ux,   hax   );                                                    ux0_1 = ADD( ux_1, hax_1 );                                                    ux0_2 = ADD( ux_2, hax_2 );                                                    ux0_3 = ADD( ux_3, hax_3 );
    uy0   = ADD( uy,   hay   );                                                    uy0_1 = ADD( uy_1, hay_1 );                                                    uy0_2 = ADD( uy_2, hay_2 );                                                    uy0_3 = ADD( uy_3, hay_3 );
    uz0   = ADD( uz,   haz   );                                                    uz0_1 = ADD( uz_1, haz_1 );                                                    uz0_2 = ADD( uz_2, haz_2 );                                                    uz0_3 = ADD( uz_3, haz_3 );
    v14   = MUL( qdt_2mc, RSQRT( ADD( one, FMA( ux0,  ux0,   FMA( uy0,  uy0,   MUL( uz0,  uz0   ) ) ) ) ) );
    v14_1 = MUL( qdt_2mc, RSQRT( ADD( one, FMA( ux0_1,ux0_1, FMA( uy0_1,uy0_1, MUL( uz0_1,uz0_1 ) ) ) ) ) );
    v14_2 = MUL( qdt_2mc, RSQRT( ADD( one, FMA( ux0_2,ux0_2, FMA( uy0_2,uy0_2, MUL( uz0_2,uz0_2 ) ) ) ) ) );
    v14_3 = MUL( qdt_2mc, RSQRT( ADD( one, FMA( ux0_3,ux0_3, FMA( uy0_3,uy0_3, MUL( uz0_3,uz0_3 ) ) ) ) ) );
    cbs   = FMA( cbx,  cbx,   FMA( cby,  cby,   MUL(cbz,  cbz  ) ) );              cbs_1 = FMA( cbx_1,cbx_1, FMA( cby_1,cby_1, MUL(cbz_1,cbz_1) ) );              cbs_2 = FMA( cbx_2,cbx_2, FMA( cby_2,cby_2, MUL(cbz_2,cbz_2) ) );              cbs_3 = FMA( cbx_3,cbx_3, FMA( cby_3,cby_3, MUL(cbz_3,cbz_3) ) ); 
    ths   = MUL( MUL( v14,  v14   ), cbs   );                                      ths_1 = MUL( MUL( v14_1,v14_1 ), cbs_1 );                                      ths_2 = MUL( MUL( v14_2,v14_2 ), cbs_2 );                                      ths_3 = MUL( MUL( v14_3,v14_3 ), cbs_3 );
    v15   = MUL(v14,  FMA(FMA(two_fifteenths,ths,  one_third ),ths,  one));        v15_1 = MUL(v14_1,FMA(FMA(two_fifteenths,ths_1,one_third ),ths_1,one));        v15_2 = MUL(v14_2,FMA(FMA(two_fifteenths,ths_2,one_third ),ths_2,one));        v15_3 = MUL(v14_3,FMA(FMA(two_fifteenths,ths_3,one_third ),ths_3,one));
    v16   = MUL( v15,   RCP( FMA( MUL(v15,  v15  ), cbs,   one ) ) );              v16_1 = MUL( v15_1, RCP( FMA( MUL(v15_1,v15_1), cbs_1, one ) ) );              v16_2 = MUL( v15_2, RCP( FMA( MUL(v15_2,v15_2), cbs_2, one ) ) );              v16_3 = MUL( v15_3, RCP( FMA( MUL(v15_3,v15_3), cbs_3, one ) ) );
    v17   = ADD( v16,   v16   );                                                   v17_1 = ADD( v16_1, v16_1 );                                                   v17_2 = ADD( v16_2, v16_2 );                                                   v17_3 = ADD( v16_3, v16_3 );
    wx0   =     FMA(FMS(uy0,  cbz,   MUL(uz0,  cby  )), v15,   ux0  );             wx0_1 =     FMA(FMS(uy0_1,cbz_1, MUL(uz0_1,cby_1)), v15_1, ux0_1);             wx0_2 =     FMA(FMS(uy0_2,cbz_2, MUL(uz0_2,cby_2)), v15_2, ux0_2);             wx0_3 =     FMA(FMS(uy0_3,cbz_3, MUL(uz0_3,cby_3)), v15_3, ux0_3);
    wy0   =     FMA(FMS(uz0,  cbx,   MUL(ux0,  cbz  )), v15,   uy0  );             wy0_1 =     FMA(FMS(uz0_1,cbx_1, MUL(ux0_1,cbz_1)), v15_1, uy0_1);             wy0_2 =     FMA(FMS(uz0_2,cbx_2, MUL(ux0_2,cbz_2)), v15_2, uy0_2);             wy0_3 =     FMA(FMS(uz0_3,cbx_3, MUL(ux0_3,cbz_3)), v15_3, uy0_3);
    wz0   =     FMA(FMS(ux0,  cby,   MUL(uy0,  cbx  )), v15,   uz0  );             wz0_1 =     FMA(FMS(ux0_1,cby_1, MUL(uy0_1,cbx_1)), v15_1, uz0_1);             wz0_2 =     FMA(FMS(ux0_2,cby_2, MUL(uy0_2,cbx_2)), v15_2, uz0_2);             wz0_3 =     FMA(FMS(ux0_3,cby_3, MUL(uy0_3,cbx_3)), v15_3, uz0_3);
    uxh   = ADD(FMA(FMS(wy0,  cbz,   MUL(wz0,  cby  )), v17,   ux0  ),hax  );      uxh_1 = ADD(FMA(FMS(wy0_1,cbz_1, MUL(wz0_1,cby_1)), v17_1, ux0_1),hax_1);      uxh_2 = ADD(FMA(FMS(wy0_2,cbz_2, MUL(wz0_2,cby_2)), v17_2, ux0_2),hax_2);      uxh_3 = ADD(FMA(FMS(wy0_3,cbz_3, MUL(wz0_3,cby_3)), v17_3, ux0_3),hax_3);
    uyh   = ADD(FMA(FMS(wz0,  cbx,   MUL(wx0,  cbz  )), v17,   uy0  ),hay  );      uyh_1 = ADD(FMA(FMS(wz0_1,cbx_1, MUL(wx0_1,cbz_1)), v17_1, uy0_1),hay_1);      uyh_2 = ADD(FMA(FMS(wz0_2,cbx_2, MUL(wx0_2,cbz_2)), v17_2, uy0_2),hay_2);      uyh_3 = ADD(FMA(FMS(wz0_3,cbx_3, MUL(wx0_3,cbz_3)), v17_3, uy0_3),hay_3);
    uzh   = ADD(FMA(FMS(wx0,  cby,   MUL(wy0,  cbx  )), v17,   uz0  ),haz  );      uzh_1 = ADD(FMA(FMS(wx0_1,cby_1, MUL(wy0_1,cbx_1)), v17_1, uz0_1),haz_1);      uzh_2 = ADD(FMA(FMS(wx0_2,cby_2, MUL(wy0_2,cbx_2)), v17_2, uz0_2),haz_2);      uzh_3 = ADD(FMA(FMS(wx0_3,cby_3, MUL(wy0_3,cbx_3)), v17_3, uz0_3),haz_3);
    p0u  = uxh;                                                                    p4u  = uxh_1;                                                                  p8u  = uxh_2;                                                                  p12u = uxh_3;
    p1u  = uyh;                                                                    p5u  = uyh_1;                                                                  p9u  = uyh_2;                                                                  p13u = uyh_3;
    p2u  = uzh;                                                                    p6u  = uzh_1;                                                                  p10u = uzh_2;                                                                  p14u = uzh_3;
    /* p3u  is unchanged */                                                        /* p7u  is unchanged */                                                        /* p11u is unchanged */                                                        /* p15u is unchanged */
    TRANSPOSE(  p0u,  p1u,  p2u,  p3u );                                           TRANSPOSE(  p4u,  p5u,  p6u,  p7u );                                           TRANSPOSE(  p8u,  p9u, p10u, p11u );                                           TRANSPOSE( p12u, p13u, p14u, p15u );
    STORE_4x1(  p0u, &p[n+ 0].ux );                                                STORE_4x1(  p4u, &p[n+ 4].ux );                                                STORE_4x1(  p8u, &p[n+ 8].ux );                                                STORE_4x1( p12u, &p[n+12].ux );
    STORE_4x1(  p1u, &p[n+ 1].ux );                                                STORE_4x1(  p5u, &p[n+ 5].ux );                                                STORE_4x1(  p9u, &p[n+ 9].ux );                                                STORE_4x1( p13u, &p[n+13].ux );
    STORE_4x1(  p2u, &p[n+ 2].ux );                                                STORE_4x1(  p6u, &p[n+ 6].ux );                                                STORE_4x1( p10u, &p[n+10].ux );                                                STORE_4x1( p14u, &p[n+14].ux );
    STORE_4x1(  p3u, &p[n+ 3].ux );                                                STORE_4x1(  p7u, &p[n+ 7].ux );                                                STORE_4x1( p11u, &p[n+11].ux );                                                STORE_4x1( p15u, &p[n+15].ux );
    
    // Update the position of inbnd particles

    rgamma   = RSQRT( ADD( one, FMA( uxh,  uxh,   FMA( uyh,  uyh,   MUL(uzh,  uzh  ) ) ) ) );
    rgamma_1 = RSQRT( ADD( one, FMA( uxh_1,uxh_1, FMA( uyh_1,uyh_1, MUL(uzh_1,uzh_1) ) ) ) );
    rgamma_2 = RSQRT( ADD( one, FMA( uxh_2,uxh_2, FMA( uyh_2,uyh_2, MUL(uzh_2,uzh_2) ) ) ) );
    rgamma_3 = RSQRT( ADD( one, FMA( uxh_3,uxh_3, FMA( uyh_3,uyh_3, MUL(uzh_3,uzh_3) ) ) ) );
    ddx      = MUL( MUL( uxh,   cdt_dx ), rgamma   );                              ddx_1    = MUL( MUL( uxh_1, cdt_dx ), rgamma_1 );                              ddx_2    = MUL( MUL( uxh_2, cdt_dx ), rgamma_2 );                              ddx_3    = MUL( MUL( uxh_3, cdt_dx ), rgamma_3 );
    ddy      = MUL( MUL( uyh,   cdt_dy ), rgamma   );                              ddy_1    = MUL( MUL( uyh_1, cdt_dy ), rgamma_1 );                              ddy_2    = MUL( MUL( uyh_2, cdt_dy ), rgamma_2 );                              ddy_3    = MUL( MUL( uyh_3, cdt_dy ), rgamma_3 );
    ddz      = MUL( MUL( uzh,   cdt_dz ), rgamma   );                              ddz_1    = MUL( MUL( uzh_1, cdt_dz ), rgamma_1 );                              ddz_2    = MUL( MUL( uzh_2, cdt_dz ), rgamma_2 );                              ddz_3    = MUL( MUL( uzh_3, cdt_dz ), rgamma_3 );
    dxh      = ADD( dx,    ddx   );                                                dxh_1    = ADD( dx_1,  ddx_1 );                                                dxh_2    = ADD( dx_2,  ddx_2 );                                                dxh_3    = ADD( dx_3,  ddx_3 );
    dyh      = ADD( dy,    ddy   );                                                dyh_1    = ADD( dy_1,  ddy_1 );                                                dyh_2    = ADD( dy_2,  ddy_2 );                                                dyh_3    = ADD( dy_3,  ddy_3 );
    dzh      = ADD( dz,    ddz   );                                                dzh_1    = ADD( dz_1,  ddz_1 );                                                dzh_2    = ADD( dz_2,  ddz_2 );                                                dzh_3    = ADD( dz_3,  ddz_3 );
    dx1      = ADD( dxh,   ddx   );                                                dx1_1    = ADD( dxh_1, ddx_1 );                                                dx1_2    = ADD( dxh_2, ddx_2 );                                                dx1_3    = ADD( dxh_3, ddx_3 );
    dy1      = ADD( dyh,   ddy   );                                                dy1_1    = ADD( dyh_1, ddy_1 );                                                dy1_2    = ADD( dyh_2, ddy_2 );                                                dy1_3    = ADD( dyh_3, ddy_3 );
    dz1      = ADD( dzh,   ddz   );                                                dz1_1    = ADD( dzh_1, ddz_1 );                                                dz1_2    = ADD( dzh_2, ddz_2 );                                                dz1_3    = ADD( dzh_3, ddz_3 );
    outbnd   = OR( OR( OR( CMPLT(dx1,  neg_one), CMPGT(dx1,  one) ), OR( CMPLT(dy1,  neg_one), CMPGT(dy1,  one) ) ), OR( CMPLT(dz1,  neg_one), CMPGT(dz1,  one) ) );
    outbnd_1 = OR( OR( OR( CMPLT(dx1_1,neg_one), CMPGT(dx1_1,one) ), OR( CMPLT(dy1_1,neg_one), CMPGT(dy1_1,one) ) ), OR( CMPLT(dz1_1,neg_one), CMPGT(dz1_1,one) ) );
    outbnd_2 = OR( OR( OR( CMPLT(dx1_2,neg_one), CMPGT(dx1_2,one) ), OR( CMPLT(dy1_2,neg_one), CMPGT(dy1_2,one) ) ), OR( CMPLT(dz1_2,neg_one), CMPGT(dz1_2,one) ) );
    outbnd_3 = OR( OR( OR( CMPLT(dx1_3,neg_one), CMPGT(dx1_3,one) ), OR( CMPLT(dy1_3,neg_one), CMPGT(dy1_3,one) ) ), OR( CMPLT(dz1_3,neg_one), CMPGT(dz1_3,one) ) );

#   define TEST_OUTBND(Q,E,o)                       \
    if( unlikely( ( gather##Q & (1<<E) ) != 0 ) ) { \
       pm[nm].dispx = EXTRACT( ddx##Q, E );         \
       pm[nm].dispx = EXTRACT( ddy##Q, E );         \
       pm[nm].dispx = EXTRACT( ddz##Q, E );         \
       pm[nm].i     = n + o;                        \
       nm++;                                        \
     }

    gather   = GATHER( outbnd   );                                                 gather_1 = GATHER( outbnd_1 );                                                 gather_2 = GATHER( outbnd_2 );                                                 gather_3 = GATHER( outbnd_3 );

    if( unlikely( gather  !=0 ) ) { TEST_OUTBND(,  0, 0); TEST_OUTBND(,  1, 1); TEST_OUTBND(,  2, 2); TEST_OUTBND(,  3, 3); }
    if( unlikely( gather_1!=0 ) ) { TEST_OUTBND(_1,0, 4); TEST_OUTBND(_1,1, 5); TEST_OUTBND(_1,2, 6); TEST_OUTBND(_1,3, 7); }
    if( unlikely( gather_2!=0 ) ) { TEST_OUTBND(_2,0, 8); TEST_OUTBND(_2,1, 9); TEST_OUTBND(_2,2,10); TEST_OUTBND(_2,3,11); }
    if( unlikely( gather_3!=0 ) ) { TEST_OUTBND(_3,0,12); TEST_OUTBND(_3,1,13); TEST_OUTBND(_3,2,14); TEST_OUTBND(_3,3,14); }

#   undef TEST_OUTBND

    p0r      = MERGE(outbnd,  dx,  dx1  );                                         p4r      = MERGE(outbnd_1,dx_1,dx1_1);                                         p8r      = MERGE(outbnd_2,dx_2,dx1_2);                                         p12r     = MERGE(outbnd_3,dx_3,dx1_3);
    p1r      = MERGE(outbnd,  dy,  dy1  );                                         p5r      = MERGE(outbnd_1,dy_1,dy1_1);                                         p9r      = MERGE(outbnd_2,dy_2,dy1_2);                                         p13r     = MERGE(outbnd_3,dy_3,dy1_3);
    p2r      = MERGE(outbnd,  dz,  dz1  );                                         p6r      = MERGE(outbnd_1,dz_1,dz1_1);                                         p10r     = MERGE(outbnd_2,dz_2,dz1_2);                                         p14r     = MERGE(outbnd_3,dz_3,dz1_3);
    /* p3r  is unchanged */                                                        /* p7r  is unchanged */                                                        /* p11r is unchanged */                                                        /* p15r is unchanged */
    TRANSPOSE(  p0r,  p1r,  p2r,  p3r );                                           TRANSPOSE(  p4r,  p5r,  p6r,  p7r );                                           TRANSPOSE(  p8r,  p9r, p10r, p11r );                                           TRANSPOSE( p12r, p13r, p14r, p15r );
    STORE_4x1(  p0r, &p[n+ 0].dx );                                                STORE_4x1(  p4r, &p[n+ 4].dx );                                                STORE_4x1(  p8r, &p[n+ 8].dx );                                                STORE_4x1( p12r, &p[n+12].dx );
    STORE_4x1(  p1r, &p[n+ 1].dx );                                                STORE_4x1(  p5r, &p[n+ 5].dx );                                                STORE_4x1(  p9r, &p[n+ 9].dx );                                                STORE_4x1( p13r, &p[n+13].dx );
    STORE_4x1(  p2r, &p[n+ 2].dx );                                                STORE_4x1(  p6r, &p[n+ 6].dx );                                                STORE_4x1( p10r, &p[n+10].dx );                                                STORE_4x1( p14r, &p[n+14].dx );
    STORE_4x1(  p3r, &p[n+ 3].dx );                                                STORE_4x1(  p7r, &p[n+ 7].dx );                                                STORE_4x1( p11r, &p[n+11].dx );                                                STORE_4x1( p15r, &p[n+15].dx );

    // Accumulate current of inbnd particles
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step

    qa    = CZERO( outbnd,   q   );                                                qa_1  = CZERO( outbnd_1, q_1 );                                                qa_2  = CZERO( outbnd_2, q_2 );                                                qa_3  = CZERO( outbnd_3, q_3 );
    ccc   = MUL( MUL( one_third, qa   ), MUL( ddx,   MUL( ddy,   ddz   ) ) );      ccc_1 = MUL( MUL( one_third, qa_1 ), MUL( ddx_1, MUL( ddy_1, ddz_1 ) ) );
    ccc_2 = MUL( MUL( one_third, qa_2 ), MUL( ddx_2, MUL( ddy_2, ddz_2 ) ) );      ccc_3 = MUL( MUL( one_third, qa_3 ), MUL( ddx_3, MUL( ddy_3, ddz_3 ) ) );

#   define ACCUMULATE_J(X,Y,Z)                                                                                                                                                                                                                                                         \
    a4##X     = MUL(qa,       dd##X    );                                          a4##X##_1 = MUL(qa_1,     dd##X##_1);                                          a4##X##_2 = MUL(qa_2,     dd##X##_2);                                          a4##X##_3 = MUL(qa_3,     dd##X##_3); \
    a1##X     = MUL(a4##X,    d##Y##h  );                                          a1##X##_1 = MUL(a4##X##_1,d##Y##h_1);                                          a1##X##_2 = MUL(a4##X##_2,d##Y##h_2);                                          a1##X##_3 = MUL(a4##X##_3,d##Y##h_3); \
    a0##X     = SUB(a4##X,    a1##X    );                                          a0##X##_1 = SUB(a4##X##_1,a1##X##_1);                                          a0##X##_2 = SUB(a4##X##_2,a1##X##_2);                                          a0##X##_3 = SUB(a4##X##_3,a1##X##_3); \
    a1##X     = ADD(a1##X,    a4##X    );                                          a1##X##_1 = ADD(a1##X##_1,a4##X##_1);                                          a1##X##_2 = ADD(a1##X##_2,a4##X##_2);                                          a1##X##_3 = ADD(a1##X##_3,a4##X##_3); \
    a4##X     = ADD(one,      d##Z##h  );                                          a4##X##_1 = ADD(one,      d##Z##h_1);                                          a4##X##_2 = ADD(one,      d##Z##h_2);                                          a4##X##_3 = ADD(one,      d##Z##h_3); \
    a2##X     = MUL(a0##X,    a4##X    );                                          a2##X##_1 = MUL(a0##X##_1,a4##X##_1);                                          a2##X##_2 = MUL(a0##X##_2,a4##X##_2);                                          a2##X##_3 = MUL(a0##X##_3,a4##X##_3); \
    a3##X     = MUL(a1##X,    a4##X    );                                          a3##X##_1 = MUL(a1##X##_1,a4##X##_1);                                          a3##X##_2 = MUL(a1##X##_2,a4##X##_2);                                          a3##X##_3 = MUL(a1##X##_3,a4##X##_3); \
    a4##X     = SUB(one,      d##Z##h  );                                          a4##X##_1 = SUB(one,      d##Z##h_1);                                          a4##X##_2 = SUB(one,      d##Z##h_2);                                          a4##X##_3 = SUB(one,      d##Z##h_3); \
    a0##X     = MUL(a0##X,    a4##X    );                                          a0##X##_1 = MUL(a0##X##_1,a4##X##_1);                                          a0##X##_2 = MUL(a0##X##_2,a4##X##_2);                                          a0##X##_3 = MUL(a0##X##_3,a4##X##_3); \
    a1##X     = MUL(a1##X,    a4##X    );                                          a1##X##_1 = MUL(a1##X##_1,a4##X##_1);                                          a1##X##_2 = MUL(a1##X##_2,a4##X##_2);                                          a1##X##_3 = MUL(a1##X##_3,a4##X##_3); \
    a0##X     = ADD(a0##X,    ccc      );                                          a0##X##_1 = ADD(a0##X##_1,ccc_1    );                                          a0##X##_2 = ADD(a0##X##_2,ccc_2    );                                          a0##X##_3 = ADD(a0##X##_3,ccc_3    ); \
    a1##X     = SUB(a1##X,    ccc      );                                          a1##X##_1 = SUB(a1##X##_1,ccc_1    );                                          a1##X##_2 = SUB(a1##X##_2,ccc_2    );                                          a1##X##_3 = SUB(a1##X##_3,ccc_3    ); \
    a2##X     = SUB(a2##X,    ccc      );                                          a2##X##_1 = SUB(a2##X##_1,ccc_1    );                                          a2##X##_2 = SUB(a2##X##_2,ccc_2    );                                          a2##X##_3 = SUB(a2##X##_3,ccc_3    ); \
    a3##X     = ADD(a3##X,    ccc      );                                          a3##X##_1 = ADD(a3##X##_1,ccc_1    );                                          a3##X##_2 = ADD(a3##X##_2,ccc_2    );                                          a3##X##_3 = ADD(a3##X##_3,ccc_3    ); \
    TRANSPOSE(a0##X,    a1##X,    a2##X,    a3##X    );                            TRANSPOSE(a0##X##_1,a1##X##_1,a2##X##_1,a3##X##_1);                            TRANSPOSE(a0##X##_2,a1##X##_2,a2##X##_2,a3##X##_2);                            TRANSPOSE(a0##X##_3,a1##X##_3,a2##X##_3,a3##X##_3)

    ACCUMULATE_J( x,y,z );
    ACCUMULATE_J( y,z,x );
    ACCUMULATE_J( z,x,y );

#   undef ACCUMULATE_J

#   define INCREMENT_ACCUMULATOR(VP,I,x,y,z) \
    ja = (accumulator_t *)EXTRACT(VP,I);     \
    INCREMENT_4x1( &ja->jx[0], x );          \
    INCREMENT_4x1( &ja->jy[0], y );          \
    INCREMENT_4x1( &ja->jz[0], z )

    vp   = ADD( v_a_cache, i   );                                                  vp_1 = ADD( v_a_cache, i_1 );                                                  vp_2 = ADD( v_a_cache, i_2 );                                                  vp_3 = ADD( v_a_cache, i_3 );

    INCREMENT_ACCUMULATOR( vp,   0,  a0x,    a0y,    a0z   );
    INCREMENT_ACCUMULATOR( vp,   1,  a1x,    a1y,    a1z   );
    INCREMENT_ACCUMULATOR( vp,   2,  a2x,    a2y,    a2z   );
    INCREMENT_ACCUMULATOR( vp,   3,  a3x,    a3y,    a3z   );

    INCREMENT_ACCUMULATOR( vp_1, 0,  a0x_1,  a0y_1,  a0z_1 );
    INCREMENT_ACCUMULATOR( vp_1, 1,  a1x_1,  a1y_1,  a1z_1 );
    INCREMENT_ACCUMULATOR( vp_1, 2,  a2x_1,  a2y_1,  a2z_1 );
    INCREMENT_ACCUMULATOR( vp_1, 3,  a3x_1,  a3y_1,  a3z_1 );

    INCREMENT_ACCUMULATOR( vp_2, 0,  a0x_2,  a0y_2,  a0z_2 );
    INCREMENT_ACCUMULATOR( vp_2, 1,  a1x_2,  a1y_2,  a1z_2 );
    INCREMENT_ACCUMULATOR( vp_2, 2,  a2x_2,  a2y_2,  a2z_2 );
    INCREMENT_ACCUMULATOR( vp_2, 3,  a3x_2,  a3y_2,  a3z_2 );

    INCREMENT_ACCUMULATOR( vp_3, 0,  a0x_3,  a0y_3,  a0z_3 );
    INCREMENT_ACCUMULATOR( vp_3, 1,  a1x_3,  a1y_3,  a1z_3 );
    INCREMENT_ACCUMULATOR( vp_3, 2,  a2x_3,  a2y_3,  a2z_3 );
    INCREMENT_ACCUMULATOR( vp_3, 3,  a3x_3,  a3y_3,  a3z_3 ); 

  }

  // Process movers

  new_nm = 0;

  for( m=0; m<nm; m++ ) {
    n = pm[m].i;
    if( unlikely( move_p_spu( &p[n], &pm[m] ) ) ) {
      pm[new_nm++]   = pm[m];
      pm[new_nm-1].i = idx + n;
    }
  }

  return new_nm;
}

///////////////////////////////////////////////////////////////////////////////
// main (workload distribution and data buffering)

// FIXME: Some util functionality is not compiled for the spu

int
main( uint64_t spu_id,
      uint64_t argp,
      uint64_t envp ) {
  particle_t       * ALIGNED(128) p_block[3];
  int idx[3];
  int np_block[3];

  particle_mover_t * ALIGNED(128) m_block[3];
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
    np_block[buffer] = NP_BLOCK;                                \
    if( unlikely( np_block[buffer]>np ) ) np_block[buffer] = np; \
    next_idx += np_block[buffer];                               \
    np       -= np_block[buffer];                               \
                                                                \
    /* If we have a block, start reading it into the buffer */  \
    if( likely( np_block[buffer]!=0 ) )                         \
      mfc_get( p_block[buffer],                                 \
               args->p0 + idx[buffer]*sizeof(particle_t),       \
               np_block[buffer]*sizeof(particle_t),             \
               3*(buffer)+0, 0, 0 );                            \
  } while(0)
  
# define END_GET_PBLOCK(buffer) do {                                    \
    /* If we have a block, stop reading it into the buffer */           \
    if( likely( np_block[buffer]!=0 ) ) {                               \
      mfc_write_tag_mask( 1<<(3*(buffer)+0) );                          \
      mfc_read_tag_status_all();                                        \
    }                                                                   \
  } while(0)
  
# define PROCESS_PBLOCK(buffer)                                         \
  nm_block[buffer] = advance_p_pipeline_spu( p_block[buffer],           \
                                             m_block[buffer],           \
                                             idx[buffer],               \
                                             np_block[buffer] )

  // FIXME: mfc list for this??
# define BEGIN_PUT_PBLOCK(buffer) do {                                  \
                                                                        \
    /* If we have a block, begin writing the buffer into the block */   \
    if( likely( np_block[buffer]!=0 ) )                                 \
      mfc_put( p_block[buffer],                                         \
               args->p0 + idx[buffer]*sizeof(particle_t),               \
               np_block[buffer]*sizeof(particle_t),                     \
               3*(buffer)+1, 0, 0 );                                    \
                                                                        \
    /* Begin writing the movers corresponding to this block */          \
    /* Ignore movers that would overflow the mover segment */           \
    itmp = seg->nm + nm_block[buffer] - seg->max_nm;                    \
    if( unlikely( itmp > 0 ) ) {                                        \
      seg->n_ignored += itmp;                                           \
      nm_block[buffer] = seg->max_nm - seg->nm;                         \
    }                                                                   \
    if( likely( nm_block[buffer]!=0 ) ) { /* Unlikely in cold bench */  \
      mfc_put( m_block[buffer],                                         \
               seg->pm + seg->nm*sizeof(particle_mover_t),              \
               nm_block[buffer]*sizeof(particle_mover_t),               \
               3*(buffer)+2, 0, 0 );                                    \
      seg->nm += nm_block[buffer];                                      \
    }                                                                   \
                                                                        \
  } while(0)

# define END_PUT_PBLOCK(buffer) do {                                    \
    uint32_t mask = 0;                                                  \
    if( likely( np_block[buffer]!=0 ) ) mask |= 1 << (3*(buffer)+1);    \
    if( likely( nm_block[buffer]!=0 ) ) mask |= 1 << (3*(buffer)+2);    \
    /* above is unlikely on cold bench */                               \
    if( likely( mask!=0 ) ) {                                           \
      mfc_write_tag_mask( mask );                                       \
      mfc_read_tag_status_all();                                        \
    }                                                                   \
    np_block[buffer] = 0;                                               \
    nm_block[buffer] = 0;                                               \
  } while(0)

  p_block[0] = p_stream;              np_block[0] = 0;
  p_block[1] = p_stream + NP_BLOCK;   np_block[1] = 0;
  p_block[2] = p_stream + NP_BLOCK*2; np_block[2] = 0;

  m_block[0] = m_stream;              nm_block[0] = 0;
  m_block[1] = m_stream + NP_BLOCK;   nm_block[1] = 0;
  m_block[2] = m_stream + NP_BLOCK*2; nm_block[2] = 0;

  voxel_cache_init();

  BEGIN_GET_PBLOCK(0);
  // BEGIN_PUT_PBLOCK( 2 );
  for( buffer=0; np_block[buffer]>0; buffer=next[buffer] ) {
    BEGIN_GET_PBLOCK( next[buffer] );
    END_GET_PBLOCK( buffer );
    PROCESS_PBLOCK( buffer );
    BEGIN_PUT_PBLOCK( buffer );
    END_PUT_PBLOCK( prev[buffer] );
  }
  END_PUT_PBLOCK( prev[buffer] );

  // Writeback all the cached accumulators

  voxel_cache_writeback();

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

