/* FIXME: independent accumulator and interpolator caches for charge
   conserving accumulate. */

typedef struct voxel_cache {

  /* Four way vectorizable software cache */

  uint32_t clock; 

  int32_t        * ALIGN voxel;        /* n_lines_per_way x n_way */
  uint32_t       * ALIGN last_access;  /* n_lines_per_way x n_way */
  interpolator_t * ALIGN interpolator; /* n_lines_per_way x n_way */
  accumulator_t  * ALIGN accumulator;  /* n_lines_per_way x n_way */

  const interpolator_t * ALIGN FAR mem_interpolator;
  accumulator_t        * ALIGN FAR mem_accumulator;
  int n_line_per_way; /* Must be a power of two */

} voxel_cache_t;

voxel_cache_t *
new_voxel_cache( const interpolator_t * ALIGN FAR mem_interpolator,
                 accumulator_t        * ALIGN FAR mem_accumulator,
                 int n_line_per_way ) {
  voxel_cache_t * cache;
  int n_line;

  if( !IS_POWER_OF_TWO( n_line_per_way ) ) return NULL;

  cache = malloc( sizeof( voxel_cache_t ) );
  if( cache==NULL ) return NULL;

  n_line = n_line_per_way * SPU_VOXEL_CACHE_N_WAY;

  cache->voxel        = malloc_aligned( n_line*sizeof(int32_t) );
  cache->last_access  = malloc_aligned( n_line*sizeof(uint32_t) );
  cache->interpolator = malloc_aligned( n_line*sizeof(interpolator_t) );
  cache->accumulator  = malloc_aligned( n_line*sizeof(accumulator_t) );

  if( cache->voxel==NULL        ||
      cache->last_access==NULL  ||
      cache->interpolator==NULL ||
      cache->accumulator==NULL  ) {
    return NULL;
  }

  memset( cache->last_access, 0, n_line*sizeof(uint32_t) );
  cache->clock = 1;

  cache->mem_interpolator = mem_interpolator;
  cache->mem_accumulator  = mem_accumulator;
  cache->n_line_per_way   = n_line_per_way;

  return cache;
}

/* Given a local voxel index, return the voxel cache line that contains
   the appropriate data. */

/* FIXME: Quoad fetch_voxel_cache could be done too! */

inline int32_t 
fetch_voxel_cache( voxel_cache_t * cache,
                   int32_t voxel ) {
  int32_t line;
  int32_t cached_voxel[4];

  /* Determine which cache lines could contain the voxel */

  line = ( voxel & ( cache->n_line_per_way-1 ) ) << 2;

  /* Quad load */

  cached_voxel[0] = cache->voxel[line+0];
  cached_voxel[1] = cache->voxel[line+1];
  cached_voxel[2] = cache->voxel[line+2];
  cached_voxel[3] = cache->voxel[line+3];

  /* Quad return index of matching element (any matching element will do) */

  way = 4;
  if( cached_voxel[0] == voxel ) way = 0;
  if( cached_voxel[1] == voxel ) way = 1;
  if( cached_voxel[2] == voxel ) way = 2;
  if( cached_voxel[3] == voxel ) way = 3;

  if( way==4 ) { /* Sigh ... cache miss */

    /* Find the least recently used cache line. */
    /* Note: Rollover handling could be done here but probably not necesary */
    /* Note: Half swizzle could be used for fetches of age and voxel */

    /* Quad load */

    last_access[0] = cache->last_access[line+0];
    last_access[1] = cache->last_access[line+1];
    last_access[2] = cache->last_access[line+2];
    last_access[3] = cache->last_access[line+3];
  
    /* Quad return smallest element */
    /* Quad return index of smallest element (any smallest element will do) */

    /**/                               way = 0, last_access = last_access[0];
    if( last_access[1] < last_access ) way = 1, last_access = last_access[1];
    if( last_access[2] < last_access ) way = 2, last_access = last_access[2];
    if( last_access[3] < last_access ) way = 3, last_access = last_access[3];

    /* Replace the least recently used cache line */

    /* Flush the LRU cache line */

    if( last_access!=0 ) { 
      /* Note: interpolator is read-only so no need to flush */
      /* Store cache->accumulator[ line + way ] into
               cache->mem_accumulator[ cached_voxel[ way ] ] */
    }

    /* Load the desired voxel into LRU cache line */

    /* Load  cache->interpolator[ line + way ] with
             external_memory->interpolator[ cached_voxel[ way ] ] 
       Load  cache->accumulator[ line + way ] with
             cache->mem_accumulator[ cached_voxel[ way ] ] */
    cache->voxel[ line + way ] = voxel;
  }

  /* line:way contains the desired voxel data.  Save the time at which
     this cache line was last updated.  last_access==0 indicates the
     cache line has never been accessed so clock may not take that
     value.  Overlapping is probably unnecessary in practice. */

  cache->last_access[ line + way ] = cache->clock;
  cache->clock++; if( cache->clock==0 ) cache->clock = 1;

  return line + way;
}

/* Flush all entires in a cache to memory */

void
flush_voxel_cache( voxel_cache_t * cache ) {
  for( line=0; line<cache->n_line_per_way*4; line++ ) {
    if( cache->last_access[line]==0 ) continue;
    /* Note: interpolator is read-only
       Store cache->accumulator[ line + way ] into
             cache->mem_accumulator[ cached_voxel[ way ] ] */
    cache->last_access[line] = 0;
  }
  cache->clock = 1;
}

void
delete_voxel_cache( voxel_cache_t * cache ) {
  free_aligned( cache->voxel );
  free_aligned( cache->last_access );
  free_aligned( cache->interpolator );
  free_aligned( cache->accumulator );
  free( cache );
}

cache = new_voxel_cache( mem_accumulator, mem_interpolator, 512 );

operton_advance_p() {
  call cell_advance_p
}

cell_advance_p() {
  for each spu,
    call spu_advance_p function on a range of particle quads,
  end for
  advance_p for stragglers
}

spu_advance_p() {

  ... triple buffered particle data streaming from memory 
  ... one dma load and one dma store occuring during particle processing
  ... advance_p uses software managed 4-way associative simple LRU caches to
  ... minimize number of interpolator, accumulator, cell connectivity and
  ... other transformations

  new voxel cache[0]

  new pbuffer[0]
  new pbuffer[1]
  new pbuffer[2]

  begin dma load pblock[0] into pbuffer[0]

  for( i=0, j=0; i<n_blocks; i++, j=(j+1)%3 ) {

    end   dma load pblock[i  ] into pbuffer[ j     ]
    begin dma load pblock[i+1] into pbuffer[(j+1)%3]

    spu_quad_advance_p on pbuffer[j]                 

    if i>0, end   dma store pbuffer[(j-1)%3] to pblock[i-1]
    /**/    begin dma store pbuffer[ j     ] to pblock[i  ]

  }

  end   dma store buffer[(j-1)%3] to block[i-1]
}

... guard list particles should be written direct to memory
... fixme: cell mesh connectivity caching!

void
spu_quad_advance_p( ) {
  for( ; ; ) {
    ... load particles into memory ...
    /* FIXME: quad_fetch_voxel_cache */
    l0 = fetch_voxel_cache( cache, p0->voxel );
    l1 = fetch_voxel_cache( cache, p1->voxel );
    l2 = fetch_voxel_cache( cache, p2->voxel );
    l3 = fetch_voxel_cache( cache, p3->voxel );
    ...  rest of regular quad advance ...
    ... fancy interpolate functions need to be written ...
  }
}

