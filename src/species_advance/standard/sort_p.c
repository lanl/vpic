#define IN_spa
#define HAS_SPU_PIPELINE
#include "spa_private.h"

#if defined(__SSE__)
#include <xmmintrin.h>
#endif

// FIXME: HOOK UP IN-PLACE / OUT-PLACE OPTIONS AGAIN
// FIXME: Add N_VOXEL convenience field to grid

// FIXME: ADD RESTRICT TO UTIL_BASE.H
#ifndef RESTRICT
#define RESTRICT __restrict
#endif

void
subsort_pipeline( sort_p_pipeline_args_t * args,
                  int pipeline_rank,
                  int n_pipeline ) {
  const particle_t * RESTRICT ALIGNED(128) p_src;
  /**/  particle_t * RESTRICT ALIGNED(128) p_dst;
  int i0, i1, v0, v1, i, j, v, sum, count;

  int * RESTRICT ALIGNED(128) partition = args->partition;
  int * RESTRICT ALIGNED(128) next      = args->next;

  if( pipeline_rank==n_pipeline ) return; // No straggler cleanup needed
 
  // This pipeline sorts particles in this voxel range in [i0,i1) in
  // the aux array.  These particles are in voxels [v0,v1).
  p_src = args->aux_p;
  p_dst = args->p;
  i0    = args->coarse_partition[pipeline_rank  ];
  i1    = args->coarse_partition[pipeline_rank+1];
  v0    = P2V( pipeline_rank,   n_pipeline, args->vl, args->vh );
  v1    = P2V( pipeline_rank+1, n_pipeline, args->vl, args->vh );
  if( pipeline_rank==0            ) v0 = 0;
  if( pipeline_rank==n_pipeline-1 ) v1 = args->n_voxel;

  // Clear fine grained count
  CLEAR( &next[v0], v1-v0 );

  // Fine grained count
  for( i=i0; i<i1; i++ ) next[ p_src[i].i ]++;

  // Compute the partitioning
  sum = i0;
  for( v=v0; v<v1; v++ ) {
    count = next[v];
    next[v] = sum;
    partition[v] = sum;
    sum += count;
  }
  partition[v1] = sum; // All threads who write this agree

  // Local fine grained sort
  for( i=i0; i<i1; i++ ) {
    v = p_src[i].i;
    j = next[v]++;
#   if defined(__SSE__)
    _mm_store_ps( &p_dst[j].dx, _mm_load_ps( &p_src[i].dx ) );
    _mm_store_ps( &p_dst[j].ux, _mm_load_ps( &p_src[i].ux ) );
#   else
    p_dst[j] = p_src[i];
#   endif
  }
}

#define VOXEL(x,y,z) INDEX_FORTRAN_3(x,y,z,0,g->nx+1,0,g->ny+1,0,g->nz+1)

void
sort_p( species_t * sp,
        const grid_t * g ) {

  static char * ALIGNED(128) scratch = NULL;
  static size_t max_scratch = 0;
  size_t sz_scratch;

  int n_pipeline = N_PIPELINE;
  int n_voxel    = (g->nx+2)*(g->ny+2)*(g->nz+2);

  DECLARE_ALIGNED_ARRAY( sort_p_pipeline_args_t, 128, args, 1 );
  DECLARE_ALIGNED_ARRAY( int, 128, coarse_partition,
                         max_subsort_pipeline*(MAX_PIPELINE+1) );

  // FIXME: TEMPORARY HACK UNTIL SPECIES_ADVANCE API INSTALLED
  if( sp->partition==NULL ) MALLOC_ALIGNED( sp->partition, n_voxel+1, 128 );

  // Insure enough scratch space is allocated for the sorting
  sz_scratch = ( sizeof(args->aux_p[0])*sp->np +            /* For aux_p */
                 128 + sizeof(args->next[0])*(n_voxel+1) ); /* For next */
  if( sz_scratch > max_scratch ) {
    FREE_ALIGNED( scratch );
    MALLOC_ALIGNED( scratch, sz_scratch, 128 );
    max_scratch = sz_scratch;
  }

  // Setup pipeline arguments
  // FIXME: NOTE THAT THE SPU PIPELINE DOESN'T ACTUALLY USE NEXT.
  args->p                  = sp->p;
  args->aux_p              = (particle_t *)scratch;
  args->coarse_partition   = coarse_partition;
  args->partition          = sp->partition;
  args->next               = (int *)( args->aux_p + sp->np );//ALIGN_PTR( int, args->aux_p + sp->np, 128 );
  args->n                  = sp->np;
  args->n_subsort_pipeline = n_pipeline;
  args->vl                 = VOXEL(1,1,1);
  args->vh                 = VOXEL(g->nx,g->ny,g->nz);
  args->n_voxel            = n_voxel;

  if( n_pipeline>1 ) {

    // Do coarse sort
    coarse_sort_p( args );

    // Convert the coarse_partitioning used durign the coarse sort
    // into the partitioning of the particle list by coarse bucket
    coarse_partition[ n_pipeline ] = sp->np;

    // Do fine grained subsorts
    EXEC_PIPELINES( subsort, args, 0 );
    WAIT_PIPELINES();

  } else {

    // Just do the subsort when single threaded.  We need to hack the
    // aux arrays and what not to make it look like coarse sorting was
    // done to the subsort pipeline.
    coarse_partition[0]          = 0;
    coarse_partition[n_pipeline] = sp->np;
    args->p     = (particle_t *)scratch;
    args->aux_p = sp->p;
    subsort_pipeline( args, 0, 1 );

    // Results ended up in the wrong place as a result of the ugly
    // hack above.  Copy it to the right place and undo the above
    // hack.  FIXME: IF WILLING TO MOVE SP->P AROUND AND DO MORE
    // MALLOCS PER STEP (I.E. HEAP FRAGMENTATION), COULD AVOID THIS
    // COPY.
    COPY( sp->p, (particle_t *)scratch, sp->np );

  }
}

#if 0 // In-place, single threaded legacy version

void
sort_p( species_t * sp,
        const grid_t * g ) {
  particle_t * ALIGNED(128) p = sp->p;
  //const int32_t * __restrict ALIGNED(128) sfc = g->sfc;
  const int np                = sp->np; 
  const int nc                = (g->nx+2)*(g->ny+2)*(g->nz+2);
  const int nc1               = nc+1;
  int * __restrict ALIGNED(128) partition = sp->partition;

  static int * __restrict ALIGNED(128) next = NULL;
  static int max_nc1 = 0;

  int i, j;

  // FIXME: TEMPORARY HACK UNTIL SPECIES_ADVANCE API INSTALLED
  if( sp->partition==NULL ) MALLOC_ALIGNED( sp->partition, nc1, 128 );
  partition = sp->partition;

  if( np==0 ) return; // Do not need to sort

  // Allocate the sorting intermediate
  // Making this into a static is done to avoid heap shredding
 
  if( max_nc1<nc1 ) {
    int * tmp = next; // Hack around __restrict__ issues
    FREE_ALIGNED(   tmp );
    MALLOC_ALIGNED( tmp, nc1, 128 );
    next    = tmp;
    max_nc1 = nc1;
  }

  // Count particles in each cell
  CLEAR( next, nc1 );
  //for( i=0; i<np; i++ ) next[ sfc[ p[i].i ] ]++;
  for( i=0; i<np; i++ ) next[ p[i].i ]++;

  // Convert the count to a partitioning (and save a copy in next)
  j=0;
  for( i=0; i<nc1; i++ ) {
    partition[i] = j;
    j += next[i];
    next[i] = partition[i];
  }

  if( sp->sort_out_of_place ) {

    // Throw down the particle array in order

    particle_t * ALIGNED(128) new_p;
    const particle_t * __restrict ALIGNED(32) in_p;
    /**/  particle_t * __restrict ALIGNED(32) out_p;

    MALLOC_ALIGNED( new_p, sp->max_np, 128 );

    in_p  = sp->p;
    out_p = new_p;
    //for( i=0; i<np; i++ ) out_p[ next[ sfc[ in_p[i].i ] ]++ ] = in_p[i];
    for( i=0; i<np; i++ ) out_p[ next[ in_p[i].i ]++ ] = in_p[i];

    FREE_ALIGNED( sp->p );
    sp->p = new_p;

  } else {

    // Run sort cycles until the list is sorted

    particle_t save_p, * ALIGNED(32) src, * ALIGNED(32) dest;

    i=0;
    while( i<nc ) {
      if( next[i]>=partition[i+1] ) i++;
      else {
        src = &p[ next[i] ];
        for(;;) {
          //dest = &p[ next[ sfc[ src->i ] ]++ ];
          dest = &p[ next[ src->i ]++ ];
          if( src==dest ) break;
          save_p = *dest;
          *dest  = *src;
          *src   = save_p;
        }
      }

    }
  }
}

#endif
