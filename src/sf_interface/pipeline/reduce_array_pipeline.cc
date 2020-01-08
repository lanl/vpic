#define IN_sf_interface

// NOTE: Experimentation shows no benefit to explicitly vectorizing array
// reductions, but maintaining a v4 version just in case.
#define HAS_V4_PIPELINE
#define HAS_V8_PIPELINE
#define HAS_V16_PIPELINE

#include "sf_interface_pipeline.h"

#include "../sf_interface_private.h"

#include "../../util/pipelines/pipelines_exec.h"

// FIXME: N_ARRAY>1 ALWAYS BUT THIS ISN'T STRICTLY NECESSARY BECAUSE
// HOST IS THREAD FOR THE SERIAL AND THREADED DISPATCHERS.  SHOULD
// PROBABLY CHANGE N_ARRAY TO
// max({serial,thread}.n_pipeline,spu.n_pipeline+1)

void
reduce_array_pipeline_scalar( reduce_pipeline_args_t * args,
	                            int pipeline_rank,
	                            int n_pipeline ) {

  int i, j, i1, r0;
  int r, nr = args->n_array-1, sr = args->s_array;

  DISTRIBUTE( args->n, 1, pipeline_rank, n_pipeline, i, i1 ); i1 += i;

  // a is broken into restricted rw and ro parts to allow the compiler
  // to do more aggresive optimizations

  /**/  float * RESTRICT ALIGNED(16) a = args->a;
  const float * RESTRICT ALIGNED(16) b = a + sr;

  float v0, v1, v2, v3, v4, v5, v6, v7, v8, v9;

# define LOOP(OP) for(j=i ; j<i1; j+=1) { OP(j); }

# define A(k)   v0 = a[k];
# define B(k,r) v##r = b[k+(r+r0)*sr];
# define C(k,v) a[k] = v
# define O1(k)A(k  )B(k,1)                                                 \
    C(k,   v0+v1)
# define O2(k)A(k  )B(k,1)B(k,2)                                           \
    C(k,  (v0+v1)+ v2)
# define O3(k)A(k  )B(k,1)B(k,2)B(k,3)                                     \
    C(k,  (v0+v1)+(v2+v3))
# define O4(k)A(k  )B(k,1)B(k,2)B(k,3)B(k,4)                               \
    C(k, ((v0+v1)+(v2+v3))+  v4)
# define O5(k)A(k  )B(k,1)B(k,2)B(k,3)B(k,4)B(k,5)                         \
    C(k, ((v0+v1)+(v2+v3))+ (v4+v5))
# define O6(k)A(k  )B(k,1)B(k,2)B(k,3)B(k,4)B(k,5)B(k,6)                   \
    C(k, ((v0+v1)+(v2+v3))+((v4+v5)+ v6))
# define O7(k)A(k  )B(k,1)B(k,2)B(k,3)B(k,4)B(k,5)B(k,6)B(k,7)             \
    C(k, ((v0+v1)+(v2+v3))+((v4+v5)+(v6+v7)))
# define O8(k)A(k  )B(k,1)B(k,2)B(k,3)B(k,4)B(k,5)B(k,6)B(k,7)B(k,8)       \
    C(k,(((v0+v1)+(v2+v3))+((v4+v5)+(v6+v7)))+   v8)
# define O9(k)A(k  )B(k,1)B(k,2)B(k,3)B(k,4)B(k,5)B(k,6)B(k,7)B(k,8)B(k,9) \
    C(k,(((v0+v1)+(v2+v3))+((v4+v5)+(v6+v7)))+  (v8+v9))

  r0 = -1;
  while( nr ) {
    switch( nr ) {
    case  0:                               break;
    case  1: LOOP(O1); nr -= 1 ; r0 += 1 ; break;
    case  2: LOOP(O2); nr -= 2 ; r0 += 2 ; break;
    case  3: LOOP(O3); nr -= 3 ; r0 += 3 ; break;
    case  4: LOOP(O4); nr -= 4 ; r0 += 4 ; break;
    case  5: LOOP(O5); nr -= 5 ; r0 += 5 ; break;
    case  6: LOOP(O6); nr -= 6 ; r0 += 6 ; break;
    case  7: LOOP(O7); nr -= 7 ; r0 += 7 ; break;
    case  8: LOOP(O8); nr -= 8 ; r0 += 8 ; break;
    default: LOOP(O9); nr -= 9 ; r0 += 9 ; break;
    }
  }

# undef O9
# undef O8
# undef O7
# undef O6
# undef O5
# undef O4
# undef O3
# undef O2
# undef O1
# undef C
# undef B
# undef A
# undef LOOP

}

#define VOX(x,y,z) VOXEL( x, y, z, aa->g->nx, aa->g->ny, aa->g->nz )

void
reduce_accumulator_array_pipeline( accumulator_array_t * RESTRICT aa )
{
  DECLARE_ALIGNED_ARRAY( reduce_pipeline_args_t, 128, args, 1 );

  int i0, na, nfloats;

  if ( !aa )
  {
    ERROR( ( "Bad args" ) );
  }

  i0 = ( VOX(1,1,1) / 2 ) * 2; // Round i0 down to even for 128B align on Cell
  na = ( ( ( VOX( aa->g->nx, aa->g->ny, aa->g->nz ) - i0 + 1 ) + 1 ) / 2 ) * 2;

  nfloats       = sizeof(accumulator_t)/sizeof(float);
  args->a       = (float *)(aa->a + i0);
  args->n       = na*nfloats;
  args->n_array = aa->n_pipeline + 1;
  args->s_array = aa->stride*nfloats;
  args->n_block = accumulators_n_block;

  EXEC_PIPELINES( reduce_array, args, 0 );

  WAIT_PIPELINES();
}

#undef VOX

void
reduce_hydro_array_pipeline( hydro_array_t * RESTRICT ha )
{
  DECLARE_ALIGNED_ARRAY( reduce_pipeline_args_t, 128, args, 1 );

  int nfloats;

  if ( !ha )
  {
    ERROR( ( "Bad args" ) );
  }

  nfloats       = sizeof(hydro_t)/sizeof(float);
  args->a       = (float *)(ha->h);
  args->n       = ha->g->nv*nfloats;
  args->n_array = ha->n_pipeline + 1;
  args->s_array = ha->stride*nfloats;
  args->n_block = hydro_n_block;

  EXEC_PIPELINES( reduce_array, args, 0 );

  WAIT_PIPELINES();
}
