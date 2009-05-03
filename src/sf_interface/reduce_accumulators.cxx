#define IN_sf_interface
#include "sf_interface_private.h"

// FIXME: THIS NEEDS TO REDUCE THE ACTUAL NUMBER OF ACCUMULATORS PRESENT
// AS THE NUMBER OF PIPELINES EXECUTED HERE MAY NOT MATCH THE NUMBER OF
// ACCUMULATORS USED DURING OTHER OPERATIONS (E.G. 8 SPU PIPELINES AND
// 2 PPU PIPELINES!)

// FIXME: N_ARRAY>1 ALWAYS BUT THIS ISN'T STRICTLY NECESSARY BECAUSE
// HOST IS THREAD FOR THE SERIAL AND THREADED DISPATCHERS.  SHOULD
// PROBABLY CHANGE N_ARRAY T
// max({serial,thread}.n_pipeline,spu.n_pipeline+1)

void
reduce_accumulators_pipeline( accumulators_pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline ) {
  int si = sizeof(accumulator_t) / sizeof(float);
  int sj = si*args->stride;   
  int nj = args->n_array - 1;                     
  int i, j, k, l, i1;

  /* a is broken into restricted rw and ro parts to allow the compiler
     to do more aggresive optimizations. */
  /**/  float * __restrict a = args->a->jx;
  const float * __restrict b = a + sj;

  DISTRIBUTE( args->stride, accumulators_n_block,
              pipeline_rank, n_pipeline, i, i1 ); i1 += i;

# if defined(V4_ACCELERATION)

  using namespace v4;

  v4float v0, v1, v2, v3, v4, v5, v6, v7;

# define LOOP(_)                                \
  for( ; i<i1; i++ ) {                          \
    k = i*si;                                   \
    _(k+ 0); _(k+ 4); _(k+ 8);                  \
  }
# define A(k)    load_4x1( &a[k], v0 );
# define B(k,j)  load_4x1( &b[k+(j-1)*sj], v##j );
# define C(k,op) store_4x1( op, &a[k] )
# define O1(_)A(_  )B(_,1)            C(_,  v0+v1                            )
# define O2(_)A(_  )B(_,1)B(_,2)      C(_, (v0+v1)+ v2                       )
# define O3(_)A(_  )B(_,1)B(_,2)B(_,3)C(_, (v0+v1)+(v2+v3)                   )
# define O4(_)A(_  )B(_,1)B(_,2)B(_,3)                                  \
    /**/      B(_,4)                  C(_,((v0+v1)+(v2+v3))+  v4             )
# define O5(_)A(_  )B(_,1)B(_,2)B(_,3)                                  \
    /**/      B(_,4)B(_,5)            C(_,((v0+v1)+(v2+v3))+ (v4+v5)         )
# define O6(_)A(_  )B(_,1)B(_,2)B(_,3)                                  \
    /**/      B(_,4)B(_,5)B(_,6)      C(_,((v0+v1)+(v2+v3))+((v4+v5)+ v6    ))
# define O7(_)A(_  )B(_,1)B(_,2)B(_,3)                                  \
    /**/      B(_,4)B(_,5)B(_,6)B(_,7)C(_,((v0+v1)+(v2+v3))+((v4+v5)+(v6+v7)))

# else

  float f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11;

# define LOOP(_)                                \
  for( ; i<i1; i++ ) {                          \
    k = i*si;                                   \
    _(k+ 0);_(k+ 1);_(k+ 2);_(k+ 3);            \
    _(k+ 4);_(k+ 5);_(k+ 6);_(k+ 7);            \
    _(k+ 8);_(k+ 9);_(k+10);_(k+11);            \
  }
# define O1(_) a[_]=  a[_     ] + b[_     ]
# define O2(_) a[_]= (a[_     ] + b[_     ]) +  b[_+  sj]
# define O3(_) a[_]= (a[_     ] + b[_     ]) + (b[_+  sj] + b[_+2*sj])
# define O4(_) a[_]=((a[_     ] + b[_     ]) + (b[_+  sj] + b[_+2*sj])) + \
    /**/              b[_+3*sj]
# define O5(_) a[_]=((a[_     ] + b[_     ]) + (b[_+  sj] + b[_+2*sj])) + \
    /**/             (b[_+3*sj] + b[_+4*sj])
# define O6(_) a[_]=((a[_     ] + b[_     ]) + (b[_+  sj] + b[_+2*sj])) + \
    /**/            ((b[_+3*sj] + b[_+4*sj]) +  b[_+5*sj]          )
# define O7(_) a[_]=((a[_     ] + b[_     ]) + (b[_+  sj] + b[_+2*sj])) + \
    /**/            ((b[_+3*sj] + b[_+4*sj]) + (b[_+5*sj] + b[_+6*sj]))

# endif
  
  switch( nj ) {        
  case 0:           break;
  case 1: LOOP(O1); break;
  case 2: LOOP(O2); break;
  case 3: LOOP(O3); break;
  case 4: LOOP(O4); break;
  case 5: LOOP(O5); break;
  case 6: LOOP(O6); break;
  case 7: LOOP(O7); break;
  default:
#   if defined(V4_ACCELERATION)
    for( ; i<i1; i++ ) {
      k = i*si;
      load_4x1(&a[k+0],v0);  load_4x1(&a[k+4],v1);  load_4x1(&a[k+8],v2);
      for( j=0; j<nj; j++ ) {
        l = k + j*sj;
        load_4x1(&b[l+0],v3);  load_4x1(&b[l+4],v4);  load_4x1(&b[l+8],v5);
        v0 += v3;              v1 += v4;              v2 += v5;
      }
      store_4x1(v0,&a[k+0]); store_4x1(v1,&a[k+4]); store_4x1(v2,&a[k+8]);
    }
#   else
    for( ; i<i1; i++ ) {
      k = i*si;
      f0  = a[k+ 0]; f1  = a[k+ 1]; f2  = a[k+ 2]; f3  = a[k+ 3];
      f4  = a[k+ 4]; f5  = a[k+ 5]; f6  = a[k+ 6]; f7  = a[k+ 7];
      f8  = a[k+ 8]; f9  = a[k+ 9]; f10 = a[k+10]; f11 = a[k+11];
      for( j=0; j<nj; j++ ) {
        l = k + j*sj;
        f0  += b[l+ 0]; f1  += b[l+ 1]; f2  += b[l+ 2]; f3  += b[l+ 3];
        f4  += b[l+ 4]; f5  += b[l+ 5]; f6  += b[l+ 6]; f7  += b[l+ 7];
        f8  += b[l+ 8]; f9  += b[l+ 9]; f10 += b[l+10]; f11 += b[l+11];
      }
      a[k+ 0] =  f0; a[k+ 1] =  f1; a[k+ 2] =  f2; a[k+ 3] =  f3;
      a[k+ 4] =  f4; a[k+ 5] =  f5; a[k+ 6] =  f6; a[k+ 7] =  f7;
      a[k+ 8] =  f8; a[k+ 9] =  f9; a[k+10] = f10; a[k+11] = f11;
    }
#   endif
    break;
  }

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

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

// SPU pipeline is defined in a different compilation unit

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "The regular pipeline is already V4 accelerated!"

#endif

void
reduce_accumulators( accumulator_t * ALIGNED(128) a,
                     const grid_t  *              g ) {
  accumulators_pipeline_args_t args[1];
  int n_array, stride;

  if( a==NULL ) ERROR(("Invalid accumulator"));
  if( g==NULL ) ERROR(("Invalid grid"));

  /**/                            n_array = serial.n_pipeline;
  if( n_array<thread.n_pipeline ) n_array = thread.n_pipeline;
# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  if( n_array<spu.n_pipeline    ) n_array = spu.n_pipeline;
# endif
  n_array++; /* n_array = 1 + max( {serial,thread,spu}.n_pipeline ) */

  stride = POW2_CEIL((g->nx+2)*(g->ny+2)*(g->nz+2),2);

  if( n_array==1 ) return; /* Nothing to do */
  /* FIXME: ADD SHORTCUT FOR WHEN STRIDE IS TOO SMALL TO BE WORTH
     THREADING PARALLELIZING */
  args->a = a, args->n_array = n_array, args->stride = stride;
  EXEC_PIPELINES( reduce_accumulators, args, 0 );
  WAIT_PIPELINES();
}
