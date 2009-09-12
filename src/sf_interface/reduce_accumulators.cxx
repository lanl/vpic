#define IN_sf_interface
#define HAS_SPU_PIPELINE
#include "sf_interface_private.h"

// FIXME: N_ARRAY>1 ALWAYS BUT THIS ISN'T STRICTLY NECESSARY BECAUSE
// HOST IS THREAD FOR THE SERIAL AND THREADED DISPATCHERS.  SHOULD
// PROBABLY CHANGE N_ARRAY TO
// max({serial,thread}.n_pipeline,spu.n_pipeline+1)

void
reduce_accumulators_pipeline( accumulators_pipeline_args_t * args,
                              int pipeline_rank,
                              int n_pipeline ) {
  int i, i1, si = sizeof(accumulator_t) / sizeof(float);
  int r, nr = args->n_array-1, sr = si*args->s_array;
  int j, k;

  DISTRIBUTE( args->n, accumulators_n_block,
              pipeline_rank, n_pipeline, i, i1 ); i1 += i;

  // a is broken into restricted rw and ro parts to allow the compiler
  // to do more aggresive optimizations

  /**/  float * RESTRICT ALIGNED(16) a = args->a->jx;
  const float * RESTRICT ALIGNED(16) b = a + sr;

# if defined(V4_ACCELERATION)

  using namespace v4;

  v4float v0, v1, v2, v3, v4, v5, v6, v7, v8, v9;

# define LOOP(OP)                               \
  for( ; i<i1; i++ ) {                          \
    k = i*si;                                   \
    OP(k   ); OP(k+ 4); OP(k+ 8);               \
  }
# define A(k)   load_4x1(  &a[k],          v0   );
# define B(k,r) load_4x1(  &b[k+(r-1)*sr], v##r );
# define C(k,v) store_4x1( v, &a[k] )
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

# else

  float f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11;

# define LOOP(OP)                               \
  for( ; i<i1; i++ ) {                          \
    k = i*si;                                   \
    OP(k   );OP(k+ 1);OP(k+ 2);OP(k+ 3);        \
    OP(k+ 4);OP(k+ 5);OP(k+ 6);OP(k+ 7);        \
    OP(k+ 8);OP(k+ 9);OP(k+10);OP(k+11);        \
  }
# define O1(k) a[k] =    a[k     ] + b[k     ]
# define O2(k) a[k] =   (a[k     ] + b[k     ]) +  b[k+  sr]
# define O3(k) a[k] =   (a[k     ] + b[k     ]) + (b[k+  sr] + b[k+2*sr])
# define O4(k) a[k] =  ((a[k     ] + b[k     ]) + (b[k+  sr] + b[k+2*sr]))  + \
    /**/                 b[k+3*sr]
# define O5(k) a[k] =  ((a[k     ] + b[k     ]) + (b[k+  sr] + b[k+2*sr]))  + \
    /**/                (b[k+3*sr] + b[k+4*sr])
# define O6(k) a[k] =  ((a[k     ] + b[k     ]) + (b[k+  sr] + b[k+2*sr]))  + \
    /**/               ((b[k+3*sr] + b[k+4*sr]) +  b[k+5*sr]          )
# define O7(k) a[k] =  ((a[k     ] + b[k     ]) + (b[k+  sr] + b[k+2*sr]))  + \
    /**/               ((b[k+3*sr] + b[k+4*sr]) + (b[k+5*sr] + b[k+6*sr]))
# define O8(k) a[k] = (((a[k     ] + b[k     ]) + (b[k+  sr] + b[k+2*sr]))  + \
    /**/               ((b[k+3*sr] + b[k+4*sr]) + (b[k+5*sr] + b[k+6*sr]))) + \
    /**/                 b[k+7*sr]
# define O9(k) a[k] = (((a[k     ] + b[k     ]) + (b[k+  sr] + b[k+2*sr]))  + \
    /**/               ((b[k+3*sr] + b[k+4*sr]) + (b[k+5*sr] + b[k+6*sr]))) + \
    /**/                (b[k+7*sr] + b[k+8*sr])

# endif
  
  switch( nr ) {        
  case 0:           break;
  case 1: LOOP(O1); break;
  case 2: LOOP(O2); break;
  case 3: LOOP(O3); break;
  case 4: LOOP(O4); break;
  case 5: LOOP(O5); break;
  case 6: LOOP(O6); break;
  case 7: LOOP(O7); break;
  case 8: LOOP(O8); break;
  case 9: LOOP(O9); break;
  default:
#   if defined(V4_ACCELERATION)
    for( ; i<i1; i++ ) {
      j = i*si;
      load_4x1(&a[j+0],v0);  load_4x1(&a[j+4],v1);  load_4x1(&a[j+8],v2);
      for( r=0; r<nr; r++ ) {
        k = j + r*sr;
        load_4x1(&b[k+0],v3);  load_4x1(&b[k+4],v4);  load_4x1(&b[k+8],v5);
        v0 += v3;              v1 += v4;              v2 += v5;
      }
      store_4x1(v0,&a[j+0]); store_4x1(v1,&a[j+4]); store_4x1(v2,&a[j+8]);
    }
#   else
    for( ; i<i1; i++ ) {
      j = i*si;
      f0  = a[j+ 0]; f1  = a[j+ 1]; f2  = a[j+ 2]; f3  = a[j+ 3];
      f4  = a[j+ 4]; f5  = a[j+ 5]; f6  = a[j+ 6]; f7  = a[j+ 7];
      f8  = a[j+ 8]; f9  = a[j+ 9]; f10 = a[j+10]; f11 = a[j+11];
      for( r=0; r<nr; r++ ) {
        k = j + r*sr;
        f0  += b[k+ 0]; f1  += b[k+ 1]; f2  += b[k+ 2]; f3  += b[k+ 3];
        f4  += b[k+ 4]; f5  += b[k+ 5]; f6  += b[k+ 6]; f7  += b[k+ 7];
        f8  += b[k+ 8]; f9  += b[k+ 9]; f10 += b[k+10]; f11 += b[k+11];
      }
      a[j+ 0] =  f0; a[j+ 1] =  f1; a[j+ 2] =  f2; a[j+ 3] =  f3;
      a[j+ 4] =  f4; a[j+ 5] =  f5; a[j+ 6] =  f6; a[j+ 7] =  f7;
      a[j+ 8] =  f8; a[j+ 9] =  f9; a[j+10] = f10; a[j+11] = f11;
    }
#   endif
    break;
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

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

// SPU pipeline is defined in a different compilation unit

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "The regular pipeline is already V4 accelerated!"

#endif

#define VOX(x,y,z) VOXEL(x,y,z, aa->g->nx,aa->g->ny,aa->g->nz)

void
reduce_accumulator_array( accumulator_array_t * RESTRICT aa ) {
  DECLARE_ALIGNED_ARRAY( accumulators_pipeline_args_t, 128, args, 1 );
  int i0;

  if( !aa ) ERROR(( "Bad args" ));

  i0 = (VOX(1,1,1)/2)*2; // Round i0 down to even for 128B align on Cell

  args->a       = aa->a + i0;
  args->n       = (((VOX(aa->g->nx,aa->g->ny,aa->g->nz) - i0 + 1 )+1)/2)*2;
  args->n_array = aa->n_pipeline + 1;
  args->s_array = aa->stride;

  EXEC_PIPELINES( reduce_accumulators, args, 0 );
  WAIT_PIPELINES();
}

