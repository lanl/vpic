#define IN_sf_interface

#include "sf_interface_pipeline.h"

#include "../sf_interface_private.h"

#if defined(V8_ACCELERATION)

using namespace v8;

void
reduce_array_pipeline_v8( reduce_pipeline_args_t * args,
                          int pipeline_rank,
                          int n_pipeline )
{
  int i, i1, j, r0;
  int nr = args->n_array - 1;
  int sr = args->s_array;

  DISTRIBUTE( args->n, 8, pipeline_rank, n_pipeline, i, i1 );

  i1 += i;

  // a is broken into restricted rw and ro parts to allow the compiler
  // to do more aggresive optimizations.

  /**/  float * RESTRICT ALIGNED(16) a = args->a;
  const float * RESTRICT ALIGNED(16) b = a + sr;

  v8float v0, v1, v2, v3, v4, v5, v6, v7, v8, v9;

  #define LOOP( OP ) for( j = i ; j < i1; j += 8 ) { OP( j ); }

  #define A( k )     load_8x1( &a[ k ], v0 );

  #define B( k, r )  load_8x1( &b[ k + ( r + r0 ) * sr ], v##r );

  #define C( k, v )  store_8x1( v, &a[ k ] )

  #define O1( k )    A( k )                        \
                     B( k, 1 )                     \
                     C( k, v0 + v1 )

  #define O2( k )    A( k )                        \
                     B( k, 1 )                     \
                     B( k, 2 )                     \
                     C( k, ( v0 + v1 ) +           \
                             v2 )

  #define O3( k )    A( k )                        \
                     B( k, 1 )                     \
                     B( k, 2 )                     \
                     B( k, 3 )                     \
                     C( k, ( v0 + v1 ) +           \
                           ( v2 + v3 ) )

  #define O4( k )    A( k )                        \
                     B( k, 1 )                     \
                     B( k, 2 )                     \
                     B( k, 3 )                     \
                     B( k, 4 )                     \
                     C( k, ( ( v0 + v1 ) +         \
                             ( v2 + v3 ) ) +       \
                               v4 )

  #define O5( k )    A( k )                        \
                     B( k, 1 )                     \
                     B( k, 2 )                     \
                     B( k, 3 )                     \
                     B( k, 4 )                     \
                     B( k, 5 )                     \
                     C( k, ( ( v0 + v1 ) +         \
                             ( v2 + v3 ) ) +       \
                             ( v4 + v5 ) )

  #define O6( k )    A( k )                        \
                     B( k, 1 )                     \
                     B( k, 2 )                     \
                     B( k, 3 )                     \
                     B( k, 4 )                     \
                     B( k, 5 )                     \
                     B( k, 6 )                     \
                     C( k, ( ( v0 + v1 ) +         \
                             ( v2 + v3 ) ) +       \
                           ( ( v4 + v5 ) +         \
                               v6 ) )

  #define O7( k )    A( k )                        \
                     B( k, 1 )                     \
                     B( k, 2 )                     \
                     B( k, 3 )                     \
                     B( k, 4 )                     \
                     B( k, 5 )                     \
                     B( k, 6 )                     \
                     B( k, 7 )                     \
                     C( k, ( ( v0 + v1 ) +         \
                             ( v2 + v3 ) ) +       \
                           ( ( v4 + v5 ) +         \
                             ( v6 + v7 ) ) )

  #define O8( k )    A( k )                        \
                     B( k, 1 )                     \
                     B( k, 2 )                     \
                     B( k, 3 )                     \
                     B( k, 4 )                     \
                     B( k, 5 )                     \
                     B( k, 6 )                     \
                     B( k, 7 )                     \
                     B( k, 8 )                     \
                     C( k, ( ( ( v0 + v1 ) +       \
                               ( v2 + v3 ) ) +     \
                             ( ( v4 + v5 ) +       \
                               ( v6 + v7 ) ) ) +   \
                                 v8 )

  #define O9( k )    A( k )                        \
                     B( k, 1 )                     \
                     B( k, 2 )                     \
                     B( k, 3 )                     \
                     B( k, 4 )                     \
                     B( k, 5 )                     \
                     B( k, 6 )                     \
                     B( k, 7 )                     \
                     B( k, 8 )                     \
                     B( k, 9 )                     \
                     C( k, ( ( ( v0 + v1 ) +       \
                               ( v2 + v3 ) ) +     \
                             ( ( v4 + v5 ) +       \
                               ( v6 + v7 ) ) ) +   \
                               ( v8 + v9 ) )

  r0 = -1;

  while( nr )
  {
    switch( nr )
    {
      case  0:                                 break;
      case  1: LOOP( O1 ); nr -= 1 ; r0 += 1 ; break;
      case  2: LOOP( O2 ); nr -= 2 ; r0 += 2 ; break;
      case  3: LOOP( O3 ); nr -= 3 ; r0 += 3 ; break;
      case  4: LOOP( O4 ); nr -= 4 ; r0 += 4 ; break;
      case  5: LOOP( O5 ); nr -= 5 ; r0 += 5 ; break;
      case  6: LOOP( O6 ); nr -= 6 ; r0 += 6 ; break;
      case  7: LOOP( O7 ); nr -= 7 ; r0 += 7 ; break;
      case  8: LOOP( O8 ); nr -= 8 ; r0 += 8 ; break;
      default: LOOP( O9 ); nr -= 9 ; r0 += 9 ; break;
    }
  }

  #undef O9
  #undef O8
  #undef O7
  #undef O6
  #undef O5
  #undef O4
  #undef O3
  #undef O2
  #undef O1
  #undef C
  #undef B
  #undef A
  #undef LOOP
}

#else

void
reduce_array_pipeline_v8( reduce_pipeline_args_t * args,
                          int pipeline_rank,
                          int n_pipeline )
{
  // No v8 implementation.
  ERROR( ( "No reduce_array_pipeline_v8 implementation." ) );
}

#endif
