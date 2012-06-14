#include <rng.h>
#include <util.h>
#include <stdio.h>

#define N 10000000

int
main( int argc,
      char ** argv ) {
  boot_services( &argc, &argv );

# if 0

  /* Output data suitable for testing against SFMT test vectors */

  do {
    rng_t * r = new_rng( 1234 );
    int n;
    for( n=0; n<1000; n++ ) {
      printf( "%10u ", uirand( r ) );
      if( n%5==4 ) printf( "\n" );
    }
    delete_rng( r );
  } while(0);

# endif


# if 0

# define TEST( type, gen, seed, N )                                     \
  do {                                                                  \
    FILE * f;                                                           \
    rng_t * r;                                                          \
    type * x;                                                           \
    int n;                                                              \
    r = new_rng( seed );                                                \
    MALLOC( x, N );                                                     \
    gen##_fill( r, x, 1, N );                                           \
    seed_rng( r, seed );                                                \
    for( n=0; n<N; n++ ) if( x[n]!=gen(r) ) break;                      \
    if( n!=N ) ERROR(( #gen " failed" ));                               \
    f = fopen( #gen ".bin", "wb" );                                     \
    if( !f ) ERROR(( #gen ".bin fopen failed" ));                       \
    if( fwrite( x, sizeof(*x)*N, 1, f )!=1 )                            \
      ERROR(( #gen ".bin fwrite failed" ));                             \
    if( fclose( f ) ) ERROR(( #gen ".bin fclose failed" ));             \
    FREE( x );                                                          \
    delete_rng( r );                                                    \
  } while(0)

  /* Test integer generators */

  TEST( char,      crand,  0, N ); TEST( unsigned char,   ucrand,  1, N );
  TEST( short,     hrand,  2, N ); TEST( unsigned short,  uhrand,  3, N );
  TEST( int,       irand,  4, N ); TEST( unsigned int,    uirand,  5, N );
  TEST( long,      lrand,  6, N ); TEST( unsigned long,   ulrand,  7, N );
  TEST( int8_t,   i8rand,  8, N ); TEST( uint8_t,         u8rand,  9, N );
  TEST( int16_t, i16rand, 10, N ); TEST( uint16_t,       u16rand, 11, N );
  TEST( int32_t, i32rand, 12, N ); TEST( uint32_t,       u32rand, 13, N );
  TEST( int64_t, i64rand, 14, N ); TEST( uint64_t,       u64rand, 14, N );

  /* Test the uniform floating point generators */

  TEST( float, frand,    16, N ); TEST( double, drand,    17, N );
  TEST( float, frand_c0, 18, N ); TEST( double, drand_c0, 19, N );
  TEST( float, frand_c1, 20, N ); TEST( double, drand_c1, 21, N );
  TEST( float, frand_c,  22, N ); TEST( double, drand_c,  23, N );

  /* Test the non-uniform floating point generators */

  TEST( float, frandn,   24, N ); TEST( double, drandn,   25, N );
  TEST( float, frande,   26, N ); TEST( double, drande,   27, N );

# undef TEST

# endif

  /* Test randperm */

  do {
    rng_t * r;
    int * x, n;
    unsigned char * h;
    r = new_rng( 28 );
    MALLOC( x, N ); randperm( r, x, N );
    MALLOC( h, N ); CLEAR( h, N );
    for( n=0; n<N; n++ ) h[x[n]]++;
    for( n=0; n<N; n++ ) if( h[n]!=1 ) break;
    if( n!=N ) ERROR(( "randperm failed" ));
    FREE( x );
    FREE( h );
    delete_rng( r );
  } while(0);

  /* Test shuffle */

# define TEST( type, seed, N )                  \
  do {                                          \
    rng_t * r;                                  \
    type * x, n;                                \
    unsigned char * h;                          \
    r = new_rng( seed );                        \
    MALLOC( x, N );                             \
    for( n=0; n<N; n++ ) x[(int)n] = n;         \
    shuffle( r, x, sizeof(*x), sizeof(*x), N ); \
    MALLOC( h, N ); CLEAR( h, N );              \
    for( n=0; n<N; n++ ) h[(int)x[(int)n]]++;      \
    for( n=0; n<N; n++ ) if( h[(int)n]!=1 ) break; \
    if( n!=N ) ERROR(( "shuffle failed" ));     \
    FREE( x );                                  \
    FREE( h );                                  \
    delete_rng( r );                            \
  } while(0)
    
  TEST( uint8_t,     29,   255 );
  TEST( uint16_t,    30, 65535 );
  TEST( uint32_t,    31,     N );
  TEST( uint64_t,    32,     N );
  TEST( long double, 33,     N );

# undef TEST
    
  halt_services();
  return 0;
}

