#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define HALF (0.5l)
#define TWO  (2.0l)

long double
f( long double x ) {
  return expl( -HALF*x*x );
}

long double
inv_f( long double y ) {
  if( y<=0 || y>=1 ) return 0;
  return sqrt(-TWO*logl(y));
}

/* Integral of tail from r to infinity */

long double
int_g( long double r ) {
  return f(r)/r;
}

long double
make_zig_aux( long double * x,
              long double * y,
              int N,
              long double r ) {
  long double v = int_g(r) + r*f(r);
  int n;
  /**/                   x[N  ] = v/f(r),                     y[N]   = f(x[N]  );
  /**/                   x[N-1] = r,                          y[N-1] = f(x[N-1]);
  for( n=N-2; n>0; n-- ) x[n  ] = inv_f( y[n+1] + v/x[n+1] ), y[n]   = f(x[n]  );
  /**/                   x[0  ] = 0,                          y[0]   = f(x[0]  );
  return v - ( x[1]-x[0] )*( y[0]-y[1] );
}

long double
make_zig( long double * x,
          long double * y,
          int N ) {
  long double a = 0, b = 10, r, dv;
  for(;;) {
    r = 0.5*(a+b);
    if( r==a || r==b ) break;
    dv = make_zig_aux( x, y, N, r );
    if( dv==0 ) break;
    if( dv>0 ) a = r;
    else       b = r;
  }

  return r;
}

int
main( int argc,
      char ** argv ) {
  
  int n, N = atoi( argv[1] );
  long double * x = malloc( (N+1)*sizeof(long double) );
  long double * y = malloc( (N+1)*sizeof(long double) );
  long double R = make_zig( x, y, N );

  printf(   "  static const int   N     = %i;\n"
            "  static const float R     = %.40Lef;\n"
            "  static const float scale = 1.f/%.40Lef;\n",
            N, R, exp2l( atof( argv[2] ) ) );
  printf(   "\n" );
  printf(   "  static const float zig_x[%i] = {\n", N+1 );
  for( n=0; n<N; n++ )
    printf( "    %.40Lef,   // %i\n", x[n], n );
  printf(   "    %.40Lef }; // %i\n", x[n], n );
  printf(   "\n" );
  printf(   "  static const float zig_y[%i] = {\n", N+1 );
  for( n=0; n<N; n++ )
    printf( "    %.40Lef,   // %i\n", y[n], n );
  printf(   "    %.40Lef }; // %i\n", y[n], n );

  return 0;
}

