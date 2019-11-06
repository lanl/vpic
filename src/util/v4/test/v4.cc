/*~--------------------------------------------------------------------------~*
 *~--------------------------------------------------------------------------~*/

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"

#include "../../util.h"
#include "../v4.h"

using namespace v4;

TEST_CASE("any", "[v4]") {
  v4int a;
  int i;

  for( i=0; i<16; i++ ) {
    a[0] = i&1, a[1] = i&2, a[2] = i&4, a[3] = i&8;
    //ASSERT_FALSE( ( i>0 && !any(a) ) || ( i==0 && any(a) ) );
    if ( i > 0 )
    {
      REQUIRE( any(a) );
    }
    else if ( i == 0 )
    {
      REQUIRE_FALSE( any(a) );
    }
  } // for
} // TEST_CASE


TEST_CASE("TEST_CASE_all", "[v4]") {
  v4int a;
  int i;
  for( i=0; i<16; i++ ) {
    a[0] = i&1, a[1] = i&2, a[2] = i&4, a[3] = i&8;
    //ASSERT_FALSE( ( i<15 && all(a) ) || ( i==15 && !all(a) ) );
    if ( i < 15 ) {
      REQUIRE_FALSE( all(a) );
    }
    else if ( i == 15 ) {
      REQUIRE_FALSE( !all(a) );
    }
  } // for
} // TEST_CASE

TEST_CASE("TEST_CASE_splat", "[v4]") {
  v4int a( 1, 2, 3, 4), b( 5, 6, 7, 8), c( 9,10,11,12), d(13,14,15,16);
  v4int e(17,18,19,20);
  b = splat<0>(a); c = splat<1>(a); d = splat<2>(a); e = splat<3>(a);
  //ASSERT_FALSE( any(a!=v4int(1,2,3,4)) || any(b!=v4int(1,1,1,1)) ||
  //any(c!=v4int(2,2,2,2)) || any(d!=v4int(3,3,3,3)) ||
  //any(e!=v4int(4,4,4,4)) );
  REQUIRE( any( a==v4int(1,2,3,4) ) );
  REQUIRE( any( b==v4int(1,1,1,1) ) );
  REQUIRE( any( c==v4int(2,2,2,2) ) );
  REQUIRE( any( d==v4int(3,3,3,3) ) );
  REQUIRE( any( e==v4int(4,4,4,4) ) );
} // TEST_CASE

TEST_CASE("TEST_CASE_shuffle", "[v4]") {
  v4int a( 0, 1, 2, 3), b( 4, 8,12,16), c( 5, 9,13,17), d( 6,10,14,18);
  v4int e( 7,11,15,19);
  b = shuffle<1,2,3,0>(a);
  c = shuffle<2,3,0,1>(a);
  d = shuffle<3,0,1,2>(a);
  e = shuffle<3,2,1,0>(a);
  //ASSERT_FALSE( any(a!=v4int(0,1,2,3)) || any(b!=v4int(1,2,3,0)) ||
  //any(c!=v4int(2,3,0,1)) || any(d!=v4int(3,0,1,2)) ||
  //any(e!=v4int(3,2,1,0)) );
  REQUIRE_FALSE( any(a!=v4int(0,1,2,3)) );
  REQUIRE_FALSE( any(b!=v4int(1,2,3,0)) );
  REQUIRE_FALSE( any(c!=v4int(2,3,0,1)) );
  REQUIRE_FALSE( any(d!=v4int(3,0,1,2)) );
  REQUIRE_FALSE( any(e!=v4int(3,2,1,0)) );
} // TEST_CASE

TEST_CASE("TEST_CASE_swap", "[v4]") {
  v4int a(1,2,3,4), b(5,6,7,8);
  swap(a,b);
  //ASSERT_FALSE( any(a!=v4int(5,6,7,8)) || any(b!=v4int(1,2,3,4)) );
  REQUIRE( any( a==v4int(5,6,7,8) ) );
  REQUIRE( any( b=v4int(1,2,3,4) ) );
} // TEST_CASE

TEST_CASE("TEST_CASE_transpose", "[v4]") {
  v4int a0( 0, 1, 2, 3), a1( 4, 5, 6, 7), a2( 8, 9,10,11), a3(12,13,14,15);
  transpose(a0,a1,a2,a3);
  //ASSERT_FALSE( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
  //any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) );
  REQUIRE( any(a0=v4int( 0, 4, 8,12)) );
  REQUIRE( any(a1==v4int( 1, 5, 9,13)) );
  REQUIRE( any(a2==v4int( 2, 6,10,14)) );
  REQUIRE( any(a3==v4int( 3, 7,11,15)) );
} // TEST_CASE

TEST_CASE("TEST_CASE_load_4x1", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0(1,0,0,0), a1(0,0,0,0), a2(0,0,0,0), a3(0,0,0,0);
  int i;
  for( i=0; i<16; i++ ) {
    mem[i] = i;
  }

  load_4x1( mem,    a0 );
  load_4x1( mem+4,  a1 );
  load_4x1( mem+8,  a2 );
  load_4x1( mem+12, a3 );

  for( i=0; i<16; i++ ) {
    if( mem[i]!=i ) break;
  }

  //ASSERT_FALSE( any(a0!=v4int( 0, 1, 2, 3)) || any(a1!=v4int( 4, 5, 6, 7)) ||
  //any(a2!=v4int( 8, 9,10,11)) || any(a3!=v4int(12,13,14,15)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 1, 2, 3)) );
  REQUIRE( any(a1==v4int( 4, 5, 6, 7)) );
  REQUIRE( any(a2==v4int( 8, 9,10,11)) );
  REQUIRE( any(a3==v4int(12,13,14,15)) );
  REQUIRE( i==16 );

} // TEST_CASE

TEST_CASE("TEST_CASE_store_4x1", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 1, 2, 3), a1( 4, 5, 6, 7), a2( 8, 9,10,11), a3(12,13,14,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x1(a0,mem);
  store_4x1(a1,mem+4);
  store_4x1(a2,mem+8);
  store_4x1(a3,mem+12);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;

  //ASSERT_FALSE( any(a0!=v4int( 0, 1, 2, 3)) || any(a1!=v4int( 4, 5, 6, 7)) ||
  //any(a2!=v4int( 8, 9,10,11)) || any(a3!=v4int(12,13,14,15)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 1, 2, 3)) );
  REQUIRE( any(a1==v4int( 4, 5, 6, 7)) );
  REQUIRE( any(a2==v4int( 8, 9,10,11)) );
  REQUIRE( any(a3==v4int(12,13,14,15)) );
  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_stream_4x1", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 1, 2, 3), a1( 4, 5, 6, 7), a2( 8, 9,10,11), a3(12,13,14,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  stream_4x1(a0,mem);
  stream_4x1(a1,mem+4);
  stream_4x1(a2,mem+8);
  stream_4x1(a3,mem+12);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;

  //ASSERT_FALSE( any(a0!=v4int( 0, 1, 2, 3)) || any(a1!=v4int( 4, 5, 6, 7)) ||
  //any(a2!=v4int( 8, 9,10,11)) || any(a3!=v4int(12,13,14,15)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 1, 2, 3)) );
  REQUIRE( any(a1==v4int( 4, 5, 6, 7)) );
  REQUIRE( any(a2==v4int( 8, 9,10,11)) );
  REQUIRE( any(a3==v4int(12,13,14,15)) );
  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_clear_4x1", "[v4]") {
  v4float vmem[4]; float * mem = (float *)vmem;
  int i;
  for(i=0; i<16; i++) mem[i] = i;
  clear_4x1( mem+8  );
  clear_4x1( mem+12 );
  for(i=8; i<16; i++) mem[i] += i;
  for(i=0; i<16; i++) if( mem[i]!=i ) break;

  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_copy_4x1", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  int i;
  for( i=0; i<8; i++ ) mem[i] = i;
  copy_4x1(mem+8, mem  );
  copy_4x1(mem+12,mem+4);
  for( i=8; i<16; i++ ) mem[i] += 8;
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;

  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_swap_4x1", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  int i;
  for( i=0; i<8; i++ ) mem[i] = i;
  copy_4x1(mem+12, mem  );
  copy_4x1(mem+8,  mem+4);
  for( i=8; i<16; i++ ) mem[i] += 8;
  swap_4x1(mem+8,  mem+12 );
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_load_4x1_tr", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0, a1, a2, a3;
  int i;
  for( i=0; i<16; i++ ) mem[i] = i;
  load_4x1_tr(mem,  mem+4,mem+8, mem+12,a0);
  load_4x1_tr(mem+1,mem+5,mem+9, mem+13,a1);
  load_4x1_tr(mem+2,mem+6,mem+10,mem+14,a2);
  load_4x1_tr(mem+3,mem+7,mem+11,mem+15,a3);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;

  //ASSERT_FALSE( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
  //any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) ||
  //i!=16 );

  REQUIRE( any(a0==v4int( 0, 4, 8,12)) );
  REQUIRE( any(a1==v4int( 1, 5, 9,13)) );
  REQUIRE( any(a2==v4int( 2, 6,10,14)) );
  REQUIRE( any(a3==v4int( 3, 7,11,15)) );
  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_load_4x2_tr", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 32 );
  v4int a0, a1, a2, a3;
  int i;
  for( i=0; i<32; i++ ) mem[i] = i;
  load_4x2_tr(mem,   mem+4, mem+8, mem+12,a0,a1);
  load_4x2_tr(mem+16,mem+20,mem+24,mem+28,a2,a3);
  for( i=0; i<32; i++ ) if( mem[i]!=i ) break;

  //ASSERT_FALSE( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
  //any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 4, 8,12)) );
  REQUIRE( any(a1==v4int( 1, 5, 9,13)) );
  REQUIRE( any(a2==v4int(16,20,24,28)) );
  REQUIRE( any(a3==v4int(17,21,25,29)) );
  REQUIRE( i==32 );

} // TEST_CASE

TEST_CASE("TEST_CASE_load_4x3_tr", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0, a1, a2;
  int i;
  for( i=0; i<16; i++ ) mem[i] = i;
  load_4x3_tr(mem,mem+4,mem+8,mem+12,a0,a1,a2);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  //ASSERT_FALSE( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
  //any(a2!=v4int( 2, 6,10,14)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 4, 8,12)) );
  REQUIRE( any(a1==v4int( 1, 5, 9,13)) );
  REQUIRE( any(a2==v4int( 2, 6,10,14)) );
  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_load_4x4_tr", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0, a1, a2, a3;
  int i;
  for( i=0; i<16; i++ ) mem[i] = i;
  load_4x4_tr(mem,mem+4,mem+8,mem+12,a0,a1,a2,a3);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  //ASSERT_FALSE( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
  //any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 4, 8,12)) );
  REQUIRE( any(a1==v4int( 1, 5, 9,13)) );
  REQUIRE( any(a2==v4int( 2, 6,10,14)) );
  REQUIRE( any(a3==v4int( 3, 7,11,15)) );
  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_store_4x1_tr", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 4, 8,12), a1( 1, 5, 9,13), a2( 2, 6,10,14), a3( 3, 7,11,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x1_tr(a0,mem,  mem+4,mem+8, mem+12);
  store_4x1_tr(a1,mem+1,mem+5,mem+9, mem+13);
  store_4x1_tr(a2,mem+2,mem+6,mem+10,mem+14);
  store_4x1_tr(a3,mem+3,mem+7,mem+11,mem+15);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  //ASSERT_FALSE( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
  //any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 4, 8,12)) );
  REQUIRE( any(a1==v4int( 1, 5, 9,13)) );
  REQUIRE( any(a2==v4int( 2, 6,10,14)) );
  REQUIRE( any(a3==v4int( 3, 7,11,15)) );
  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_store_4x2_tr", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 4, 8,12), a1( 1, 5, 9,13), a2( 2, 6,10,14), a3( 3, 7,11,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x2_tr(a0,a1,mem,  mem+4,mem+8, mem+12);
  store_4x2_tr(a2,a3,mem+2,mem+6,mem+10,mem+14);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  //ASSERT_FALSE( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
  //any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 4, 8,12)) );
  REQUIRE( any(a1==v4int( 1, 5, 9,13)) );
  REQUIRE( any(a2==v4int( 2, 6,10,14)) );
  REQUIRE( any(a3==v4int( 3, 7,11,15)) );
  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_store_4x3_tr", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 4, 8,12), a1( 1, 5, 9,13), a2( 2, 6,10,14), a3( 3, 7,11,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x3_tr(a0,a1,a2,mem,  mem+4,mem+8, mem+12);
  for( i=0; i<16; i++ )
    if( ( (i&3)!=3 && mem[i]!=i ) || ( (i&3)==3 && mem[i]!=0 ) ) break;
  //ASSERT_FALSE( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
  //any(a2!=v4int( 2, 6,10,14)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 4, 8,12)) );
  REQUIRE( any(a1==v4int( 1, 5, 9,13)) );
  REQUIRE( any(a2==v4int( 2, 6,10,14)) );
  REQUIRE( i==16 );
} // TEST_CASE

TEST_CASE("TEST_CASE_store_4x4_tr", "[v4]") {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 4, 8,12), a1( 1, 5, 9,13), a2( 2, 6,10,14), a3( 3, 7,11,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x4_tr(a0,a1,a2,a3,mem,  mem+4,mem+8, mem+12);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  //ASSERT_FALSE( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
  //any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 );

  REQUIRE( any(a0==v4int( 0, 4, 8,12)) );
  REQUIRE( any(a1==v4int( 1, 5, 9,13)) );
  REQUIRE( any(a2==v4int( 2, 6,10,14)) );
  REQUIRE( any(a3==v4int( 3, 7,11,15)) );
  REQUIRE( i==16 );
} // TEST_CASE

/*~-------------------------------------------------------------------------~-*
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
