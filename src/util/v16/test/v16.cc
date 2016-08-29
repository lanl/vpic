#include <gtest/gtest.h>

#include "../../util.h"
#include "../v16.h"

#include <iostream>

using namespace v16;

TEST(v16, test_transpose)
{
  v16int a00(  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15);
  v16int a01( 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31);
  v16int a02( 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47);
  v16int a03( 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63);
  v16int a04( 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79);
  v16int a05( 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95);
  v16int a06( 96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111);
  v16int a07(112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127);
  v16int a08(128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143);
  v16int a09(144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159);
  v16int a10(160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175);
  v16int a11(176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191);
  v16int a12(192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207);
  v16int a13(208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223);
  v16int a14(224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239);
  v16int a15(240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255);

  transpose( a00, a01, a02, a03, a04, a05, a06, a07, a08, a09, a10, a11, a12, a13, a14, a15 );

  ASSERT_FALSE( any( a00 != v16int(  0, 16, 32, 48, 64, 80,  96, 112, 128, 144, 160, 176, 192, 208, 224, 240 ) ) ||
		any( a01 != v16int(  1, 17, 33, 49, 65, 81,  97, 113, 129, 145, 161, 177, 193, 209, 225, 241 ) ) ||
		any( a02 != v16int(  2, 18, 34, 50, 66, 82,  98, 114, 130, 146, 162, 178, 194, 210, 226, 242 ) ) ||
		any( a03 != v16int(  3, 19, 35, 51, 67, 83,  99, 115, 131, 147, 163, 179, 195, 211, 227, 243 ) ) ||
		any( a04 != v16int(  4, 20, 36, 52, 68, 84, 100, 116, 132, 148, 164, 180, 196, 212, 228, 244 ) ) ||
		any( a05 != v16int(  5, 21, 37, 53, 69, 85, 101, 117, 133, 149, 165, 181, 197, 213, 229, 245 ) ) ||
		any( a06 != v16int(  6, 22, 38, 54, 70, 86, 102, 118, 134, 150, 166, 182, 198, 214, 230, 246 ) ) ||
		any( a07 != v16int(  7, 23, 39, 55, 71, 87, 103, 119, 135, 151, 167, 183, 199, 215, 231, 247 ) ) ||
		any( a08 != v16int(  8, 24, 40, 56, 72, 88, 104, 120, 136, 152, 168, 184, 200, 216, 232, 248 ) ) ||
		any( a09 != v16int(  9, 25, 41, 57, 73, 89, 105, 121, 137, 153, 169, 185, 201, 217, 233, 249 ) ) ||
		any( a10 != v16int( 10, 26, 42, 58, 74, 90, 106, 122, 138, 154, 170, 186, 202, 218, 234, 250 ) ) ||
		any( a11 != v16int( 11, 27, 43, 59, 75, 91, 107, 123, 139, 155, 171, 187, 203, 219, 235, 251 ) ) ||
		any( a12 != v16int( 12, 28, 44, 60, 76, 92, 108, 124, 140, 156, 172, 188, 204, 220, 236, 252 ) ) ||
		any( a13 != v16int( 13, 29, 45, 61, 77, 93, 109, 125, 141, 157, 173, 189, 205, 221, 237, 253 ) ) ||
		any( a14 != v16int( 14, 30, 46, 62, 78, 94, 110, 126, 142, 158, 174, 190, 206, 222, 238, 254 ) ) ||
		any( a15 != v16int( 15, 31, 47, 63, 79, 95, 111, 127, 143, 159, 175, 191, 207, 223, 239, 255 ) ) );
}

#if 0

TEST(v16, test_any)
{
  v16int a;
  int i;

  for( i=0; i < 256; i++ )
  {
    a[0] = i&1,  a[1] = i&2,  a[2] = i&4,  a[3] = i&8;
    a[4] = i&16, a[5] = i&32, a[6] = i&64, a[7] = i&128;

    ASSERT_FALSE( ( i>0 && !any(a) ) || ( i==0 && any(a) ) );
  }
}

TEST(v16, test_all)
{
  v16int a;
  int i;
  for( i=0; i < 256; i++ )
  {
    a[0] = i&1,  a[1] = i&2,  a[2] = i&4,  a[3] = i&8;
    a[4] = i&16, a[5] = i&32, a[6] = i&64, a[7] = i&128;

    ASSERT_FALSE( ( i < 255 && all(a) ) || ( i == 255 && !all(a) ) );
  }
}

TEST(v16, test_splat)
{
  v16int a( 1, 2, 3, 4, 5, 6, 7, 8);
  v16int b( 9,10,11,12,13,14,15,16);
  v16int c(17,18,19,20,21,22,23,24);
  v16int d(25,26,27,28,29,30,31,32);
  v16int e(33,34,35,36,37,38,39,40);
  v16int f(41,42,43,44,45,46,47,48);
  v16int g(49,50,51,52,53,54,55,56);
  v16int h(57,58,59,60,61,62,63,64);
  v16int i(65,66,67,68,69,70,71,72);

  b = splat<0>(a);
  c = splat<1>(a);
  d = splat<2>(a);
  e = splat<3>(a);
  f = splat<4>(a);
  g = splat<5>(a);
  h = splat<6>(a);
  i = splat<7>(a);

  ASSERT_FALSE( any(a!=v16int(1,2,3,4,5,6,7,8)) ||
		any(b!=v16int(1,1,1,1,1,1,1,1)) ||
		any(c!=v16int(2,2,2,2,2,2,2,2)) ||
		any(d!=v16int(3,3,3,3,3,3,3,3)) ||
		any(e!=v16int(4,4,4,4,4,4,4,4)) ||
		any(f!=v16int(5,5,5,5,5,5,5,5)) ||
		any(g!=v16int(6,6,6,6,6,6,6,6)) ||
		any(h!=v16int(7,7,7,7,7,7,7,7)) ||
		any(i!=v16int(8,8,8,8,8,8,8,8)) );
}

TEST(v16, test_shuffle)
{
  v16int a( 0, 1, 2, 3, 4, 5, 6, 7);
  v16int b( 9,10,11,12,13,14,15,16);
  v16int c(17,18,19,20,21,22,23,24);
  v16int d(25,26,27,28,29,30,31,32);
  v16int e(33,34,35,36,37,38,39,40);
  v16int f(41,42,43,44,45,46,47,48);
  v16int g(49,50,51,52,53,54,55,56);
  v16int h(57,58,59,60,61,62,63,64);
  v16int i(65,66,67,68,69,70,71,72);

  b = shuffle<1,2,3,4,5,6,7,0>(a);
  c = shuffle<2,3,4,5,6,7,0,1>(a);
  d = shuffle<3,4,5,6,7,0,1,2>(a);
  e = shuffle<4,5,6,7,0,1,2,3>(a);
  f = shuffle<5,6,7,0,1,2,3,4>(a);
  g = shuffle<6,7,0,1,2,3,4,5>(a);
  h = shuffle<7,0,1,2,3,4,5,6>(a);
  i = shuffle<7,6,5,4,3,2,1,0>(a);

  ASSERT_FALSE( any(a!=v16int(0,1,2,3,4,5,6,7)) ||
		any(b!=v16int(1,2,3,4,5,6,7,0)) ||
		any(c!=v16int(2,3,4,5,6,7,0,1)) ||
		any(d!=v16int(3,4,5,6,7,0,1,2)) ||
		any(e!=v16int(4,5,6,7,0,1,2,3)) ||
		any(f!=v16int(5,6,7,0,1,2,3,4)) ||
		any(g!=v16int(6,7,0,1,2,3,4,5)) ||
		any(h!=v16int(7,0,1,2,3,4,5,6)) ||
		any(i!=v16int(7,6,5,4,3,2,1,0)) );
}

// #endif

TEST(v16, test_swap)
{
  v16int a( 1, 2, 3, 4, 5, 6, 7, 8);
  v16int b( 9,10,11,12,13,14,15,16);

  swap(a,b);

  ASSERT_FALSE( any( a != v16int( 9,10,11,12,13,14,15,16) ) ||
		any( b != v16int( 1, 2, 3, 4, 5, 6, 7, 8) ) );
}

TEST(v16, test_load_8x1)
{
  DECLARE_ALIGNED_ARRAY( int, 32, mem, 32 );

  v16int a0(1,0,0,0,0,0,0,0);
  v16int a1(0,0,0,0,0,0,0,0);
  v16int a2(0,0,0,0,0,0,0,0);
  v16int a3(0,0,0,0,0,0,0,0);

  int i;

  for( i=0; i < 32; i++ ) mem[i] = i;

  load_8x1( mem,    a0 );
  load_8x1( mem+8,  a1 );
  load_8x1( mem+16, a2 );
  load_8x1( mem+24, a3 );

  for( i=0; i < 32; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int( 0, 1, 2, 3, 4, 5, 6, 7) ) ||
		any( a1 != v16int( 8, 9,10,11,12,13,14,15) ) ||
		any( a2 != v16int(16,17,18,19,20,21,22,23) ) ||
		any( a3 != v16int(24,25,26,27,28,29,30,31) ) ||
		i != 32 );
}

TEST(v16, test_store_8x1)
{
  DECLARE_ALIGNED_ARRAY( int, 32, mem, 32 );

  v16int a0( 0, 1, 2, 3, 4, 5, 6, 7);
  v16int a1( 8, 9,10,11,12,13,14,15);
  v16int a2(16,17,18,19,20,21,22,23);
  v16int a3(24,25,26,27,28,29,30,31);

  int i;

  for( i=0; i < 32; i++ ) mem[i] = 0;

  store_8x1( a0, mem      );
  store_8x1( a1, mem +  8 );
  store_8x1( a2, mem + 16 );
  store_8x1( a3, mem + 24 );

  for( i=0; i < 32; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int( 0, 1, 2, 3, 4, 5, 6, 7) ) ||
		any( a1 != v16int( 8, 9,10,11,12,13,14,15) ) ||
		any( a2 != v16int(16,17,18,19,20,21,22,23) ) ||
		any( a3 != v16int(24,25,26,27,28,29,30,31) ) ||
		i != 32 );
}

TEST(v16, test_stream_8x1)
{
  DECLARE_ALIGNED_ARRAY( int, 32, mem, 32 );

  v16int a0( 0, 1, 2, 3, 4, 5, 6, 7);
  v16int a1( 8, 9,10,11,12,13,14,15);
  v16int a2(16,17,18,19,20,21,22,23);
  v16int a3(24,25,26,27,28,29,30,31);

  int i;

  for( i=0; i < 32; i++ ) mem[i] = 0;

  stream_8x1( a0, mem      );
  stream_8x1( a1, mem +  8 );
  stream_8x1( a2, mem + 16 );
  stream_8x1( a3, mem + 24 );

  for( i=0; i < 32; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int( 0, 1, 2, 3, 4, 5, 6, 7) ) ||
		any( a1 != v16int( 8, 9,10,11,12,13,14,15) ) ||
		any( a2 != v16int(16,17,18,19,20,21,22,23) ) ||
		any( a3 != v16int(24,25,26,27,28,29,30,31) ) ||
		i != 32 );
}

TEST(v16, test_clear_8x1)
{
  v16float vmem[4]; float * mem = (float *)vmem;

  int i;

  for(i=0; i < 32; i++) mem[i] = i;

  clear_8x1( mem + 16 );
  clear_8x1( mem + 24 );

  for(i=16; i < 32; i++) mem[i] += i;

  for(i=0; i < 32; i++) if( mem[i] != i ) break;

  ASSERT_FALSE( i != 32 );
}

TEST(v16, test_copy_8x1)
{
  DECLARE_ALIGNED_ARRAY( int, 32, mem, 32 );

  int i;

  for( i=0; i < 16; i++ ) mem[i] = i;

  copy_8x1( mem + 16, mem     );
  copy_8x1( mem + 24, mem + 8 );

  for( i=16; i < 32; i++ ) mem[i] += 16;

  for( i=0; i < 32; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( i != 32 );
}

TEST(v16, test_swap_8x1)
{
  DECLARE_ALIGNED_ARRAY( int, 32, mem, 32 );

  int i;

  for( i=0; i < 16; i++ ) mem[i] = i;

  copy_8x1( mem + 24, mem     );
  copy_8x1( mem + 16, mem + 8 );

  for( i=16; i < 32; i++ ) mem[i] += 16;

  swap_8x1( mem + 16, mem + 24 );

  for( i=0; i < 32; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( i != 32 );
}

TEST(v16, test_load_8x1_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0, a1, a2, a3, a4, a5, a6, a7;

  int i;

  for( i=0; i < 64; i++ ) mem[i] = i;

  load_8x1_tr( mem,   mem+8,  mem+16, mem+24, mem+32, mem+40, mem+48, mem+56, a0 );
  load_8x1_tr( mem+1, mem+9,  mem+17, mem+25, mem+33, mem+41, mem+49, mem+57, a1 );
  load_8x1_tr( mem+2, mem+10, mem+18, mem+26, mem+34, mem+42, mem+50, mem+58, a2 );
  load_8x1_tr( mem+3, mem+11, mem+19, mem+27, mem+35, mem+43, mem+51, mem+59, a3 );
  load_8x1_tr( mem+4, mem+12, mem+20, mem+28, mem+36, mem+44, mem+52, mem+60, a4 );
  load_8x1_tr( mem+5, mem+13, mem+21, mem+29, mem+37, mem+45, mem+53, mem+61, a5 );
  load_8x1_tr( mem+6, mem+14, mem+22, mem+30, mem+38, mem+46, mem+54, mem+62, a6 );
  load_8x1_tr( mem+7, mem+15, mem+23, mem+31, mem+39, mem+47, mem+55, mem+63, a7 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		any( a3 != v16int(3,11,19,27,35,43,51,59) ) ||
		any( a4 != v16int(4,12,20,28,36,44,52,60) ) ||
		any( a5 != v16int(5,13,21,29,37,45,53,61) ) ||
		any( a6 != v16int(6,14,22,30,38,46,54,62) ) ||
		any( a7 != v16int(7,15,23,31,39,47,55,63) ) ||
		i != 64 );
}

TEST(v16, test_load_8x2_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0, a1, a2, a3, a4, a5, a6, a7;

  int i;

  for( i=0; i < 64; i++ ) mem[i] = i;

  load_8x2_tr( mem,   mem+8,  mem+16, mem+24, mem+32, mem+40, mem+48, mem+56, a0, a1 );
  load_8x2_tr( mem+2, mem+10, mem+18, mem+26, mem+34, mem+42, mem+50, mem+58, a2, a3 );
  load_8x2_tr( mem+4, mem+12, mem+20, mem+28, mem+36, mem+44, mem+52, mem+60, a4, a5 );
  load_8x2_tr( mem+6, mem+14, mem+22, mem+30, mem+38, mem+46, mem+54, mem+62, a6, a7 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		any( a3 != v16int(3,11,19,27,35,43,51,59) ) ||
		any( a4 != v16int(4,12,20,28,36,44,52,60) ) ||
		any( a5 != v16int(5,13,21,29,37,45,53,61) ) ||
		any( a6 != v16int(6,14,22,30,38,46,54,62) ) ||
		any( a7 != v16int(7,15,23,31,39,47,55,63) ) ||
		i != 64 );
}

TEST(v16, test_load_8x2_tr_a)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0, a1, a2, a3, a4, a5, a6, a7;

  int i;

  for( i=0; i < 64; i++ ) mem[i] = i;

  load_8x2_tr( mem,    mem+2,  mem+4,  mem+6,  mem+8,  mem+10, mem+12, mem+14, a0, a1 );
  load_8x2_tr( mem+16, mem+18, mem+20, mem+22, mem+24, mem+26, mem+28, mem+30, a2, a3 );
  load_8x2_tr( mem+32, mem+34, mem+36, mem+38, mem+40, mem+42, mem+44, mem+46, a4, a5 );
  load_8x2_tr( mem+48, mem+50, mem+52, mem+54, mem+56, mem+58, mem+60, mem+62, a6, a7 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int( 0, 2, 4, 6, 8,10,12,14) ) ||
		any( a1 != v16int( 1, 3, 5, 7, 9,11,13,15) ) ||
		any( a2 != v16int(16,18,20,22,24,26,28,30) ) ||
		any( a3 != v16int(17,19,21,23,25,27,29,31) ) ||
		any( a4 != v16int(32,34,36,38,40,42,44,46) ) ||
		any( a5 != v16int(33,35,37,39,41,43,45,47) ) ||
		any( a6 != v16int(48,50,52,54,56,58,60,62) ) ||
		any( a7 != v16int(49,51,53,55,57,59,61,63) ) ||
		i != 64 );
}

TEST(v16, test_load_8x3_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0, a1, a2;

  int i;

  for( i=0; i < 64; i++ ) mem[i] = i;

  load_8x3_tr( mem, mem+8, mem+16, mem+24, mem+32, mem+40, mem+48, mem+56,
	       a0, a1, a2 );

  for( i=0; i < 64; i++ ) if( mem[i]!=i ) break;

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		i != 64 );
}

TEST(v16, test_load_8x4_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0, a1, a2, a3;

  int i;

  for( i=0; i < 64; i++ ) mem[i] = i;

  load_8x4_tr( mem, mem+8, mem+16, mem+24, mem+32, mem+40, mem+48, mem+56,
	       a0, a1, a2, a3 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		any( a3 != v16int(3,11,19,27,35,43,51,59) ) ||
		i != 64 );
}

TEST(v16, test_load_8x4_tr_a)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0, a1, a2, a3;

  int i;

  for( i=0; i < 64; i++ ) mem[i] = i;

  load_8x4_tr( mem, mem+4, mem+8, mem+12, mem+16, mem+20, mem+24, mem+28,
	       a0, a1, a2, a3 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int(0,4, 8,12,16,20,24,28) ) ||
		any( a1 != v16int(1,5, 9,13,17,21,25,29) ) ||
		any( a2 != v16int(2,6,10,14,18,22,26,30) ) ||
		any( a3 != v16int(3,7,11,15,19,23,27,31) ) ||
		i != 64 );
}

TEST(v16, test_load_8x8_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0, a1, a2, a3, a4, a5, a6, a7;

  int i;

  for( i=0; i < 64; i++ ) mem[i] = i;

  load_8x8_tr( mem, mem+8, mem+16, mem+24, mem+32, mem+40, mem+48, mem+56,
	       a0, a1, a2, a3, a4, a5, a6, a7 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		any( a3 != v16int(3,11,19,27,35,43,51,59) ) ||
		any( a4 != v16int(4,12,20,28,36,44,52,60) ) ||
		any( a5 != v16int(5,13,21,29,37,45,53,61) ) ||
		any( a6 != v16int(6,14,22,30,38,46,54,62) ) ||
		any( a7 != v16int(7,15,23,31,39,47,55,63) ) ||
		i != 64 );
}

TEST(v16, test_store_8x1_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0(0, 8,16,24,32,40,48,56);
  v16int a1(1, 9,17,25,33,41,49,57);
  v16int a2(2,10,18,26,34,42,50,58);
  v16int a3(3,11,19,27,35,43,51,59);
  v16int a4(4,12,20,28,36,44,52,60);
  v16int a5(5,13,21,29,37,45,53,61);
  v16int a6(6,14,22,30,38,46,54,62);
  v16int a7(7,15,23,31,39,47,55,63);

  int i;

  for( i=0; i < 64; i++ ) mem[i] = 0;

  store_8x1_tr( a0, mem,   mem+ 8, mem+16, mem+24, mem+32, mem+40, mem+48, mem+56 );
  store_8x1_tr( a1, mem+1, mem+ 9, mem+17, mem+25, mem+33, mem+41, mem+49, mem+57 );
  store_8x1_tr( a2, mem+2, mem+10, mem+18, mem+26, mem+34, mem+42, mem+50, mem+58 );
  store_8x1_tr( a3, mem+3, mem+11, mem+19, mem+27, mem+35, mem+43, mem+51, mem+59 );
  store_8x1_tr( a4, mem+4, mem+12, mem+20, mem+28, mem+36, mem+44, mem+52, mem+60 );
  store_8x1_tr( a5, mem+5, mem+13, mem+21, mem+29, mem+37, mem+45, mem+53, mem+61 );
  store_8x1_tr( a6, mem+6, mem+14, mem+22, mem+30, mem+38, mem+46, mem+54, mem+62 );
  store_8x1_tr( a7, mem+7, mem+15, mem+23, mem+31, mem+39, mem+47, mem+55, mem+63 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		any( a3 != v16int(3,11,19,27,35,43,51,59) ) ||
		any( a4 != v16int(4,12,20,28,36,44,52,60) ) ||
		any( a5 != v16int(5,13,21,29,37,45,53,61) ) ||
		any( a6 != v16int(6,14,22,30,38,46,54,62) ) ||
		any( a7 != v16int(7,15,23,31,39,47,55,63) ) ||
		i != 64 );
}

TEST(v16, test_store_8x2_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0(0, 8,16,24,32,40,48,56);
  v16int a1(1, 9,17,25,33,41,49,57);
  v16int a2(2,10,18,26,34,42,50,58);
  v16int a3(3,11,19,27,35,43,51,59);
  v16int a4(4,12,20,28,36,44,52,60);
  v16int a5(5,13,21,29,37,45,53,61);
  v16int a6(6,14,22,30,38,46,54,62);
  v16int a7(7,15,23,31,39,47,55,63);

  int i;

  for( i=0; i < 64; i++ ) mem[i] = 0;

  store_8x2_tr( a0, a1, mem,   mem+ 8, mem+16, mem+24, mem+32, mem+40, mem+48, mem+56 );
  store_8x2_tr( a2, a3, mem+2, mem+10, mem+18, mem+26, mem+34, mem+42, mem+50, mem+58 );
  store_8x2_tr( a4, a5, mem+4, mem+12, mem+20, mem+28, mem+36, mem+44, mem+52, mem+60 );
  store_8x2_tr( a6, a7, mem+6, mem+14, mem+22, mem+30, mem+38, mem+46, mem+54, mem+62 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		any( a3 != v16int(3,11,19,27,35,43,51,59) ) ||
		any( a4 != v16int(4,12,20,28,36,44,52,60) ) ||
		any( a5 != v16int(5,13,21,29,37,45,53,61) ) ||
		any( a6 != v16int(6,14,22,30,38,46,54,62) ) ||
		any( a7 != v16int(7,15,23,31,39,47,55,63) ) ||
		i != 64 );
}

TEST(v16, test_store_8x2_tr_a)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0( 0, 2, 4, 6, 8,10,12,14);
  v16int a1( 1, 3, 5, 7, 9,11,13,15);
  v16int a2(16,18,20,22,24,26,28,30);
  v16int a3(17,19,21,23,25,27,29,31);
  v16int a4(32,34,36,38,40,42,44,46);
  v16int a5(33,35,37,39,41,43,45,47);
  v16int a6(48,50,52,54,56,58,60,62);
  v16int a7(49,51,53,55,57,59,61,63);

  int i;

  for( i=0; i < 64; i++ ) mem[i] = 0;

  store_8x2_tr( a0, a1, mem,    mem+ 2, mem+ 4, mem+ 6, mem+ 8, mem+10, mem+12, mem+14 );
  store_8x2_tr( a2, a3, mem+16, mem+18, mem+20, mem+22, mem+24, mem+26, mem+28, mem+30 );
  store_8x2_tr( a4, a5, mem+32, mem+34, mem+36, mem+38, mem+40, mem+42, mem+44, mem+46 );
  store_8x2_tr( a6, a7, mem+48, mem+50, mem+52, mem+54, mem+56, mem+58, mem+60, mem+62 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int( 0, 2, 4, 6, 8,10,12,14) ) ||
		any( a1 != v16int( 1, 3, 5, 7, 9,11,13,15) ) ||
		any( a2 != v16int(16,18,20,22,24,26,28,30) ) ||
		any( a3 != v16int(17,19,21,23,25,27,29,31) ) ||
		any( a4 != v16int(32,34,36,38,40,42,44,46) ) ||
		any( a5 != v16int(33,35,37,39,41,43,45,47) ) ||
		any( a6 != v16int(48,50,52,54,56,58,60,62) ) ||
		any( a7 != v16int(49,51,53,55,57,59,61,63) ) ||
		i != 64 );
}

TEST(v16, test_store_8x3_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0(0, 8,16,24,32,40,48,56);
  v16int a1(1, 9,17,25,33,41,49,57);
  v16int a2(2,10,18,26,34,42,50,58);
  v16int a3(3,11,19,27,35,43,51,59);
  v16int a4(4,12,20,28,36,44,52,60);
  v16int a5(5,13,21,29,37,45,53,61);
  v16int a6(6,14,22,30,38,46,54,62);
  v16int a7(7,15,23,31,39,47,55,63);

  int i, j;

  for( i=0; i < 64; i++ ) mem[i] = 0;

  store_8x3_tr( a0, a1, a2,
		mem, mem+8, mem+16, mem+24, mem+32, mem+40, mem+48, mem+56 );

  j = 0;
  for( i=0; i < 8; i++ )
  {
    if( ( i <  3 && mem[i]    != i    ) ||
	( i >= 3 && mem[i]    != 0    ) )
      break;
    else
      j++;

    if( ( i <  3 && mem[i+ 8] != i+ 8 ) ||
	( i >= 3 && mem[i+ 8] != 0   ) )
      break;
    else
      j++;

    if( ( i <  3 && mem[i+16] != i+16 ) ||
	( i >= 3 && mem[i+16] != 0    ) )
      break;
    else
      j++;

    if( ( i <  3 && mem[i+24] != i+24 ) ||
	( i >= 3 && mem[i+24] != 0    ) )
      break;
    else
      j++;

    if( ( i <  3 && mem[i+32] != i+32 ) ||
	( i >= 3 && mem[i+32] != 0    ) )
      break;
    else
      j++;

    if( ( i <  3 && mem[i+40] != i+40 ) ||
	( i >= 3 && mem[i+40] != 0    ) )
      break;
    else
      j++;

    if( ( i <  3 && mem[i+48] != i+48 ) ||
	( i >= 3 && mem[i+48] != 0    ) )
      break;
    else
      j++;

    if( ( i <  3 && mem[i+56] != i+56 ) ||
	( i >= 3 && mem[i+56] != 0    ) )
      break;
    else
      j++;
  }

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		any( a3 != v16int(3,11,19,27,35,43,51,59) ) ||
		any( a4 != v16int(4,12,20,28,36,44,52,60) ) ||
		any( a5 != v16int(5,13,21,29,37,45,53,61) ) ||
		any( a6 != v16int(6,14,22,30,38,46,54,62) ) ||
		any( a7 != v16int(7,15,23,31,39,47,55,63) ) ||
		j != 64 );
}

TEST(v16, test_store_8x4_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0(0, 8,16,24,32,40,48,56);
  v16int a1(1, 9,17,25,33,41,49,57);
  v16int a2(2,10,18,26,34,42,50,58);
  v16int a3(3,11,19,27,35,43,51,59);
  v16int a4(4,12,20,28,36,44,52,60);
  v16int a5(5,13,21,29,37,45,53,61);
  v16int a6(6,14,22,30,38,46,54,62);
  v16int a7(7,15,23,31,39,47,55,63);

  int i;

  for( i=0; i < 64; i++ ) mem[i] = 0;

  store_8x4_tr( a0, a1, a2, a3,
		mem,   mem+ 8, mem+16, mem+24, mem+32, mem+40, mem+48, mem+56 );

  store_8x4_tr( a4, a5, a6, a7,
		mem+4, mem+12, mem+20, mem+28, mem+36, mem+44, mem+52, mem+60 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		any( a3 != v16int(3,11,19,27,35,43,51,59) ) ||
		any( a4 != v16int(4,12,20,28,36,44,52,60) ) ||
		any( a5 != v16int(5,13,21,29,37,45,53,61) ) ||
		any( a6 != v16int(6,14,22,30,38,46,54,62) ) ||
		any( a7 != v16int(7,15,23,31,39,47,55,63) ) ||
		i != 64 );
}

TEST(v16, test_store_8x4_tr_a)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0( 0, 4, 8,12,16,20,24,28);
  v16int a1( 1, 5, 9,13,17,21,25,29);
  v16int a2( 2, 6,10,14,18,22,26,30);
  v16int a3( 3, 7,11,15,19,23,27,31);
  v16int a4(32,36,40,44,48,52,56,60);
  v16int a5(33,37,41,45,49,53,57,61);
  v16int a6(34,38,42,46,50,54,58,62);
  v16int a7(35,39,43,47,51,55,59,63);

  int i;

  for( i=0; i < 64; i++ ) mem[i] = 0;

  store_8x4_tr( a0, a1, a2, a3,
		mem,    mem+ 4, mem+ 8, mem+12, mem+16, mem+20, mem+24, mem+28 );

  store_8x4_tr( a4, a5, a6, a7,
		mem+32, mem+36, mem+40, mem+44, mem+48, mem+52, mem+56, mem+60 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int( 0, 4, 8,12,16,20,24,28) ) ||
		any( a1 != v16int( 1, 5, 9,13,17,21,25,29) ) ||
		any( a2 != v16int( 2, 6,10,14,18,22,26,30) ) ||
		any( a3 != v16int( 3, 7,11,15,19,23,27,31) ) ||
		any( a4 != v16int(32,36,40,44,48,52,56,60) ) ||
		any( a5 != v16int(33,37,41,45,49,53,57,61) ) ||
		any( a6 != v16int(34,38,42,46,50,54,58,62) ) ||
		any( a7 != v16int(35,39,43,47,51,55,59,63) ) ||
		i != 64 );
}

TEST(v16, test_store_8x8_tr)
{
  DECLARE_ALIGNED_ARRAY( int, 64, mem, 64 );

  v16int a0(0, 8,16,24,32,40,48,56);
  v16int a1(1, 9,17,25,33,41,49,57);
  v16int a2(2,10,18,26,34,42,50,58);
  v16int a3(3,11,19,27,35,43,51,59);
  v16int a4(4,12,20,28,36,44,52,60);
  v16int a5(5,13,21,29,37,45,53,61);
  v16int a6(6,14,22,30,38,46,54,62);
  v16int a7(7,15,23,31,39,47,55,63);

  int i;

  for( i=0; i < 64; i++ ) mem[i] = 0;

  store_8x8_tr( a0, a1, a2, a3, a4, a5, a6, a7,
		mem, mem+ 8, mem+16, mem+24, mem+32, mem+40, mem+48, mem+56 );

  for( i=0; i < 64; i++ ) if( mem[i] != i ) break;

  ASSERT_FALSE( any( a0 != v16int(0, 8,16,24,32,40,48,56) ) ||
		any( a1 != v16int(1, 9,17,25,33,41,49,57) ) ||
		any( a2 != v16int(2,10,18,26,34,42,50,58) ) ||
		any( a3 != v16int(3,11,19,27,35,43,51,59) ) ||
		any( a4 != v16int(4,12,20,28,36,44,52,60) ) ||
		any( a5 != v16int(5,13,21,29,37,45,53,61) ) ||
		any( a6 != v16int(6,14,22,30,38,46,54,62) ) ||
		any( a7 != v16int(7,15,23,31,39,47,55,63) ) ||
		i != 64 );
}

#endif
