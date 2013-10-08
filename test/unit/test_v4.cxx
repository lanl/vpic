#include <util.h>

#ifndef V4_ACCELERATION

int
main( int argc, char ** argv ) {
  boot_services( &argc, &argv );
  MESSAGE(( "No V4_ACCELERATION specified" ));
  halt_services();
  return 0;
}

#else

using namespace v4;

#define TEST_ASSIGN(type,name,op,ai,bi,af,bf)                           \
  void test_##type##_assign_##name(void) {                              \
    v4##type a ai, b bi;                                                \
    a op b;                                                             \
    if( any(a!=v4##type af) || any(b!=v4##type bf) )                    \
      ERROR(( ""#type"_assign_"#name": FAIL" ));                        \
    MESSAGE(( ""#type"_assign_"#name": pass" ));			\
  }

#define TEST_PREFIX_UNARY(type,name,op,ai,bi,af,bf)                     \
  void test_##type##_prefix_unary_##name( void ) {                      \
    v4##type a ai, b bi;                                                \
    b = op a;                                                           \
    if( any(a!=v4##type af) || any(b!=v4##type bf) )                    \
      ERROR(( ""#type"_prefix_unary_"#name": FAIL" ));                  \
    MESSAGE(( ""#type"_prefix_unary_"#name": pass" ));			\
  }

#define TEST_POSTFIX_UNARY(type,name,op,ai,bi,af,bf)                    \
  void test_##type##_postfix_unary_##name( void ) {                     \
    v4##type a ai, b bi;                                                \
    b = a op;                                                           \
    if( any(a!=v4##type af) || any(b!=v4##type bf) )                    \
      ERROR(( ""#type"_postfix_unary_"#name": FAIL" ));			\
    MESSAGE(( ""#type"_postfix_unary_"#name": pass" ));			\
  }

#define TEST_BINARY(type,name,op,ai,bi,ci,cf)                           \
  void test_##type##_binary_##name(void) {                              \
    v4##type a ai, b bi, c ci;                                          \
    c = a op b;                                                         \
    if( any(a!=v4##type ai) || any(b!=v4##type bi) || any(c!=v4##type cf) ) \
      ERROR(( ""#type"_binary_"#name": FAIL" ));                        \
    MESSAGE(( ""#type"_binary_"#name": pass" ));			\
  }

#define TEST_LOGICAL(type,name,op,ai,bi,ci,cf)                          \
  void test_##type##_binary_##name(void) {                              \
    v4##type a ai, b bi;                                                \
    v4int c ci;                                                         \
    c = a op b;                                                         \
    if( any(a!=v4##type ai) || any(b!=v4##type bi) || any(c!=v4int cf) ) \
      ERROR(( ""#type"_logical_"#name": FAIL" ));			\
    MESSAGE(( ""#type"_logical_"#name": pass" ));			\
  }

///////////////////////////////////////////////////////////////////////////////
// class v4 tests

void test_any(void) {
  v4int a;
  int i;
  for( i=0; i<16; i++ ) {
    a[0] = i&1, a[1] = i&2, a[2] = i&4, a[3] = i&8;
    if( ( i>0 && !any(a) ) || ( i==0 && any(a) ) ) ERROR(( "any: FAIL" ));
  }
  MESSAGE(( "any: pass" ));
}

void test_all(void) {
  v4int a;
  int i;
  for( i=0; i<16; i++ ) {
    a[0] = i&1, a[1] = i&2, a[2] = i&4, a[3] = i&8;
    if( ( i<15 && all(a) ) || ( i==15 && !all(a) ) ) ERROR(( "all: FAIL" ));
  }
  MESSAGE(( "all: pass" ));
}

void test_splat(void) {
  v4int a( 1, 2, 3, 4), b( 5, 6, 7, 8), c( 9,10,11,12), d(13,14,15,16);
  v4int e(17,18,19,20);
  b = splat<0>(a); c = splat<1>(a); d = splat<2>(a); e = splat<3>(a);
  if( any(a!=v4int(1,2,3,4)) || any(b!=v4int(1,1,1,1)) ||
      any(c!=v4int(2,2,2,2)) || any(d!=v4int(3,3,3,3)) ||
      any(e!=v4int(4,4,4,4)) )
    ERROR(( "splat: FAIL" ));
  MESSAGE(( "splat: pass" ));
}

void test_shuffle(void) {
  v4int a( 0, 1, 2, 3), b( 4, 8,12,16), c( 5, 9,13,17), d( 6,10,14,18);
  v4int e( 7,11,15,19);
  b = shuffle<1,2,3,0>(a);
  c = shuffle<2,3,0,1>(a);
  d = shuffle<3,0,1,2>(a);
  e = shuffle<3,2,1,0>(a);
  if( any(a!=v4int(0,1,2,3)) || any(b!=v4int(1,2,3,0)) ||
      any(c!=v4int(2,3,0,1)) || any(d!=v4int(3,0,1,2)) ||
      any(e!=v4int(3,2,1,0)) )
    ERROR(( "shuffle: FAIL" ));
  MESSAGE(( "shuffle: pass" ));
}

void test_swap(void) {
  v4int a(1,2,3,4), b(5,6,7,8);
  swap(a,b);
  if( any(a!=v4int(5,6,7,8)) || any(b!=v4int(1,2,3,4)) )
    ERROR(( "swap: FAIL" ));
  MESSAGE(( "swap: pass" ));
}

void test_transpose(void) {
  v4int a0( 0, 1, 2, 3), a1( 4, 5, 6, 7), a2( 8, 9,10,11), a3(12,13,14,15);
  transpose(a0,a1,a2,a3);
  if( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
      any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) )
    ERROR(( "transpose: FAIL" ));
  MESSAGE(( "transpose: pass" ));
}

void test_load_4x1(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0(1,0,0,0), a1(0,0,0,0), a2(0,0,0,0), a3(0,0,0,0);
  int i;
  for( i=0; i<16; i++ ) mem[i] = i;
  load_4x1( mem,    a0 );
  load_4x1( mem+4,  a1 );
  load_4x1( mem+8,  a2 );
  load_4x1( mem+12, a3 );
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 1, 2, 3)) || any(a1!=v4int( 4, 5, 6, 7)) ||
      any(a2!=v4int( 8, 9,10,11)) || any(a3!=v4int(12,13,14,15)) || i!=16 )
    ERROR(( "load_4x1: FAIL" ));
  MESSAGE(( "load_4x1: pass" ));
}

void test_store_4x1(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 1, 2, 3), a1( 4, 5, 6, 7), a2( 8, 9,10,11), a3(12,13,14,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x1(a0,mem);
  store_4x1(a1,mem+4);
  store_4x1(a2,mem+8);
  store_4x1(a3,mem+12);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 1, 2, 3)) || any(a1!=v4int( 4, 5, 6, 7)) ||
      any(a2!=v4int( 8, 9,10,11)) || any(a3!=v4int(12,13,14,15)) || i!=16 )
    ERROR(( "store_4x1: FAIL" ));
  MESSAGE(( "store_4x1: pass" ));
}

void test_stream_4x1(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 1, 2, 3), a1( 4, 5, 6, 7), a2( 8, 9,10,11), a3(12,13,14,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  stream_4x1(a0,mem);
  stream_4x1(a1,mem+4);
  stream_4x1(a2,mem+8);
  stream_4x1(a3,mem+12);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 1, 2, 3)) || any(a1!=v4int( 4, 5, 6, 7)) ||
      any(a2!=v4int( 8, 9,10,11)) || any(a3!=v4int(12,13,14,15)) || i!=16 )
    ERROR(( "stream_4x1: FAIL" ));
  MESSAGE(( "stream_4x1: pass" ));
}

void test_clear_4x1(void) {
  v4float vmem[4]; float * mem = (float *)vmem;
  int i;
  for(i=0; i<16; i++) mem[i] = i;
  clear_4x1( mem+8  );
  clear_4x1( mem+12 );
  for(i=8; i<16; i++) mem[i] += i;
  for(i=0; i<16; i++) if( mem[i]!=i ) break;
  if( i!=16 )
    ERROR(( "clear_4x1: FAIL" ));
  MESSAGE(( "clear_4x1: pass" ));
}

void test_copy_4x1(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  int i;
  for( i=0; i<8; i++ ) mem[i] = i;
  copy_4x1(mem+8, mem  );
  copy_4x1(mem+12,mem+4);
  for( i=8; i<16; i++ ) mem[i] += 8;
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( i!=16 )
    ERROR(( "copy_4x1: FAIL" ));
  MESSAGE(( "copy_4x1: pass" ));
}

void test_swap_4x1(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  int i;
  for( i=0; i<8; i++ ) mem[i] = i;
  copy_4x1(mem+12, mem  );
  copy_4x1(mem+8,  mem+4);
  for( i=8; i<16; i++ ) mem[i] += 8;
  swap_4x1(mem+8,  mem+12 );
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( i!=16 )
    ERROR(( "swap_4x1: FAIL" ));
  MESSAGE(( "swap_4x1: pass" ));
}

void test_load_4x1_tr(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0, a1, a2, a3;
  int i;
  for( i=0; i<16; i++ ) mem[i] = i;
  load_4x1_tr(mem,  mem+4,mem+8, mem+12,a0);
  load_4x1_tr(mem+1,mem+5,mem+9, mem+13,a1);
  load_4x1_tr(mem+2,mem+6,mem+10,mem+14,a2);
  load_4x1_tr(mem+3,mem+7,mem+11,mem+15,a3);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
      any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) ||
      i!=16 )
    ERROR(( "load_4x1_tr: FAIL" ));
  MESSAGE(( "load_4x1_tr: pass" ));
}

void test_load_4x2_tr(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0, a1, a2, a3;
  int i;
  for( i=0; i<16; i++ ) mem[i] = i;
  load_4x2_tr(mem,  mem+4,mem+8, mem+12,a0,a1);
  load_4x2_tr(mem+2,mem+6,mem+10,mem+14,a2,a3);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
      any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 )
    ERROR(( "load_4x2_tr: FAIL" ));
  MESSAGE(( "load_4x2_tr: pass" ));
}

void test_load_4x3_tr(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0, a1, a2;
  int i;
  for( i=0; i<16; i++ ) mem[i] = i;
  load_4x3_tr(mem,mem+4,mem+8,mem+12,a0,a1,a2);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
      any(a2!=v4int( 2, 6,10,14)) || i!=16 )
    ERROR(( "load_4x3_tr: FAIL" ));
  MESSAGE(( "load_4x3_tr: pass" ));
}

void test_load_4x4_tr(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0, a1, a2, a3;
  int i;
  for( i=0; i<16; i++ ) mem[i] = i;
  load_4x4_tr(mem,mem+4,mem+8,mem+12,a0,a1,a2,a3);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
      any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 )
    ERROR(( "load_4x4_tr: FAIL" ));
  MESSAGE(( "load_4x4_tr: pass" ));
}

void test_store_4x1_tr(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 4, 8,12), a1( 1, 5, 9,13), a2( 2, 6,10,14), a3( 3, 7,11,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x1_tr(a0,mem,  mem+4,mem+8, mem+12);
  store_4x1_tr(a1,mem+1,mem+5,mem+9, mem+13);
  store_4x1_tr(a2,mem+2,mem+6,mem+10,mem+14);
  store_4x1_tr(a3,mem+3,mem+7,mem+11,mem+15);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
      any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 )
    ERROR(( "store_4x1_tr: FAIL" ));
  MESSAGE(( "store_4x1_tr: pass" ));
}

void test_store_4x2_tr(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 4, 8,12), a1( 1, 5, 9,13), a2( 2, 6,10,14), a3( 3, 7,11,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x2_tr(a0,a1,mem,  mem+4,mem+8, mem+12);
  store_4x2_tr(a2,a3,mem+2,mem+6,mem+10,mem+14);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
      any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 )
    ERROR(( "store_4x2_tr: FAIL" ));
  MESSAGE(( "store_4x2_tr: pass" ));
}

void test_store_4x3_tr(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 4, 8,12), a1( 1, 5, 9,13), a2( 2, 6,10,14), a3( 3, 7,11,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x3_tr(a0,a1,a2,mem,  mem+4,mem+8, mem+12);
  for( i=0; i<16; i++ )
    if( ( (i&3)!=3 && mem[i]!=i ) || ( (i&3)==3 && mem[i]!=0 ) ) break;
  if( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
      any(a2!=v4int( 2, 6,10,14)) || i!=16 )
    ERROR(( "store_4x3_tr: FAIL" ));
  MESSAGE(( "store_4x3_tr: pass" ));
}

void test_store_4x4_tr(void) {
  DECLARE_ALIGNED_ARRAY( int, 16, mem, 16 );
  v4int a0( 0, 4, 8,12), a1( 1, 5, 9,13), a2( 2, 6,10,14), a3( 3, 7,11,15);
  int i;
  for( i=0; i<16; i++ ) mem[i] = 0;
  store_4x4_tr(a0,a1,a2,a3,mem,  mem+4,mem+8, mem+12);
  for( i=0; i<16; i++ ) if( mem[i]!=i ) break;
  if( any(a0!=v4int( 0, 4, 8,12)) || any(a1!=v4int( 1, 5, 9,13)) ||
      any(a2!=v4int( 2, 6,10,14)) || any(a3!=v4int( 3, 7,11,15)) || i!=16 )
    ERROR(( "store_4x4_tr: FAIL" ));
  MESSAGE(( "store_4x4_tr: pass" ));
}

///////////////////////////////////////////////////////////////////////////////
// class v4int tests

void test_int_scalar_constructor(void) {
  v4int a0(0), a1(1), a2(2), a3(3);
  if( any(a0!=v4int( 0, 0, 0, 0)) || any(a1!=v4int( 1, 1, 1, 1)) ||
      any(a2!=v4int( 2, 2, 2, 2)) || any(a3!=v4int( 3, 3, 3, 3)) )
    ERROR(( "int_scalar_constructor: FAIL" ));
  MESSAGE(( "int_scalar_constructor: pass" ));
}

TEST_ASSIGN(int,eq,    =, (0,1,2,3), (4,5,6,7),(4,5,6,7),    (4,5,6,7))
TEST_ASSIGN(int,addeq,+=, (0,1,2,3), (4,5,6,7),(4,6,8,10),   (4,5,6,7))
TEST_ASSIGN(int,subeq,-=, (0,1,2,3), (4,5,6,7),(-4,-4,-4,-4),(4,5,6,7))
TEST_ASSIGN(int,muleq,*=, (0,1,2,3), (4,5,6,7),(0,5,12,21),  (4,5,6,7))
TEST_ASSIGN(int,diveq,/=, (5,6,7,8), (1,2,3,4),(5,3,2,2),    (1,2,3,4))
TEST_ASSIGN(int,remeq,%=, (5,6,7,8), (1,2,3,4),(0,0,1,0),    (1,2,3,4))
TEST_ASSIGN(int,xoreq,^=, (0,0,1,1), (0,1,0,1),(0,1,1,0),    (0,1,0,1))
TEST_ASSIGN(int,andeq,&=, (0,0,1,1), (0,1,0,1),(0,0,0,1),    (0,1,0,1))
TEST_ASSIGN(int,oreq, |=, (0,0,1,1), (0,1,0,1),(0,1,1,1),    (0,1,0,1))
TEST_ASSIGN(int,lsheq,<<=,(1,1,1,1), (1,2,3,4),(2,4,8,16),   (1,2,3,4))
TEST_ASSIGN(int,rsheq,>>=,(2,4,8,16),(1,2,3,4),(1,1,1,1),    (1,2,3,4))

// FIXME: This test should also use the access operator to modify a v4int
void test_int_access( void ) {
  v4int a(0,1,2,3);
  int i;
  for( i=0; i<4; i++ ) if( a(i)!=i ) break;
  for( i=0; i<4; i++ ) if( a[i]!=i ) break;
  if( i!=4 )
    ERROR(( "int_access: FAIL" ));
  MESSAGE(( "int_access: pass" ));
}

TEST_PREFIX_UNARY(int,plus,+,(0,1,2,3),(1,0,0,0),(0,1,2,3),(0,1,2,3))
TEST_PREFIX_UNARY(int,neg,-,(0,1,2,3),(1,0,0,0),(0,1,2,3),(0,-1,-2,-3))
TEST_PREFIX_UNARY(int,bnot,~,(0,1,2,3),(0,0,0,0),(0,1,2,3),(-1,-2,-3,-4))
TEST_PREFIX_UNARY(int,lnot,!,(0,1,2,3),(1,2,3,0),(0,1,2,3),(-1,0,0,0))

TEST_PREFIX_UNARY(int,inc,++,(0,1,2,3),(1,0,0,0),(1,2,3,4),(1,2,3,4))
TEST_PREFIX_UNARY(int,dec,--,(1,2,3,4),(0,0,0,0),(0,1,2,3),(0,1,2,3))

TEST_POSTFIX_UNARY(int,inc,++,(0,1,2,3),(0,0,0,0),(1,2,3,4),(0,1,2,3))
TEST_POSTFIX_UNARY(int,dec,--,(1,2,3,4),(1,0,0,0),(0,1,2,3),(1,2,3,4))

TEST_BINARY(int,add,+,(1,2,3,4),(4,3,2,1),(0,0,0,0),(5,5,5,5))
TEST_BINARY(int,sub,-,(1,2,3,4),(4,3,2,1),(0,0,0,0),(-3,-1,1,3))
TEST_BINARY(int,mul,*,(1,2,3,4),(4,3,2,1),(1,1,1,1),(4,6,6,4))
TEST_BINARY(int,div,/,(1,2,3,4),(4,3,2,1),(5,5,5,5),(0,0,1,4))
TEST_BINARY(int,rem,%,(1,2,3,4),(4,3,2,1),(3,3,3,3),(1,2,1,0))
TEST_BINARY(int,bxor,^,(0,0,1,1),(0,1,0,1),(2,2,2,2),(0,1,1,0))
TEST_BINARY(int,band,&,(0,0,1,1),(0,1,0,1),(2,2,2,2),(0,0,0,1))
TEST_BINARY(int,bor, |,(0,0,1,1),(0,1,0,1),(2,2,2,2),(0,1,1,1))
TEST_BINARY(int,lsh,<<,( 1,1,1,1),(4,3,2,1),(0,0,0,0),(16,8,4,2))
TEST_BINARY(int,rsh,>>,(16,8,4,2),(4,3,2,1),(0,0,0,0),(1,1,1,1))

TEST_LOGICAL(int,lt,<,  (1,2,3,4),(3,2,2,1),(3,0,6,-7),(-1,0,0,0))
TEST_LOGICAL(int,gt,>,  (1,2,3,4),(3,2,2,1),(3,0,6,-7),(0,0,-1,-1))
TEST_LOGICAL(int,eq,==, (1,2,3,4),(3,2,2,1),(3,0,6,-7),(0,-1,0,0))
TEST_LOGICAL(int,ne,!=, (1,2,3,4),(3,2,2,1),(3,0,6,-7),(-1,0,-1,-1))
TEST_LOGICAL(int,le,<=, (1,2,3,4),(3,2,2,1),(3,0,6,-7),(-1,-1,0,0))
TEST_LOGICAL(int,ge,>=, (1,2,3,4),(3,2,2,1),(3,0,6,-7),(0,-1,-1,-1))
TEST_LOGICAL(int,land,&&,(0,0,3,4),(0,2,0,1),(3,0,6,-7),(0,0,0,-1))
TEST_LOGICAL(int,lor,||, (0,0,3,4),(0,2,0,1),(3,0,6,-7),(0,-1,-1,-1))

void test_abs(void) {
  v4int a(0,-1,2,-3), b(4,4,4,4);
  b = abs(a);
  if( any(a!=v4int(0,-1,2,-3)) || any(b!=v4int(0, 1,2, 3)) )
    ERROR(( "abs: FAIL" ));
  MESSAGE(( "abs: pass" ));
}

void test_czero(void) {
  v4int a( 1, 2, 3, 4), b( 5, 6, 7, 8), c(-1, 0,-1, 0);
  b=czero(c,a);
  if( any(a!=v4int( 1,2, 3,4)) || any(b!=v4int( 0,2, 0,4)) ||
      any(c!=v4int(-1,0,-1,0)) )
    ERROR(( "czero: FAIL" ));
  MESSAGE(( "czero: pass" ));
}

void test_notczero(void) {
  v4int a( 1, 2, 3, 4), b( 5, 6, 7, 8), c(-1, 0,-1, 0);
  b=notczero(c,a);
  if( any(a!=v4int( 1,2, 3,4)) || any(b!=v4int( 1,0, 3,0)) ||
      any(c!=v4int(-1,0,-1,0)) )
    ERROR(( "notczero: FAIL" ));
  MESSAGE(( "notczero: pass" ));
}

void test_merge(void) {
  v4int a( 0, 1, 2, 3), b( 4, 5, 6, 7), c(-1, 0,-1, 0), d( 8, 9,10,11);
  d = merge(c,b,a);
  if( any(a!=v4int(0,1,2,3)  ) || any(b!=v4int(4,5,6,7)  ) ||
      any(c!=v4int(-1,0,-1,0)) || any(d!=v4int(4,1,6,3)  ) )
    ERROR(( "merge: FAIL" ));
  MESSAGE(( "merge: pass" ));
}

///////////////////////////////////////////////////////////////////////////////
// v4float class

void test_float_scalar_constructor(void) {
  v4float a0(0), a1(1), a2(2), a3(3);
  if( any(a0!=v4float( 0, 0, 0, 0)) || any(a1!=v4float( 1, 1, 1, 1)) ||
      any(a2!=v4float( 2, 2, 2, 2)) || any(a3!=v4float( 3, 3, 3, 3)) )
    ERROR(( "float_scalar_constructor: FAIL" ));
  MESSAGE(( "float_scalar_constructor: pass" ));
}

TEST_ASSIGN(float,eq,    =, (0,1,2,3), (4,5,6,7),(4,5,6,7),    (4,5,6,7))
TEST_ASSIGN(float,addeq,+=, (0,1,2,3), (4,5,6,7),(4,6,8,10),   (4,5,6,7))
TEST_ASSIGN(float,subeq,-=, (0,1,2,3), (4,5,6,7),(-4,-4,-4,-4),(4,5,6,7))
TEST_ASSIGN(float,muleq,*=, (0,1,2,3), (4,5,6,7),(0,5,12,21),  (4,5,6,7))
TEST_ASSIGN(float,diveq,/=, (5,6,7,8), (1,2,4,4),(5,3,1.75,2), (1,2,4,4))

// FIXME: THIS SHOULD TEST MODIFYING A V4FLOAT TOO
void test_float_access( void ) {
  v4float a(0,1,2,3);
  int i;
  for( i=0; i<4; i++ ) if( a(i)!=i ) break;
  for( i=0; i<4; i++ ) if( a[i]!=i ) break;
  if( i!=4 )
    ERROR(( "float_access: FAIL" ));
  MESSAGE(( "float_access: pass" ));
}

TEST_PREFIX_UNARY(float,plus,+,(0,1,2,3),(1,0,0,0),(0,1,2,3),(0,1,2,3))
TEST_PREFIX_UNARY(float,neg,-,(1,2,3,4),(1,0,0,0),(1,2,3,4),(-1,-2,-3,-4))

void test_float_prefix_unary_lnot( void ) {
  v4float a(0,1,2,3); v4int b(1,2,3,0);
  b = !a;
  if( any(a!=v4float(0,1,2,3)) || any(b!=v4int(-1,0,0,0)) )
    ERROR(( "float_prefix_unary_lnot: FAIL" ));
  MESSAGE(( "float_prefix_unary_lnot: pass" ));
}

TEST_PREFIX_UNARY(float,inc,++,(0,1,2,3),(1,0,0,0),(1,2,3,4),(1,2,3,4))
TEST_PREFIX_UNARY(float,dec,--,(1,2,3,4),(0,0,0,0),(0,1,2,3),(0,1,2,3))

TEST_POSTFIX_UNARY(float,inc,++,(0,1,2,3),(1,0,0,0),(1,2,3,4),(0,1,2,3))
TEST_POSTFIX_UNARY(float,dec,--,(1,2,3,4),(0,0,0,0),(0,1,2,3),(1,2,3,4))

TEST_BINARY(float,add,+,(1,2,3,4),(4,3,2,1),(0,0,0,0),(5,5,5,5))
TEST_BINARY(float,sub,-,(1,2,3,4),(4,3,2,1),(0,0,0,0),(-3,-1,1,3))
TEST_BINARY(float,mul,*,(1,2,3,4),(4,3,2,1),(1,1,1,1),(4,6,6,4))
TEST_BINARY(float,div,/,(1,2,3,4),(4,2,2,1),(5,5,5,5),(0.25,1.,1.5,4))

TEST_LOGICAL(float,lt,<,  (1,2,3,4),(3,2,2,1),(3,0,6,-7),(-1,0,0,0))
TEST_LOGICAL(float,gt,>,  (1,2,3,4),(3,2,2,1),(3,0,6,-7),(0,0,-1,-1))
TEST_LOGICAL(float,eq,==, (1,2,3,4),(3,2,2,1),(3,0,6,-7),(0,-1,0,0))
TEST_LOGICAL(float,ne,!=, (1,2,3,4),(3,2,2,1),(3,0,6,-7),(-1,0,-1,-1))
TEST_LOGICAL(float,le,<=, (1,2,3,4),(3,2,2,1),(3,0,6,-7),(-1,-1,0,0))
TEST_LOGICAL(float,ge,>=, (1,2,3,4),(3,2,2,1),(3,0,6,-7),(0,-1,-1,-1))
TEST_LOGICAL(float,land,&&,(0,0,3,4),(0,2,0,1),(3,0,6,-7),(0,0,0,-1))
TEST_LOGICAL(float,lor,||, (0,0,3,4),(0,2,0,1),(3,0,6,-7),(0,-1,-1,-1))

// FIXME: MATH LIB TESTS GO HERE!

void test_float_fabs(void) {
  v4float a(1,-2,3,-4), b(5,6,7,8);
  b = fabs(a);
  if( any(a!=v4float(1,-2,3,-4)) || any(b!=v4float(1, 2,3, 4)) )
    ERROR(( "float_fabs: FAIL" ));
  MESSAGE(( "float_fabs: pass" ));
}

void test_float_sqrt(void) {
  v4float a(0.01,0.5,4,10000), b(0.1,0.707106781,2,100);
  b = fabs( ( b - sqrt(a) ) / b );
  if( any(a!=v4float(0.01,0.5,4,10000)) || any(b>v4float(1e-6)) )
    ERROR(( "float_sqrt: FAIL" ));
  MESSAGE(( "float_sqrt: pass" ));
}

void test_float_copysign(void) {
  v4float a(1,-2, 3,-4), b(5,6,-7,-8), c(9,10,11,12);
  c = copysign(a,b);
  if( any(a!=v4float(1,-2, 3,-4)) || any(b!=v4float(5, 6,-7,-8)) ||
      any(c!=v4float(1, 2,-3,-4)) )
    ERROR(( "float_copysign: FAIL" ));
  MESSAGE(( "float_copysign: pass" ));
}

void test_float_rsqrt_approx(void) {
  v4float a(0.01,0.5,4,10000), b(10,1.41421356,0.5,0.01);
  b = fabs( ( b - rsqrt_approx(a) ) / b );
  if( any(a!=v4float(0.01,0.5,4,10000)) || any(b>v4float(1e-3)) )
    ERROR(( "float_rsqrt_approx: FAIL" ));
  MESSAGE(( "float_rsqrt_approx: pass" ));
}

void test_float_rsqrt(void) {
  v4float a(0.01,0.5,4,10000), b(10,1.41421356,0.5,0.01);
  b = fabs( ( b - rsqrt(a) ) / b );
  if( any(a!=v4float(0.01,0.5,4,10000)) || any(b>v4float(1e-6)) )
    ERROR(( "float_rsqrt: FAIL" ));
  MESSAGE(( "float_rsqrt: pass" ));
}

void test_float_rcp_approx(void) {
  v4float a(0.01,0.5,4,10000), b(100,2,0.25,0.0001);
  b = fabs( ( b - rcp_approx(a) ) / b );
  if( any(a!=v4float(0.01,0.5,4,10000)) || any(b>v4float(1e-3)) )
    ERROR(( "float_rcp_approx: FAIL" ));
  MESSAGE(( "float_rcp_approx: pass" ));
}

void test_float_rcp(void) {
  v4float a(0.01,0.5,4,10000), b(100,2,0.25,0.0001);
  b = fabs( ( b - rcp(a) ) / b );
  if( any(a!=v4float(0.01,0.5,4,10000)) || any(b>v4float(1e-6)) )
    ERROR(( "float_rcp: FAIL" ));
  MESSAGE(( "float_rcp: pass" ));
}

void test_float_fma(void) {
  v4float a(1,2,3,4), b(5,6,7,8), c(9,10,11,12), d(13,14,15,16);
  d = fma(a,b,c);
  if( any(a!=v4float( 1, 2, 3, 4)) || any(b!=v4float( 5, 6, 7, 8)) ||
      any(c!=v4float( 9,10,11,12)) || any(d!=v4float(14,22,32,44)) )
    ERROR(( "float_fma: FAIL" ));
  MESSAGE(( "float_fma: pass" ));
}

void test_float_fms(void) {
  v4float a(1,2,3,4), b(5,6,7,8), c(9,10,11,12), d(13,14,15,16);
  d = fms(a,b,c);
  if( any(a!=v4float( 1, 2, 3, 4)) || any(b!=v4float( 5, 6, 7, 8)) ||
      any(c!=v4float( 9,10,11,12)) || any(d!=v4float(-4, 2,10,20)) )
    ERROR(( "float_fms: FAIL" ));
  MESSAGE(( "float_fms: pass" ));
}

void test_float_fnms(void) {
  v4float a(1,2,3,4), b(5,6,7,8), c(9,10,11,12), d(13,14,15,16);
  d = fnms(a,b,c);
  if( any(a!=v4float( 1, 2,  3,  4)) || any(b!=v4float( 5, 6,  7,  8)) ||
      any(c!=v4float( 9,10, 11, 12)) || any(d!=v4float( 4,-2,-10,-20)) )
    ERROR(( "float_fnms: FAIL" ));
  MESSAGE(( "float_fnms: pass" ));
}

void test_float_clear_bits(void) {
  v4int   a(1<<31,1<<31,-1,0); v4float b(-1,-2,-3,-4), c(5,6,7,8);
  c = clear_bits(a,b);
  if( any(a!=v4int(1<<31,1<<31,-1,0)) || any(b!=v4float(-1,-2,-3,-4)) ||
      any(c!=v4float( 1, 2, 0,-4)) )
    ERROR(( "float_clear_bits: FAIL" ));
  MESSAGE(( "float_clear_bits: pass" ));
}

void test_float_set_bits(void) {
  v4int   a(1<<31,1<<31,0,0); v4float b(1,2,3,4), c(5,6,7,8);
  c = set_bits(a,b);
  if( any(a!=v4int(1<<31,1<<31,0,0)) || any(b!=v4float( 1, 2, 3, 4))   ||
      any(c!=v4float(-1,-2, 3, 4)) )
    ERROR(( "float_set_bits: FAIL" ));
  MESSAGE(( "float_set_bits: pass" ));
}

void test_float_toggle_bits(void) {
  v4int   a(1<<31,1<<31,0,0); v4float b(1,-2,3,-4), c(5,6,7,8);
  c = toggle_bits(a,b);
  if( any(a!=v4int(1<<31,1<<31,0,0)) || any(b!=v4float( 1,-2, 3,-4))   ||
      any(c!=v4float(-1, 2, 3,-4)) )
    ERROR(( "float_toggle_bits: FAIL" ));
  MESSAGE(( "float_toggle_bits: pass" ));
}

void test_float_increment_4x1(void) {
  DECLARE_ALIGNED_ARRAY( float, 16, mem, 4 );
  v4float a(1,2,3,4);
  mem[0] = 5; mem[1] = 6; mem[2] = 7; mem[3] = 8;
  increment_4x1( mem, a );
  if( any(a!=v4float(1,2,3,4)) ||
      mem[0]!=6  || mem[1]!=8 || mem[2]!=10 || mem[3]!=12 )
    ERROR(( "float_increment_4x1: FAIL" ));
  MESSAGE(( "float_increment_4x1: pass" ));
}

void test_float_decrement_4x1(void) {
  DECLARE_ALIGNED_ARRAY( float, 16, mem, 4 );
  v4float a(1,2,3,4);
  mem[0] = 5; mem[1] = 6; mem[2] = 7; mem[3] = 8;
  decrement_4x1( mem, a );
  if( any(a!=v4float(1,2,3,4)) ||
      mem[0]!=4 || mem[1]!=4 || mem[2]!=4 || mem[3]!=4 )
    ERROR(( "float_decrement_4x1: FAIL" ));
  MESSAGE(( "float_decrement_4x1: pass" ));
}

void test_float_scale_4x1(void) {
  DECLARE_ALIGNED_ARRAY( float, 16, mem, 4 );
  v4float a(1,2,3,4);
  mem[0] = 5; mem[1] = 6; mem[2] = 7; mem[3] = 8;
  scale_4x1( mem, a );
  if( any(a!=v4float(1,2,3,4)) ||
      mem[0]!=5  || mem[1]!=12 || mem[2]!=21 || mem[3]!=32 )
    ERROR(( "float_scale_4x1: FAIL" ));
  MESSAGE(( "float_scale_4x1: pass" ));
}

void test_float_trilinear(void) {
  v4float wl( -0.5, 0.25, -0.125, 0 ), wh(  1, 2, 3, 4 );
  trilinear(wl,wh);
  if( any(wl!=v4float(81./64., 27./64., 135./64., 45./64.)) ||
      any(wh!=v4float(63./64., 21./64., 105./64., 35./64.)) )
    ERROR(( "float_trilienar: FAIL" ));
  MESSAGE(( "float_trilinear: pass" ));
}

#undef TEST_ASSIGN
#undef TEST_PREFIX_UNARY
#undef TEST_POSTFIX_UNARY
#undef TEST_BINARY
#undef TEST_LOGICAL

int
main( int argc,
      char **argv ) {
  boot_services( &argc, &argv );

  test_any();
  test_all();
  test_splat();
  test_shuffle();
  test_swap();
  test_transpose();
  
  test_load_4x1();
  test_store_4x1();
  test_stream_4x1();
  test_clear_4x1();
  test_copy_4x1();
  test_swap_4x1();

  test_load_4x1_tr();
  test_load_4x2_tr();
  test_load_4x3_tr();
  test_load_4x4_tr();

  test_store_4x1_tr();
  test_store_4x2_tr();
  test_store_4x3_tr();
  test_store_4x4_tr();

  test_int_scalar_constructor();

  test_int_assign_eq();
  test_int_assign_addeq();
  test_int_assign_subeq();
  test_int_assign_muleq();
  test_int_assign_diveq();
  test_int_assign_remeq();
  test_int_assign_xoreq();
  test_int_assign_andeq();
  test_int_assign_oreq();
  test_int_assign_lsheq();
  test_int_assign_rsheq();

  test_int_access();

  test_int_prefix_unary_plus();
  test_int_prefix_unary_neg();
  test_int_prefix_unary_bnot();
  test_int_prefix_unary_lnot();

  test_int_prefix_unary_inc();
  test_int_prefix_unary_dec();

  test_int_postfix_unary_inc();
  test_int_postfix_unary_dec();

  test_int_binary_add();
  test_int_binary_sub();
  test_int_binary_mul();
  test_int_binary_div();
  test_int_binary_rem();
  test_int_binary_bxor();
  test_int_binary_band();
  test_int_binary_bor();
  test_int_binary_lsh();
  test_int_binary_rsh();

  test_int_binary_lt();
  test_int_binary_gt();
  test_int_binary_eq();
  test_int_binary_ne();
  test_int_binary_le();
  test_int_binary_ge();
  test_int_binary_land();
  test_int_binary_lor();

  test_abs();
  test_czero();
  test_notczero();
  test_merge();

  test_float_scalar_constructor();

  test_float_assign_eq();
  test_float_assign_addeq();
  test_float_assign_subeq();
  test_float_assign_muleq();
  test_float_assign_diveq();

  test_float_access();

  test_float_prefix_unary_plus();
  test_float_prefix_unary_neg();
  test_float_prefix_unary_lnot();

  test_float_prefix_unary_inc();
  test_float_prefix_unary_dec();

  test_float_postfix_unary_inc();
  test_float_postfix_unary_dec();

  test_float_binary_add();
  test_float_binary_sub();
  test_float_binary_mul();
  test_float_binary_div();

  test_float_binary_lt();
  test_float_binary_gt();
  test_float_binary_eq();
  test_float_binary_ne();
  test_float_binary_le();
  test_float_binary_ge();
  test_float_binary_land();
  test_float_binary_lor();

  test_float_fabs();
  test_float_sqrt();
  test_float_copysign();

  test_float_rsqrt_approx();
  test_float_rsqrt();

  test_float_rcp_approx();
  test_float_rcp();

  test_float_fma();
  test_float_fms();
  test_float_fnms();

  test_float_clear_bits();
  test_float_set_bits();
  test_float_toggle_bits();

  test_float_increment_4x1();
  test_float_decrement_4x1();
  test_float_scale_4x1();

  test_float_trilinear();

  halt_services();
  return 0;
}

#endif

