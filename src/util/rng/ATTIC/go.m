close all
clear all

fno = 1;
feps = 0.5/2^23;
deps = 0.5/2^52;

analyze(   'crand.bin',  'int8',  [ 0, 2^7-1  ], 1, fno ); fno = fno + 1;
analyze(   'hrand.bin',  'int16', [ 0, 2^15-1 ], 1, fno ); fno = fno + 1;
analyze(   'irand.bin',  'int32', [ 0, 2^31-1 ], 1, fno ); fno = fno + 1;
analyze(   'lrand.bin',  'int64', [ 0, 2^63-1 ], 1, fno ); fno = fno + 1;

analyze(  'i8rand.bin',  'int8',  [ 0, 2^7-1  ], 1, fno ); fno = fno + 1;
analyze( 'i16rand.bin',  'int16', [ 0, 2^15-1 ], 1, fno ); fno = fno + 1;
analyze( 'i32rand.bin',  'int32', [ 0, 2^31-1 ], 1, fno ); fno = fno + 1;
analyze( 'i64rand.bin',  'int64', [ 0, 2^63-1 ], 1, fno ); fno = fno + 1;

analyze(  'ucrand.bin', 'uint8',  [ 0, 2^8-1  ], 1, fno ); fno = fno + 1;
analyze(  'uhrand.bin', 'uint16', [ 0, 2^16-1 ], 1, fno ); fno = fno + 1;
analyze(  'uirand.bin', 'uint32', [ 0, 2^32-1 ], 1, fno ); fno = fno + 1;
analyze(  'ulrand.bin', 'uint64', [ 0, 2^64-1 ], 1, fno ); fno = fno + 1;

analyze(  'u8rand.bin', 'uint8',  [ 0, 2^8-1  ], 1, fno ); fno = fno + 1;
analyze( 'u16rand.bin', 'uint16', [ 0, 2^16-1 ], 1, fno ); fno = fno + 1;
analyze( 'u32rand.bin', 'uint32', [ 0, 2^32-1 ], 1, fno ); fno = fno + 1;
analyze( 'u64rand.bin', 'uint64', [ 0, 2^64-1 ], 1, fno ); fno = fno + 1;

analyze( 'frand.bin',    'float',  [ feps, 1-feps  ], 0, fno ); fno = fno + 1;
analyze( 'frand_c0.bin', 'float',  [ 0,    1-feps  ], 0, fno ); fno = fno + 1;
analyze( 'frand_c1.bin', 'float',  [ feps, 1       ], 0, fno ); fno = fno + 1;
analyze( 'frand_c.bin',  'float',  [ 0,    1       ], 0, fno ); fno = fno + 1;

analyze( 'drand.bin',    'double',  [ deps, 1-deps  ], 0, fno ); fno = fno + 1;
analyze( 'drand_c0.bin', 'double',  [ 0,    1-deps  ], 0, fno ); fno = fno + 1;
analyze( 'drand_c1.bin', 'double',  [ deps, 1       ], 0, fno ); fno = fno + 1;
analyze( 'drand_c.bin',  'double',  [ 0,    1       ], 0, fno ); fno = fno + 1;

analyze( 'frandn.bin', 'float',  [ -inf, inf ], 0, fno ); fno = fno + 1;
analyze( 'drandn.bin', 'double', [ -inf, inf ], 0, fno ); fno = fno + 1;

analyze( 'frande.bin', 'float',  [ 0, inf ], 0, fno ); fno = fno + 1;
analyze( 'drande.bin', 'double', [ 0, inf ], 0, fno ); fno = fno + 1;

