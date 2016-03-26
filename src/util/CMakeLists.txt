#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

set(util_HEADERS
  bitfield.h
  checksum.h
  swap.h
  system.h
  util.h
  util_base.h
  checkpt/checkpt.h
  checkpt/checkpt_io.h
  checkpt/checkpt_private.h
  io/FileIO.h
  io/FileIOData.h
  io/FileUtils.h
  io/P2PIOPolicy.h
  io/P2PUtilsPolicy.h
  io/StandardIOPolicy.h
  io/StandardUtilsPolicy.h
  mp/DMPPolicy.h
  mp/mp.h
  mp/MPWrapper.h
  mp/RelayPolicy.h
  pipelines/pipelines.h
  profile/profile.h
  rng/drandn_table.h
  rng/frandn_table.h
  rng/rng.h
  rng/rng_private.h
  v4/v4.h
  v4/v4_altivec.h
  v4/v4_portable.h
  v4/v4_sse.h
  v8/v8.h
  v8/v8_avx.h
  v8/v8_avx2.h
  v8/v8_portable.h
  PARENT_SCOPE
)

set(util_SOURCES
  checkpt/checkpt.c
  checkpt/checkpt_io.cc
  mp/mp.cc
  boot.c
  util_base.c
  pipelines/pipelines_serial.c
  pipelines/pipelines_thread.c
  profile/profile.c
  rng/drandn_table.c
  rng/frandn_table.c
  rng/rng.c
  rng/rng_pool.c
  PARENT_SCOPE
)

if(ENABLE_UNIT_TESTS)

  if(USE_V4)
    cinch_add_unit(v4
      SOURCES v4/test/v4.cc)
  endif(USE_V4)

  if(USE_V8)
    cinch_add_unit(v8
      SOURCES v8/test/v8.cc)
  endif(USE_V8)
  
  cinch_add_unit(rng
    SOURCES rng/test/rng.cc
    LIBRARIES vpic dl
    POLICY VPIC
    THREADS 1)

#  cinch_add_unit(casename
#    SOURCES testfile.cc
#    LIBRARIES list
#    INCLUDES list
#    POLICY MPI
#    THREADS 1 2 4
#  )

endif(ENABLE_UNIT_TESTS)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
