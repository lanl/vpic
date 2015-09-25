#~----------------------------------------------------------------------------~#
# VPIC project configuration
#~----------------------------------------------------------------------------~#

project(vpic)

#------------------------------------------------------------------------------#
# Set project defaults
#------------------------------------------------------------------------------#

set(ENABLE_MPI True CACHE BOOL "Enable MPI" FORCE)
set(ENABLE_MPI_CXX_BINDINGS False CACHE BOOL "Enable MPI C++ Bindings" FORCE)

#------------------------------------------------------------------------------#
# Add build options
#------------------------------------------------------------------------------#

option(ENABLE_INTEGRATED_TESTS "enable integrated tests" OFF)
if(ENABLE_INTEGRATED_TESTS)
  enable_testing()
  add_subdirectory(test/integrated)
endif(ENABLE_INTEGRATED_TESTS)

option(USE_OPENMP "use openmp" OFF)
if(USE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)

set(USE_V4)
option(USE_V4_ALTIVEC "Enable V4 Altivec" OFF)
if(USE_V4_ALTIVEC)
  add_definitions(-DUSE_V4_ALTIVEC)
  set(USE_V4 True)
endif(USE_V4_ALTIVEC)

option(USE_V4_PORTABLE "Enable V4 Portable" OFF)
if(USE_V4_PORTABLE)
  add_definitions(-DUSE_V4_PORTABLE)
  set(USE_V4 True)
endif(USE_V4_PORTABLE)

option(USE_V4_SSE "Enable V4 SSE" OFF)
if(USE_V4_SSE)
  add_definitions(-DUSE_V4_SSE)
  set(USE_V4 True)
endif(USE_V4_SSE)

#------------------------------------------------------------------------------#
# Add library target
#------------------------------------------------------------------------------#

cinch_add_library_target(vpic src)

#------------------------------------------------------------------------------#
# Set header suffix
#------------------------------------------------------------------------------#

set(CINCH_HEADER_SUFFIXES "\\.h")

#~---------------------------------------------------------------------------~-#
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
