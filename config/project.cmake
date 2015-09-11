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
