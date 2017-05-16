#~----------------------------------------------------------------------------~#
# VPIC project configuration
#~----------------------------------------------------------------------------~#

cinch_minimum_required(1.0)

project(vpic)

#------------------------------------------------------------------------------#
# Set project defaults
#------------------------------------------------------------------------------#

set(ENABLE_MPI True CACHE BOOL "Enable MPI" FORCE)
set(ENABLE_MPI_CXX_BINDINGS False CACHE BOOL "Enable MPI C++ Bindings" FORCE)

#------------------------------------------------------------------------------#
# Add build options
#------------------------------------------------------------------------------#

option(ENABLE_INTEGRATED_TESTS "Enable integrated tests" OFF)

option(USE_OPENMP "Use OpenMP" OFF)

option(USE_PTHREADS "Use Pthreads" ON)

option(USE_V4_ALTIVEC "Enable V4 Altivec" OFF)

option(USE_V4_PORTABLE "Enable V4 Portable" OFF)

option(USE_V4_SSE "Enable V4 SSE" OFF)

option(USE_V4_AVX "Enable V4 AVX" OFF)

option(USE_V4_AVX2 "Enable V4 AVX2" OFF)

option(USE_V8_PORTABLE "Enable V8 Portable" OFF)

option(USE_V8_AVX "Enable V8 AVX" OFF)

option(USE_V8_AVX2 "Enable V8 AVX2" OFF)

option(USE_V16_PORTABLE "Enable V16 Portable" OFF)

option(USE_V16_AVX "Enable V16 AVX" OFF)

option(USE_V16_AVX2 "Enable V16 AVX2" OFF)

option(USE_ADVANCE_P_AUTOVEC "Enable Explicit Autovec" OFF)

option(ENABLE_OPENSSL "Enable OpenSSL support for checksums" OFF)

option(VPIC_PRINT_MORE_DIGITS "Print more digits in VPIC timer info" OFF)

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
