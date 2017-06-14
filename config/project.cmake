#~----------------------------------------------------------------------------~#
# VPIC project configuration
#~----------------------------------------------------------------------------~#

project(vpic)

#------------------------------------------------------------------------------#
# Set project defaults
#------------------------------------------------------------------------------#

# As this is an MPI code, we need MPI on.
set(ENABLE_MPI True CACHE BOOL "Enable MPI" FORCE)
set(ENABLE_MPI_CXX_BINDINGS False CACHE BOOL "Enable MPI C++ Bindings" FORCE)

# Hide the ENABLE_MPI option from users, so they don't try and disable it
mark_as_advanced(ENABLE_MPI)

#------------------------------------------------------------------------------#
# Add build options
#------------------------------------------------------------------------------#

option(ENABLE_INTEGRATED_TESTS "enable integrated tests" OFF)

option(USE_OPENMP "Use OpenMP" OFF)

option(USE_V4_ALTIVEC "Enable V4 Altivec" OFF)

option(USE_V4_PORTABLE "Enable V4 Portable" OFF)

option(USE_V4_SSE "Enable V4 SSE" OFF)

option(ENABLE_OPENSSL "Enable OpenSSL support for checksums" OFF)

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
