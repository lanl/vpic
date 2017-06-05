#~----------------------------------------------------------------------------~#
# VPIC project configuration
#~----------------------------------------------------------------------------~#

cinch_minimum_required(2.0)

project(vpic)

#------------------------------------------------------------------------------#
# If a C++11 compiler is available, then set the appropriate flags
#------------------------------------------------------------------------------#

include(cxx11)

check_for_cxx11_compiler(CXX11_COMPILER)

if(CXX11_COMPILER)
  enable_cxx11()
else()
  message(FATAL_ERROR "C++11 compatible compiler not found")
endif()

#------------------------------------------------------------------------------#
# Set C flags
#------------------------------------------------------------------------------#

include(c99)

check_for_c99_compiler(CXX11_COMPILER)

if(CXX11_COMPILER)
  enable_c99()
else()
  message(FATAL_ERROR "C99 compatible compiler not found")
endif()

#------------------------------------------------------------------------------#
# Set project defaults
#
# These have to occur before cinch_load_extras() because they
# override some of cinch's default options 
#------------------------------------------------------------------------------#

set(ENABLE_MPI True CACHE BOOL "Enable MPI" FORCE)
set(ENABLE_MPI_CXX_BINDINGS False CACHE BOOL "Enable MPI C++ Bindings" FORCE)

#------------------------------------------------------------------------------#
# Set header suffix
#------------------------------------------------------------------------------#

set(CINCH_HEADER_SUFFIXES "\\.h")

#------------------------------------------------------------------------------#
# Load extra functionality in cinch
#------------------------------------------------------------------------------#

cinch_load_extras()

#------------------------------------------------------------------------------#
# Add build options
#------------------------------------------------------------------------------#

option(ENABLE_INTEGRATED_TESTS "enable integrated tests" OFF)

option(USE_V4_ALTIVEC "Enable V4 Altivec" OFF)

option(USE_V4_PORTABLE "Enable V4 Portable" OFF)

option(USE_V4_SSE "Enable V4 SSE" OFF)

option(ENABLE_OPENSSL "Enable OpenSSL support for checksums" OFF)

if ( USE_OPENMP ) 
  message( FATAL_ERROR "\
    This option has moved to ENABLE_OPENMP \
    to avoid duplication \
  " )
endif()

#------------------------------------------------------------------------------#
# Add MPI includes
#------------------------------------------------------------------------------#

# FIXME: this should be target-specific, but this will require a change
# to Cinch's add_library_target function.
include_directories(${MPI_C_INCLUDE_PATH})

#------------------------------------------------------------------------------#
# Create include and link aggregates
#
# NOTE: These must be set before creating the compile scripts below.
#------------------------------------------------------------------------------#

set(VPIC_CPPFLAGS)
if(MPI_CPPFLAGS)
  string(REPLACE ";" " " string_cppflags "${MPI_CPPFLAGS}")
  set(VPIC_CPPFLAGS "${string_cppflags}")
endif(MPI_CPPFLAGS)

string(REPLACE ";" " -I" string_includes "${MPI_C_INCLUDE_PATH}")
if(NOT ${string_includes} STREQUAL "")
  set(VPIC_CXX_FLAGS "-I${string_includes} ${MPI_C_LINK_FLAGS}")
endif(NOT ${string_includes} STREQUAL "")

# Add Debug flags to VPIC_CXX_FLAGS
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")

# Add RelWithDebInfo flags to VPIC_CXX_FLAGS
if("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")

string(REPLACE ";" " " string_libraries "${MPI_C_LIBRARIES}")
set(VPIC_CXX_LIBRARIES "${string_libraries}")

#------------------------------------------------------------------------------#
# OpenSSL
#------------------------------------------------------------------------------#

if(ENABLE_OPENSSL)
  find_package(OpenSSL REQUIRED)

  include_directories(${OPENSSL_INCLUDE_DIR})
  string(REPLACE ";" " " string_libraries "${OPENSSL_LIBRARIES}")
  set(VPIC_CXX_LIBRARIES "${VPIC_CXX_LIBRARIES} ${string_libraries}")
endif(ENABLE_OPENSSL)

#------------------------------------------------------------------------------#
# Act on build options set in project.cmake
#------------------------------------------------------------------------------#

set(USE_V4)
if(USE_V4_ALTIVEC)
  add_definitions(-DUSE_V4_ALTIVEC)
  set(USE_V4 True)
endif(USE_V4_ALTIVEC)

if(USE_V4_PORTABLE)
  add_definitions(-DUSE_V4_PORTABLE)
  set(USE_V4 True)
endif(USE_V4_PORTABLE)

if(USE_V4_SSE)
  add_definitions(-DUSE_V4_SSE)
  set(USE_V4 True)
endif(USE_V4_SSE)

if(ENABLE_OPENSSL)
  add_definitions(-DENABLE_OPENSSL)
endif(ENABLE_OPENSSL)

#------------------------------------------------------------------------------#
# Handle vpic compile script last
#------------------------------------------------------------------------------#

if(BUILD_SHARED_LIBS)
    set(VPIC_CXX_FLAGS "-rdynamic ${VPIC_CXX_FLAGS}")
endif(BUILD_SHARED_LIBS)

if(ENABLE_COVERAGE_BUILD)
    set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} --coverage")
endif(ENABLE_COVERAGE_BUILD)

# process Makefile.run.in to get a simple Makefile.run for a run. Points to
# local built exe wrapper, and has example deck/platform.
configure_file(${CMAKE_SOURCE_DIR}/sample/Makefile.run.in
  ${CMAKE_BINARY_DIR}/bin/Makefile.run)

# install script
configure_file(${CMAKE_SOURCE_DIR}/bin/vpic.in
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vpic-install)
install(FILES ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vpic-install
  DESTINATION bin
  RENAME vpic
  PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
    )

install(FILES ${CMAKE_SOURCE_DIR}/deck/main.cc
  DESTINATION share/vpic)
install(FILES ${CMAKE_SOURCE_DIR}/deck/wrapper.cc
  DESTINATION share/vpic)

# local script
configure_file(${CMAKE_SOURCE_DIR}/bin/vpic-local.in
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vpic)

file(COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vpic
  DESTINATION ${CMAKE_BINARY_DIR}/bin
  FILE_PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
)

#------------------------------------------------------------------------------#
# Add library target
#------------------------------------------------------------------------------#

cinch_add_library_target(vpic src)
target_compile_options( vpic PRIVATE ${MPI_C_COMPILE_FLAGS} )

#------------------------------------------------------------------------------#
# Add VPIC integrated test mechanism
#------------------------------------------------------------------------------#

if(ENABLE_INTEGRATED_TESTS)
  enable_testing()
  add_subdirectory(test/integrated)
endif(ENABLE_INTEGRATED_TESTS)

#~---------------------------------------------------------------------------~-#
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
