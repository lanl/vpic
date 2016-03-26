#~----------------------------------------------------------------------------~#
# VPIC package configuration
#~----------------------------------------------------------------------------~#

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
# Add VPIC unit test policy
#------------------------------------------------------------------------------#

cinch_add_unit_execution_policy(VPIC
  ${CMAKE_SOURCE_DIR}/utilities/gtest-vpic.cc
  FLAGS ${MPI_${MPI_LANGUAGE}_COMPILE_FLAGS}
  INCLUDES ${MPI_${MPI_LANGUAGE}_INCLUDE_PATH}
  LIBRARIES ${MPI_${MPI_LANGUAGE}_LIBRARIES}
  EXEC ${MPIEXEC}
  EXEC_THREADS ${MPIEXEC_NUMPROC_FLAG})

#------------------------------------------------------------------------------#
# Add VPIC integrated test mechanism
#------------------------------------------------------------------------------#

if(ENABLE_INTEGRATED_TESTS)
  enable_testing()
  add_subdirectory(test/integrated)
endif(ENABLE_INTEGRATED_TESTS)

#------------------------------------------------------------------------------#
# Act on build options set in project.cmake
#------------------------------------------------------------------------------#
if(USE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)

set(USE_V4)

if(USE_V4_PORTABLE)
  add_definitions(-DUSE_V4_PORTABLE)
  set(USE_V4 True)
endif(USE_V4_PORTABLE)

if(USE_V4_SSE)
  add_definitions(-DUSE_V4_SSE)
  set(USE_V4 True)
endif(USE_V4_SSE)

if(USE_V4_AVX)
  add_definitions(-DUSE_V4_AVX)
  set(USE_V4 True)
endif(USE_V4_AVX)

if(USE_V4_AVX2)
  add_definitions(-DUSE_V4_AVX2)
  set(USE_V4 True)
endif(USE_V4_AVX2)

if(USE_V4_ALTIVEC)
  add_definitions(-DUSE_V4_ALTIVEC)
  set(USE_V4 True)
endif(USE_V4_ALTIVEC)

set(USE_V8)

if(USE_V8_PORTABLE)
  add_definitions(-DUSE_V8_PORTABLE)
  set(USE_V8 True)
endif(USE_V8_PORTABLE)

if(USE_V8_AVX)
  add_definitions(-DUSE_V8_AVX)
  set(USE_V8 True)
endif(USE_V8_AVX)

if(USE_V8_AVX2)
  add_definitions(-DUSE_V8_AVX2)
  set(USE_V8 True)
endif(USE_V8_AVX2)

if(ENABLE_OPENSSL)
  add_definitions(-DENABLE_OPENSSL)
endif(ENABLE_OPENSSL)

#------------------------------------------------------------------------------#
# include cmake hacks
#------------------------------------------------------------------------------#

#include(config/hacks.cmake)

#------------------------------------------------------------------------------#
# Handle vpic compile script last
#------------------------------------------------------------------------------#

if(BUILD_SHARED_LIBS)
    set(VPIC_CXX_FLAGS "-rdynamic ${VPIC_CXX_FLAGS}")
endif(BUILD_SHARED_LIBS)

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

#~---------------------------------------------------------------------------~-#
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
