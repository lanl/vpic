#------------------------------------------------------------------------------#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.1)

#~----------------------------------------------------------------------------~#
# VPIC project configuration
#~----------------------------------------------------------------------------~#

project(vpic)

#------------------------------------------------------------------------------#
# Quiet newer versions of cmake about MPI_ROOT
#------------------------------------------------------------------------------#
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

#------------------------------------------------------------------------------#
# If a C++11 compiler is available, then set the appropriate flags
#------------------------------------------------------------------------------#

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#------------------------------------------------------------------------------#
# Set C flags
#------------------------------------------------------------------------------#

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

#------------------------------------------------------------------------------#
# MPI detection
#------------------------------------------------------------------------------#

find_package(MPI REQUIRED)
list(APPEND MPI_CXX_COMPILE_FLAGS -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX )
include_directories(${MPI_C_INCLUDE_PATH})

#------------------------------------------------------------------------------#
# Add build options
#------------------------------------------------------------------------------#

option(ENABLE_INTEGRATED_TESTS "Enable integrated tests" OFF)

option(ENABLE_UNIT_TESTS "Enable unit tests" OFF)

option(ENABLE_LONG_TESTS "Enable longer running tests, often skipped for most CI builds" OFF)

option(ENABLE_PERFORMANCE_TESTS "Enable performance tests" OFF)

option(USE_OPENMP "Use OpenMP" OFF)

option(USE_PTHREADS "Use Pthreads" ON)

option(USE_V4_ALTIVEC "Enable V4 Altivec" OFF)

option(USE_V4_NEON "Enable V4 NEON" OFF)

option(USE_V4_PORTABLE "Enable V4 Portable" OFF)

option(USE_V4_SSE "Enable V4 SSE" OFF)

option(USE_V4_AVX "Enable V4 AVX" OFF)

option(USE_V4_AVX2 "Enable V4 AVX2" OFF)

option(USE_V8_PORTABLE "Enable V8 Portable" OFF)

option(USE_V8_AVX "Enable V8 AVX" OFF)

option(USE_V8_AVX2 "Enable V8 AVX2" OFF)

option(USE_V16_PORTABLE "Enable V16 Portable" OFF)

option(USE_V16_AVX512 "Enable V16 AVX512" OFF)

option(USE_LEGACY_SORT "Enable Legacy Sort Implementation" OFF)

#option(USE_ADVANCE_P_AUTOVEC "Enable Explicit Autovec" OFF)

option(VPIC_PRINT_MORE_DIGITS "Print more digits in VPIC timer info" OFF)

option(ENABLE_OPENSSL "Enable OpenSSL support for checksums" OFF)

option(DISABLE_DYNAMIC_RESIZING "Prevent particle arrays from dynamically resizing during a run" OFF)

option(EXIT_ON_LOST_MOVER "EXIT if we run out of mover space during the particle push (the default is WARN)" OFF)

# option to set minimum number of particles
set(SET_MIN_NUM_PARTICLES AUTO CACHE STRING "Select minimum number of particles to use, if using dynamic particle array resizing")


#------------------------------------------------------------------------------#
# Create include and link aggregates
#
# NOTE: These must be set before creating the compile scripts below.
#------------------------------------------------------------------------------#

string(REPLACE ";" " -I" string_includes "${MPI_C_INCLUDE_PATH}")
if(NOT ${string_includes} STREQUAL "")
  set(VPIC_CXX_FLAGS "-I${string_includes} ${MPI_C_LINK_FLAGS}")
endif(NOT ${string_includes} STREQUAL "")
set(VPIC_CXX_FLAGS "${CMAKE_CXX11_STANDARD_COMPILE_OPTION} ${VPIC_CXX_FLAGS}")

# Add Debug flags to VPIC_CXX_FLAGS
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")

# Add RelWithDebInfo flags to VPIC_CXX_FLAGS
if("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")

string(REPLACE ";" " " string_libraries "${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES}")
set(VPIC_CXX_LIBRARIES "${string_libraries}")

if(EXIT_ON_LOST_MOVER)
    add_definitions(-DEXIT_ON_LOST_MOVER)
endif(EXIT_ON_LOST_MOVER)

if(DISABLE_DYNAMIC_RESIZING)
  add_definitions(-DDISABLE_DYNAMIC_RESIZING)
endif(DISABLE_DYNAMIC_RESIZING)

if(NOT SET_MIN_NUM_PARTICLES STREQUAL "AUTO")
    add_definitions(-DMIN_NP=${SET_MIN_NUM_PARTICLES})
endif()

#------------------------------------------------------------------------------#
# OpenSSL
#------------------------------------------------------------------------------#

if(ENABLE_OPENSSL)
  find_package(OpenSSL REQUIRED)

  include_directories(${OPENSSL_INCLUDE_DIR})
  string(REPLACE ";" " " string_libraries "${OPENSSL_LIBRARIES}")
  set(VPIC_CXX_LIBRARIES "${VPIC_CXX_LIBRARIES} ${string_libraries}")
endif(ENABLE_OPENSSL)

find_package(Threads REQUIRED)

#------------------------------------------------------------------------------#
# Act on build options set in project.cmake
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Add options for building with the legacy particle sort implementation.
#------------------------------------------------------------------------------#

if(USE_LEGACY_SORT)
  add_definitions(-DVPIC_USE_LEGACY_SORT)
    set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DVPIC_USE_LEGACY_SORT")
endif(USE_LEGACY_SORT)

#------------------------------------------------------------------------------#
# Add options for building with a threading model.
#------------------------------------------------------------------------------#

# We don't want both PTHREADS and OpenMP
if ((USE_PTHREADS) AND (USE_OPENMP))
     message( FATAL_ERROR "Only one threading model can be selected" )
endif()


if(USE_PTHREADS)
  add_definitions(-DVPIC_USE_PTHREADS)
    set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DVPIC_USE_PTHREADS")
endif(USE_PTHREADS)

if(USE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DVPIC_USE_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DVPIC_USE_OPENMP")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)

#------------------------------------------------------------------------------#
# Add options for building with v4 simd vector support.
#------------------------------------------------------------------------------#

set(USE_V4)

if(USE_V4_PORTABLE)
  add_definitions(-DUSE_V4_PORTABLE)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V4_PORTABLE")
  set(USE_V4 True)
endif(USE_V4_PORTABLE)

if(USE_V4_SSE)
  add_definitions(-DUSE_V4_SSE)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V4_SSE")
  set(USE_V4 True)
endif(USE_V4_SSE)

if(USE_V4_AVX)
  add_definitions(-DUSE_V4_AVX)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V4_AVX")
  set(USE_V4 True)
endif(USE_V4_AVX)

if(USE_V4_AVX2)
  add_definitions(-DUSE_V4_AVX2)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V4_AVX2")
  set(USE_V4 True)
endif(USE_V4_AVX2)

if(USE_V4_ALTIVEC)
  add_definitions(-DUSE_V4_ALTIVEC)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V4_ALTIVEC")
  set(USE_V4 True)
endif(USE_V4_ALTIVEC)

if(USE_V4_NEON)
  add_definitions(-DUSE_V4_NEON)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V4_NEON")
  set(USE_V4 True)
endif(USE_V4_NEON)

#------------------------------------------------------------------------------#
# Add options for building with v8 simd vector support.
#------------------------------------------------------------------------------#

set(USE_V8)

if(USE_V8_PORTABLE)
  add_definitions(-DUSE_V8_PORTABLE)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V8_PORTABLE")
  set(USE_V8 True)
endif(USE_V8_PORTABLE)

if(USE_V8_AVX)
  add_definitions(-DUSE_V8_AVX)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V8_AVX")
  set(USE_V8 True)
endif(USE_V8_AVX)

if(USE_V8_AVX2)
  add_definitions(-DUSE_V8_AVX2)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V8_AVX2")
  set(USE_V8 True)
endif(USE_V8_AVX2)

#------------------------------------------------------------------------------#
# Add options for building with v16 simd vector support.
#------------------------------------------------------------------------------#

set(USE_V16)

if(USE_V16_PORTABLE)
  add_definitions(-DUSE_V16_PORTABLE)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V16_PORTABLE")
  set(USE_V16 True)
endif(USE_V16_PORTABLE)

if(USE_V16_AVX512)
  add_definitions(-DUSE_V16_AVX512)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DUSE_V16_AVX512")
  set(USE_V16 True)
endif(USE_V16_AVX512)

# TODO: Can we improve the way this is done so it's detection of a positive not
# multiple negatives?
if (NOT USE_V4 AND NOT USE_V8 AND NOT USE_V16)
    # This option is intended to indicate that we're not using any of the
    # specialized intrinsics code paths, and we can expect numerical answer to
    # be the same as the "reference" (serial/auto-vec) implementation. This can
    # be especially useful when trying to test small differences in expressions
    # and we want bitwise answers to a known implementation
    set(NO_EXPLICIT_VECTOR True)
endif()

# Detect if the user enabled V8/V16 but not v8/v4. If they're not enabled, tell
# the user instead of giving them a confusing compile error inside a macro
if (USE_V16 AND (NOT USE_V8 OR NOT USE_V4))
    message(FATAL_ERROR
        "You requested V16 but did not also enable both V4 and V8, aborting.")
endif()
if (USE_V8 AND NOT USE_V4)
    message(FATAL_ERROR "You requested V8 but did not also enable V4, aborting.")
endif()


#------------------------------------------------------------------------------#
# Add options for building with explicit autovec support.
#------------------------------------------------------------------------------#

#if(USE_ADVANCE_P_AUTOVEC)
  #add_definitions(-DADVANCE_P_AUTOVEC)
  #set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DADVANCE_P_AUTOVEC")
  #set(ADVANCE_P_AUTOVEC True)
#endif(USE_ADVANCE_P_AUTOVEC)

#------------------------------------------------------------------------------#
# Miscellaneous options.
#------------------------------------------------------------------------------#

if(ENABLE_OPENSSL)
  add_definitions(-DENABLE_OPENSSL)
endif(ENABLE_OPENSSL)

if(VPIC_PRINT_MORE_DIGITS)
  add_definitions(-DVPIC_PRINT_MORE_DIGITS)
  set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} -DVPIC_PRINT_MORE_DIGITS")
endif(VPIC_PRINT_MORE_DIGITS)

#------------------------------------------------------------------------------#
# Handle vpic compile script last.
#------------------------------------------------------------------------------#

# Allow files to be references from the root, (eg #include "src/vpic/vpic.h")
include_directories(${PROJECT_SOURCE_DIR})

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

# Append all defines to VPIC_DEFINES, so it can be seen during input deck building
get_directory_property(ALL_DEFINES DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
#string(REPLACE ";" " -D" EEK "${ALL_DEFINES}")
foreach(d ${ALL_DEFINES})
    set(VPIC_DEFINES "${VPIC_DEFINES} -D${d}")
endforeach()

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

install(FILES deck/main.cc deck/wrapper.cc DESTINATION share/vpic)
install(FILES deck/wrapper.h DESTINATION include/vpic)
install(DIRECTORY src/ DESTINATION include/vpic
         FILES_MATCHING PATTERN "*.h")


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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VPIC_CXX_FLAGS}")

file(GLOB_RECURSE VPIC_SRC src/*.c src/*.cc)
file(GLOB_RECURSE VPIC_NOT_SRC
  src/util/v4/test/v4.cc
  src/util/v8/test/v8.cc
  src/util/v16/test/v16.cc
  src/util/rng/test/rng.cc)
list(REMOVE_ITEM VPIC_SRC ${VPIC_NOT_SRC})
option(NO_LIBVPIC "Don't build a libvpic, but all in one" OFF)
if(NO_LIBVPIC)
  set(VPIC_EXPOSE "INTERFACE")
  add_library(vpic INTERFACE)
else()
  set(VPIC_EXPOSE "PRIVATE")
  add_library(vpic ${VPIC_SRC})
  set(VPIC_SRC)
  install(TARGETS vpic LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
endif()
target_include_directories(vpic INTERFACE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(vpic ${VPIC_EXPOSE} ${MPI_CXX_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${OPENSSL_LIBRARIES} ${CMAKE_DL_LIBS})
target_compile_options(vpic ${VPIC_EXPOSE} ${MPI_C_COMPILE_FLAGS})

macro(build_a_vpic name deck)
  if(NOT EXISTS ${deck})
    message(FATAL_ERROR "Could not find deck '${deck}'")
  endif()
  add_executable(${name} ${CMAKE_SOURCE_DIR}/deck/main.cc ${CMAKE_SOURCE_DIR}/deck/wrapper.cc ${VPIC_SRC})
  target_link_libraries(${name} vpic)
  target_compile_definitions(${name} PRIVATE INPUT_DECK=${deck})
endmacro(build_a_vpic)

if(USER_DECKS)
  foreach(DECK ${USER_DECK})
    get_filename_component(NAME "${DECK}" NAME_WE)
    get_filename_component(DECKDIR "${DECK}" DIRECTORY)
    if(NOT DECKDIR)
      set(DECK "${CMAKE_CURRENT_BINARY_DIR}/${DECK}")
    endif()
    build_a_vpic(${NAME} ${DECK})
  endforeach()
endif()

#------------------------------------------------------------------------------#
# Add VPIC integrated test mechanism
#------------------------------------------------------------------------------#

if(ENABLE_INTEGRATED_TESTS)
  enable_testing()
  add_subdirectory(test/integrated)
endif(ENABLE_INTEGRATED_TESTS)

# Add includes for Catch2
set(CATCH_DIR "./test/include")

if(ENABLE_UNIT_TESTS)
  enable_testing()

  include_directories(${CATCH_DIR})

  # Add tests individually
  if (USE_V4)
    add_executable(v4 src/util/v4/test/v4.cc)
    add_test(NAME v4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./v4)
  endif(USE_V4)

  # RNG tests
  add_executable(rng src/util/rng/test/rng.cc)
  target_link_libraries(rng vpic)
  add_test(NAME rng COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./rng)

  add_subdirectory(test/unit)

endif(ENABLE_UNIT_TESTS)

if(ENABLE_PERFORMANCE_TESTS)
  enable_testing()
  include_directories(${CATCH_DIR})
  add_subdirectory(test/performance)
endif(ENABLE_PERFORMANCE_TESTS)
#~---------------------------------------------------------------------------~-#
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
