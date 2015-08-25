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
set(VPIC_CXX_FLAGS "-I${string_includes} ${MPI_C_LINK_FLAGS}")

string(REPLACE ";" " " string_libraries "${MPI_C_LIBRARIES}")
set(VPIC_CXX_LIBRARIES "${string_libraries}")

#------------------------------------------------------------------------------#
# Handle vpic compile script
#------------------------------------------------------------------------------#

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
