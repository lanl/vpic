project(VPICAdaptor)
cmake_minimum_required(VERSION 2.6)
#include_directories(${CoProcessor_SOURCE_DIR})

# ACBAUER a hack since I haven't looked at how to include directories
# properly
include_directories(../../vpic)
include_directories(../../grid)
include_directories(../../util)
include_directories(../../sf_interface)

set(VPICAdaptor_SRC
  VPICAdaptor.cxx
  #vtkVPICPipeline.cxx
  #vtkVPICContourPipeline.cxx
  #vtkVPICDataSetPipeline.cxx
  #vtkVPICPistonContourPipeline.cxx
)
add_library(VPICAdaptor ${VPICAdaptor_SRC})
target_link_libraries(
  VPICAdaptor vtkPVPythonCatalyst
)

# Install
#----------------------------------------------------------------------------------------------------
# Libraries
#IF(NOT PV_INSTALL_NO_LIBRARIES)

#  INSTALL(TARGETS VPICAdaptor
#      RUNTIME DESTINATION ${PV_INSTALL_BIN_DIR} COMPONENT Runtime
#      LIBRARY DESTINATION ${PV_INSTALL_LIB_DIR} COMPONENT Runtime
#      ARCHIVE DESTINATION ${PV_INSTALL_LIB_DIR} COMPONENT Development)
#ENDIF(NOT PV_INSTALL_NO_LIBRARIES)
#
# Development
#IF (NOT PV_INSTALL_NO_DEVELOPMENT)
#  GLOB_INSTALL_DEVELOPMENT(
#      ${CMAKE_CURRENT_SOURCE_DIR}
#      ${PV_INSTALL_INCLUDE_DIR}
#      "*.h;*.hxx;*.txx")
#ENDIF (NOT PV_INSTALL_NO_DEVELOPMENT)

