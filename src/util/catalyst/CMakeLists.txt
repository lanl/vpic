project(VPICAdaptor)
cmake_minimum_required(VERSION 2.6)

# probably a hack since I haven't looked at how to include directories
# properly
include_directories(../../vpic)
include_directories(../../grid)
include_directories(../../util)
include_directories(../../sf_interface)

include_directories( ${VTK_INCLUDE_DIRS})

# Test to figure out if we have the new vtkDataArray
# API or the old one. VPIC+Catalyst will work properly
# with both but the old API has been deprecated and
# will cause compiler warnings.
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/NewVTKArrayAPI.cxx
"
// Code for a CMake try_compile to see if we have newer
// VTK features

#include <vtkFloatArray.h>

int main()
{
  vtkFloatArray* array = vtkFloatArray::New();
  array->SetNumberOfTuples(1);
  float one = 1;
  array->SetTypedTuple(0, &one);
  array->Delete();
  return 0;
}"
)
try_compile(NEW_VTK_ARRAY_API
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/NewVTKArrayAPI.cxx
  CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${VTK_INCLUDE_DIRS}"
  LINK_LIBRARIES vtkCommonCore
  OUTPUT_VARIABLE ov
)
if(NEW_VTK_ARRAY_API)
  add_definitions(-DNEW_VTK_ARRAY_API)
endif(NEW_VTK_ARRAY_API)
# Done with the VTK's newer array API

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

install(TARGETS VPICAdaptor
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
