# ----------------------------------------------------------------------
# Cray systems (Cielo, Cielito, Trinitite, Trinity)
#
# This file is courtesy of Kelly Thompson (KT), kgt@lanl.gov
#
# cmake -C /path/to/vpic/arch/CrayConfig.cmake [options] /path/to/vpic
# or
# ccmake -C /path/to/vpic/arch/CrayConfig.cmake /path/to/vpic
# ----------------------------------------------------------------------

#
# Compiler settings
#
# Option: Add -craype-verbose to see actual intel compiler invocation
# with all flags.

foreach( lang C CXX Fortran )
  # forcing -dynamic allows cmake's intitial compiler check to pass.
  # adding -craype-verbose shows compiler flags when % make VERBOSE=1.
  set( CMAKE_${lang}_FLAGS "-dynamic -craype-verbose" CACHE STRING
    "Use -dynamic when checking if the compiler works. Add -craype-verbose." )
endforeach()
set( CMAKE_EXE_LINKER_FLAGS "-dynamic" CACHE STRING
  "Extra flags for linking executables")

set( CMAKE_C_COMPILER       cc  CACHE STRING "C compiler" )
set( CMAKE_CXX_COMPILER     CC  CACHE STRING "C++ compiler" )
set( CMAKE_Fortran_COMPILER ftn CACHE STRING "Fortran compiler" )

#
# MPI settings
#

# Setup MPI for Cray MPT.
set( MPIEXEC "aprun" CACHE STRING "Application Level Placement Scheduler (mpirun)")
set( MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING
  "mpirun flag used to specify the number of processors to use")

# Setting these prevents FindMPI from autodetecting the mpi compiler
# wrappers and compile flags.  These are taken care of by the cray
# compiler wrappers set above.
foreach( lang C CXX Fortran )
  set( MPI_${lang}_LIBRARIES    "" CACHE STRING "Let the compile wrapper choose.")
  set( MPI_${lang}_INCLUDE_PATH "" CACHE STRING "Let the compile wrapper choose.")
  set( MPI_${lang}_COMPILER     "" CACHE STRING "Use compile wrapper.")
endforeach()
