#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

set(grid_HEADERS
  grid.h
  PARENT_SCOPE
)

set(grid_SOURCES
  grid_comm.c
  grid_structors.c
  ops.c
  partition.c
  PARENT_SCOPE
)

#if(ENABLE_UNIT_TESTS)
#  cinch_add_unit(casename
#    SOURCES testfile.cc
#    LIBRARIES list
#    INCLUDES list
#    POLICY MPI
#    THREADS 1 2 4
#  )
#endif(ENABLE_UNIT_TESTS)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#