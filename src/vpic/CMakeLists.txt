#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

set(vpic_HEADERS
  dumpmacros.h
  vpic.h
  PARENT_SCOPE
)

set(vpic_SOURCES
  advance.cc
  diagnostics.cc
  dump.cc
  initialize.cc
  misc.cc
  vpic.cc
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