#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

set(sf_interface_HEADERS
  sf_interface.h
  sf_interface_private.h
  PARENT_SCOPE
)

set(sf_interface_SOURCES
  interpolator_array.cc
  reduce_accumulators.cc
  unload_accumulator.cc
  accumulator_array.c
  clear_accumulators.c
  hydro_array.c
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