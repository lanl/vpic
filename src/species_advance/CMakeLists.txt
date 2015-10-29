#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

set(species_advance_HEADERS
  species_advance.h
  standard/spa_private.h
  PARENT_SCOPE
)

set(species_advance_SOURCES
  standard/advance_p.cc
  standard/center_p.cc
  standard/energy_p.cc
  standard/move_p.cc
  standard/rho_p.cc
  standard/uncenter_p.cc
  species_advance.c
  standard/hydro_p.c
  standard/sort_p.c
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