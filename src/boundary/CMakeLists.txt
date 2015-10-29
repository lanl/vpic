#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

set(boundary_HEADERS
  boundary.h
  boundary_private.h
  PARENT_SCOPE
)

set(boundary_SOURCES
  absorb_tally.c
  boundary.c
  boundary_p.cc
  link.c
  maxwellian_reflux.c
  PARENT_SCOPE
)

#if(ENABLE_UNIT_TESTS)
#    cinch_add_unit(casename testfile.cc)
#endif(ENABLE_UNIT_TESTS)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
