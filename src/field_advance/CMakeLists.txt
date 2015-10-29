#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

set(field_advance_HEADERS
  field_advance.h
  field_advance_private.h
  standard/sfa_private.h
  PARENT_SCOPE
)

set(field_advance_SOURCES
  standard/advance_b.cc
  standard/advance_e.cc
  standard/clean_div_b.cc
  standard/compute_curl_b.cc
  standard/compute_div_b_err.cc
  standard/vacuum_advance_e.cc
  standard/vacuum_compute_curl_b.cc
  field_advance.c
  standard/clean_div_e.c
  standard/compute_div_e_err.c
  standard/compute_rhob.c
  standard/compute_rms_div_b_err.c
  standard/compute_rms_div_e_err.c
  standard/energy_f.c
  standard/local.c
  standard/remote.c
  standard/sfa.c
  standard/vacuum_clean_div_e.c
  standard/vacuum_compute_div_e_err.c
  standard/vacuum_compute_rhob.c
  standard/vacuum_energy_f.c
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