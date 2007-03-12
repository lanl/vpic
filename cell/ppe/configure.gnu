#! /bin/bash
################################################################################
# shell script wrapper to call configure for cross compilation
################################################################################

# just call configure with the correct cross-compile flag and passing
# the input given by the main configure script
`echo ${!#} | sed 's/--srcdir=//g'`/configure --host=powerpc "$@"
