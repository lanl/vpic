################################################################################
# SPE configure.ac
################################################################################
AC_INIT(VPIC_SPE, 3.1.2.1, bergen@lanl.gov)
AC_CONFIG_AUX_DIR(config)

################################################################################
# automake and libtool initialization
################################################################################
AM_INIT_AUTOMAKE

################################################################################
# list targets for output
################################################################################
AC_CONFIG_FILES(Makefile bin/Makefile)

#AC_CONFIG_COMMANDS(setup, [make setup])

################################################################################
# append prefix with spe
################################################################################
prefix=$prefix/spe

################################################################################
# set absolute paths
################################################################################
host_build_directory=`(cd ../..; pwd)`
host_ac_directory=`(cd $srcdir/../..; pwd)`
AC_SUBST(HOST_BUILDDIR, [$host_build_directory])
AC_SUBST(SPE_BUILDDIR, [$host_build_directory/cell/spe])
AC_SUBST(HOST_ACDIR, [$host_ac_directory])
AC_SUBST(SPE_ACDIR, [$host_ac_directory/cell/spe])

################################################################################
# machine file specification
################################################################################
CCS_WITH_MACHINE("$host_ac_directory/machine")

################################################################################
# with options
################################################################################
CCS_WITH_SPU_COMPILER
CCS_WITH_OPT_LEVEL
CCS_WITH_IBM_SIMDMATH
CCS_WITH_ADDRESSING
#CCS_WITH_LIBMISC_SPU

################################################################################
# C/C++ configuration
################################################################################
AC_PROG_CC
CCS_TUNE_SCC
AC_PROG_CXX
CCS_TUNE_SCXX

################################################################################
# Flags
################################################################################
EXTRA_CPPFLAGS="$IBM_SIMDMATH_CPPFLAGS $LIBMISC_SPU_CPPFLAGS"
AC_SUBST(EXTRA_CPPFLAGS)
EXTRA_LDFLAGS="$IBM_SIMDMATH_LDFLAGS $LIBMISC_SPU_LDFLAGS"
AC_SUBST(EXTRA_LDFLAGS)
EXTRA_LIBS="$IBM_SIMDMATH_LIBS $LIBMISC_SPU_LIBS"
AC_SUBST(EXTRA_LIBS)

################################################################################
# start output
################################################################################
AC_OUTPUT
