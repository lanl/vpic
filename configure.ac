################################################################################
# host configure.ac
################################################################################
AC_INIT(VPIC, 3.1.2.1, bergen@lanl.gov)
AC_CONFIG_AUX_DIR(config)
VERSION_INFORMATION="-release 3.1.2.1 "
AC_SUBST(VERSION_INFORMATION)

################################################################################
# automake and libtool initialization
################################################################################
AM_INIT_AUTOMAKE
AC_PROG_LIBTOOL

################################################################################
# list targets for output and commands to run
################################################################################
AC_CONFIG_FILES(Makefile \
	buildscript \
	lib/Makefile \
	mprelay/Makefile \
	accelerator/Makefile \
	cell/Makefile)
#	bin/Makefile \
AC_CONFIG_COMMANDS(setup, [make setup])

################################################################################
# need this so that ppe build will create build script correctly
################################################################################
AC_SUBST(HOST_ACDIR, [$srcdir])
AC_SUBST(HOST_BUILDDIR, [`pwd`])

################################################################################
# machine file specification
################################################################################
CCS_WITH_MACHINE("$srcdir/machine")

################################################################################
# enable options
################################################################################
CCS_ENABLE_CELL
CCS_ENABLE_HOST

################################################################################
# with options
################################################################################
CCS_WITH_BUILDSTYLE
CCS_WITH_ADDRESSING
CCS_WITH_CPUMODEL
CCS_WITH_EXTENSION
CCS_WITH_OPT_LEVEL


################################################################################
# extra flags
################################################################################
CCS_EXTRA_FLAGS

################################################################################
# C/C++ configuration
################################################################################
AC_PROG_CXX
AC_PROG_CC
CCS_SIMD_ISA
CCS_TUNE_CXX
CCS_TUNE_CC

################################################################################
# determine build style
################################################################################
case "$buildstyle" in
	standard)
		echo "standard build"
		AC_SUBST(HOSTDIRS, ["lib"])
		AC_SUBST(ARCH_BUILD_PATH, [.])
	;;
	test_relay)
		echo "test_relay build"
		AC_SUBST(HOSTDIRS, ["lib mprelay accelerator"])
		AC_SUBST(ARCH_BUILD_PATH, [.])
		AC_DEFINE(USE_MPRELAY, [1])
	;;
	hybrid_relay_cell)
		AC_SUBST(CELLDIR, [cell])
		AC_SUBST(CELLSUBDIRS, ["ppe"])
		AC_CONFIG_SUBDIRS([cell/ppe])
		AC_SUBST(HOSTDIRS, [mprelay])
		AC_SUBST(ARCH_BUILD_PATH, [cell/ppe])
		AC_DEFINE(USE_MPRELAY, [1])
	;;
esac

# this does auto-detection of mpi libraries using
# the C compiler and possibly the openmpi mpicc program
CCS_WITH_MPI(C)

EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $MPI_CPPFLAGS"
EXTRA_LDFLAGS="$EXTRA_LDFLAGS $MPI_LDFLAGS"
EXTRA_LIBS="$EXTRA_LIBS $MPI_LIBS"

# define needed for openmpi
AC_DEFINE(OMPI_SKIP_MPICXX)

AC_SUBST(EXTRA_CPPFLAGS)
AC_SUBST(EXTRA_LDFLAGS)
AC_SUBST(EXTRA_LIBS)

################################################################################
# start output
################################################################################
AC_OUTPUT
