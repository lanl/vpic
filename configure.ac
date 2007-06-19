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
	build.script \
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
CCS_ENABLE_DMP
CCS_ENABLE_SCALAR

################################################################################
# with options
################################################################################
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
# enable/disable cell build
################################################################################
if test "$ENABLE_CELL" = "1" ; then
	AC_SUBST(CELLDIR, [cell])
#	AC_SUBST(CELLSUBDIRS, ["spe ppe"])
	AC_SUBST(CELLSUBDIRS, ["ppe"])
#	AC_CONFIG_SUBDIRS([cell/spe cell/ppe])
	AC_CONFIG_SUBDIRS([cell/ppe])
	AC_SUBST(HOSTDIRS, [mprelay])
	AC_SUBST(ARCH_BUILD_PATH, [cell/ppe])
else
#	AC_SUBST(HOSTDIRS, [lib])
	AC_SUBST(HOSTDIRS, ["lib mprelay accelerator"])
	AC_SUBST(ARCH_BUILD_PATH, [.])
fi

# enable/disable distributed memory parallelism
# this is also used for the openmpi host/cell interface
if test "$ENABLE_DMP" = "1" ; then
	# this does auto-detection of mpi libraries using
	# the C compiler and possibly the openmpi mpicc program
	CCS_WITH_MPI(C)

	AC_SUBST(PARALLEL, [parallel])

	EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $MPI_CPPFLAGS"
	EXTRA_LDFLAGS="$EXTRA_LDFLAGS $MPI_LDFLAGS"
	EXTRA_LIBS="$EXTRA_LIBS $MPI_LIBS"

	# define needed for openmpi
	AC_DEFINE(OMPI_SKIP_MPICXX)
fi

AC_SUBST(EXTRA_CPPFLAGS)
AC_SUBST(EXTRA_LDFLAGS)
AC_SUBST(EXTRA_LIBS)

################################################################################
# start output
################################################################################
AC_OUTPUT
