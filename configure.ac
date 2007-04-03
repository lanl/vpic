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
	lib/Makefile \
	bin/Makefile \
	cell/Makefile)
AC_CONFIG_COMMANDS(setup, [make setup])

################################################################################
# machine file specification
################################################################################
CCS_WITH_MACHINE

################################################################################
# enable options
################################################################################
CCS_ENABLE_CELL
CCS_ENABLE_HOST
CCS_ENABLE_OPT_LEVEL
CCS_ENABLE_DMP
CCS_ENABLE_SCALAR

################################################################################
# with options
################################################################################
CCS_WITH_ADDRESSING
CCS_WITH_CELL_INSTDIR

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
	AC_SUBST(CELLSUBDIRS, ["spe ppe"])
	AC_CONFIG_SUBDIRS([cell/spe cell/ppe])
fi

# enable/disable distributed memory parallelism
# this is also used for the openmpi host/cell interface
if test "$ENABLE_DMP" = "1" ; then
	CCS_WITH_MPI

	AC_SUBST(PARALLEL, [parallel])

	EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $MPI_CPPFLAGS"
	EXTRA_LDFLAGS="$EXTRA_LDFLAGS $MPI_LDFLAGS"
	EXTRA_LIBS="$EXTRA_LIBS $MPI_LIBS"

	AC_SUBST(EXTRA_CPPFLAGS)
	AC_SUBST(EXTRA_LDFLAGS)
	AC_SUBST(EXTRA_LIBS)

	# define needed for openmpi
	AC_DEFINE(OMPI_SKIP_MPICXX)
fi

################################################################################
# enable/disable host build
################################################################################
if test "$ENABLE_HOST" = "1" ; then
	AC_SUBST(HOSTDIRS, ["lib bin"])
fi

################################################################################
# start output
################################################################################
AC_OUTPUT
