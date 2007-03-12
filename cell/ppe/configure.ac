################################################################################
# PPE configure.ac
################################################################################
AC_INIT(VPIC_PPE, 3.1.2.1, bergen@lanl.gov)
AC_CONFIG_AUX_DIR(config)
VERSION_INFORMATION="-release 3.1.2.1 "
AC_SUBST(VERSION_INFORMATION)

################################################################################
# automake and libtool initialization
################################################################################
AM_INIT_AUTOMAKE
AC_PROG_LIBTOOL

################################################################################
# list targets for output
################################################################################
AC_CONFIG_FILES(Makefile \
	config/rules.mk \
	lib/Makefile \
	bin/Makefile)

################################################################################
# I haven't figured out a unique tag defined by the ppu compiler
# so I made one myself
################################################################################
AC_DEFINE([__PPU__])

################################################################################
# append prefix with ppe
################################################################################
prefix=$prefix/ppe

################################################################################
# set absolute path to spe build directory
################################################################################
host_build_directory=`(cd ../..; pwd)`
host_ac_directory=`(cd $srcdir/../..; pwd)`
AC_SUBST(SPE_BUILDDIR, [$host_build_directory/cell/spe])
AC_SUBST(HOST_ACDIR, [$host_ac_directory])
echo "srcdir $srcdir"

################################################################################
# enable options
################################################################################
CCS_ENABLE_DMP

################################################################################
# with options
################################################################################
CCS_WITH_PPU_COMPILER
CCS_WITH_TOOLCHAIN
CCS_WITH_LIBSPE2
CCS_WITH_ADDRESSING
CCS_WITH_CELL_INSTDIR

################################################################################
# C/C++ configuration
################################################################################
AC_PROG_CC
CCS_TUNE_PCC
AC_PROG_CXX
CCS_TUNE_PCXX

################################################################################
# setup EXTRA
################################################################################
EXTRA_CPPFLAGS="$LIBSPE2_CPPFLAGS $TOOLCHAIN_CPPFLAGS"
EXTRA_LDFLAGS="$LIBSPE2_LDFLAGS"
EXTRA_LIBS="$LIBSPE2_LIBS"
AC_SUBST(EXTRA_CPPFLAGS)
AC_SUBST(EXTRA_LDFLAGS)
AC_SUBST(EXTRA_LIBS)

################################################################################
# enable/disable distributed memory parallelism
################################################################################
if test "$ENABLE_DMP" = "1" ; then
	CCS_WITH_PPU_MPI

	EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $PPU_MPI_CPPFLAGS"
	EXTRA_LDFLAGS="$EXTRA_LDFLAGS $PPU_MPI_LDFLAGS"
	EXTRA_LIBS="$EXTRA_LIBS $PPU_MPI_LIBS"

	AC_SUBST(EXTRA_CPPFLAGS)
	AC_SUBST(EXTRA_LDFLAGS)
	AC_SUBST(EXTRA_LIBS)

	# define needed for openmpi
	AC_DEFINE(OMPI_SKIP_MPICXX)
fi

################################################################################
# enable/disable hosted cell build
################################################################################
#if test "$ENABLE_HOST" = "1" ; then
#	AC_SUBST(PPEDIRS, ["bin"])
AC_SUBST(PPEDIRS, ["lib bin"])
#AC_SUBST(PPEDIRS, ["bin"])
#fi

################################################################################
# start output
################################################################################
AC_OUTPUT
