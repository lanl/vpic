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
	mprelay/Makefile \
	accelerator/Makefile \
	../../buildscript)

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
AC_SUBST(HOST_BUILDDIR, [$host_build_directory])
AC_SUBST(HOST_ACDIR, [$host_ac_directory])
echo "host directory $host_ac_directory"

################################################################################
# define to get PPU code
################################################################################
#AC_DEFINE(CELL_PPU_BUILD)

################################################################################
# machine file specification
################################################################################
CCS_WITH_MACHINE("$host_ac_directory/machine")

################################################################################
# with options
################################################################################
CCS_WITH_BUILDSTYLE
CCS_WITH_PPU_COMPILER
CCS_WITH_OPT_LEVEL
CCS_WITH_LIBSPE2
CCS_WITH_ADDRESSING
CCS_WITH_CELL_INSTDIR
CCS_WITH_EXTENSION
CCS_WITH_SPU_MODULES

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

case "$buildstyle" in

	standard_cell)
		CCS_WITH_PPU_MPI

		AC_SUBST(PPEDIRS, ["lib"])

		EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $PPU_MPI_CPPFLAGS"
		EXTRA_LDFLAGS="$EXTRA_LDFLAGS $PPU_MPI_LDFLAGS"
		EXTRA_LIBS="$EXTRA_LIBS $PPU_MPI_LIBS"

		AC_DEFINE(CELL_PPU_BUILD)
		AC_DEFINE(USE_CELL_SPUS)

		# define needed for openmpi
		AC_DEFINE(OMPI_SKIP_MPICXX)

		for file in "$SPU_MODULES" ; do
			EXTRA_O_FILES=${EXTRA_O_FILES}`echo emb-$file.o`
			EXTRA_LO_FILES=${EXTRA_LO_FILES}`echo emb-$file.lo`
		done
	;;

	ompi_relay_cell)
		CCS_WITH_PPU_MPI

		AC_SUBST(PPEDIRS, ["lib mprelay accelerator"])

		EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $PPU_MPI_CPPFLAGS"
		EXTRA_LDFLAGS="$EXTRA_LDFLAGS $PPU_MPI_LDFLAGS"
		EXTRA_LIBS="$EXTRA_LIBS $PPU_MPI_LIBS"

		AC_DEFINE(USE_MPRELAY)
		AC_DEFINE(USE_MPI_P2P)

		AC_DEFINE(CELL_PPU_BUILD)
		AC_DEFINE(USE_CELL_SPUS)

		# define needed for openmpi
		AC_DEFINE(OMPI_SKIP_MPICXX)

		for file in "$SPU_MODULES" ; do
			EXTRA_O_FILES=${EXTRA_O_FILES}`echo emb-$file.o`
			EXTRA_LO_FILES=${EXTRA_LO_FILES}`echo emb-$file.lo`
		done
	;;

	ompi_relay_hybrid_cell)
		CCS_WITH_PPU_MPI

		AC_SUBST(PPEDIRS, ["lib accelerator"])

		EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $PPU_MPI_CPPFLAGS"
		EXTRA_LDFLAGS="$EXTRA_LDFLAGS $PPU_MPI_LDFLAGS"
		EXTRA_LIBS="$EXTRA_LIBS $PPU_MPI_LIBS"

		AC_DEFINE(USE_MPRELAY)
		AC_DEFINE(USE_MPI_P2P)

		AC_DEFINE(CELL_PPU_BUILD)
		AC_DEFINE(USE_CELL_SPUS)

		# define needed for openmpi
		AC_DEFINE(OMPI_SKIP_MPICXX)

		for file in "$SPU_MODULES" ; do
			EXTRA_O_FILES=${EXTRA_O_FILES}`echo emb-$file.o`
			EXTRA_LO_FILES=${EXTRA_LO_FILES}`echo emb-$file.lo`
		done
	;;

	ompi_relay_hybrid_cell_ppe)
		CCS_WITH_PPU_MPI

		AC_SUBST(PPEDIRS, ["lib accelerator"])

		EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $PPU_MPI_CPPFLAGS"
		EXTRA_LDFLAGS="$EXTRA_LDFLAGS $PPU_MPI_LDFLAGS"
		EXTRA_LIBS="$EXTRA_LIBS $PPU_MPI_LIBS"

		AC_DEFINE(USE_MPRELAY)
		AC_DEFINE(USE_MPI_P2P)

		AC_DEFINE(CELL_PPU_BUILD)

		# define needed for openmpi
		AC_DEFINE(OMPI_SKIP_MPICXX)
	;;

	dacs_relay_hybrid_cell)
		CCS_WITH_PPU_DACS(hybrid)

		AC_SUBST(PPEDIRS, ["lib accelerator"])

		EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $PPU_DACS_CPPFLAGS"
		EXTRA_LDFLAGS="$EXTRA_LDFLAGS $PPU_DACS_LDFLAGS"
		EXTRA_LIBS="$EXTRA_LIBS $PPU_DACS_LIBS"

		AC_DEFINE(USE_MPRELAY)
		AC_DEFINE(USE_DACS_P2P)

		AC_DEFINE(CELL_PPU_BUILD)
		AC_DEFINE(USE_SPU_PIPELINE)
	;;

esac

AC_SUBST(EXTRA_CPPFLAGS)
AC_SUBST(EXTRA_LDFLAGS)
AC_SUBST(EXTRA_LIBS)

AC_SUBST(EXTRA_O_FILES)
AC_SUBST(EXTRA_LO_FILES)

################################################################################
# start output
################################################################################
AC_OUTPUT
