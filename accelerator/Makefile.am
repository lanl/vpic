################################################################################
# Server Makefile.am
################################################################################

bin_PROGRAMS = accelerator

INCLUDES = -I${top_srcdir}/src/include \
	@EXTRA_CPPFLAGS@

accelerator_SOURCES = @abs_top_srcdir@/accelerator/accelerator.cxx
accelerator_LDFLAGS = @EXTRA_LDFLAGS@
accelerator_LDADD = @EXTRA_LIBS@
