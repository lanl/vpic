################################################################################
# Server Makefile.am
################################################################################

bin_PROGRAMS = mprelay

INCLUDES = -I${top_srcdir}/src/include \
	@EXTRA_CPPFLAGS@
AM_CPPFLAGS = -DHOST_BUILD

mprelay_SOURCES = @abs_top_srcdir@/mprelay/mprelay.cxx
mprelay_LDFLAGS = @EXTRA_LDFLAGS@
mprelay_LDADD = @EXTRA_LIBS@
