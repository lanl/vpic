################################################################################
# Server Makefile.am
################################################################################

bin_PROGRAMS = mprelay

INCLUDES = -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@
AM_CPPFLAGS = -DHOST_BUILD

mprelay_SOURCES = @HOST_ACDIR@/mprelay/mprelay.cxx
mprelay_LDFLAGS = @EXTRA_LDFLAGS@
mprelay_LDADD = @EXTRA_LIBS@
