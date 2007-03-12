################################################################################
# SPE bin directory Makefile.am
################################################################################

bin_PROGRAMS = test

INCLUDES = -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

test_SOURCES = @HOST_ACDIR@/src/cell/spe/simple/test.cc
test_LDFLAGS = @EXTRA_LDFLAGS@
test_LDADD = @EXTRA_LIBS@
