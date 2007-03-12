################################################################################
# bin Makefile.am
################################################################################

bin_PROGRAMS = example

INCLUDES = -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

example_SOURCES = @HOST_ACDIR@/src/main.cpp \
	@HOST_ACDIR@/src/deck_wrapper.cpp
example_CPPFLAGS = -DINPUT_DECK=@HOST_ACDIR@/sample/bench/advance_p
example_LDFLAGS = @EXTRA_LDFLAGS@
example_LDADD = libtest.o @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

include @top_builddir@/config/rules.mk
