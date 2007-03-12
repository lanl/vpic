################################################################################
# bin Makefile.am
################################################################################

bin_PROGRAMS = advance

INCLUDES = -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

advance_SOURCES = @HOST_ACDIR@/src/main.cpp \
	@HOST_ACDIR@/src/deck_wrapper.cpp
advance_CPPFLAGS = -DINPUT_DECK=@HOST_ACDIR@/sample/bench/advance_p
advance_LDFLAGS = @EXTRA_LDFLAGS@
advance_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

