################################################################################
# bin Makefile.am
################################################################################

bin_PROGRAMS = advance

INCLUDES = @EXTRA_CPPFLAGS@ \
	-I@top_srcdir@/src/include

advance_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
advance_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/bench/advance_p
advance_LDFLAGS = @EXTRA_LDFLAGS@
advance_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@
