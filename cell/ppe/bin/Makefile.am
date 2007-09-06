################################################################################
# bin Makefile.am
################################################################################

bin_PROGRAMS = advance_p

INCLUDES = -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

advance_p_SOURCES = @HOST_ACDIR@/tmp/advance_p/advance_p.c
advance_p_LDFLAGS = @EXTRA_LDFLAGS@
advance_p_DEPENDENCIES = @EXTRA_O_FILES@
advance_p_LDADD = @EXTRA_LO_FILES@ ../lib/libvpic.la \
	@EXTRA_LIBS@

include @top_builddir@/config/rules.mk
