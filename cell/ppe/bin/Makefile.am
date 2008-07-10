################################################################################
# bin Makefile.am
################################################################################

#bin_PROGRAMS = advance_p simple
bin_PROGRAMS = simple

INCLUDES = -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

#advance_p_SOURCES = @HOST_ACDIR@/tmp/advance_p/advance_p.c
#advance_p_LDFLAGS = @EXTRA_LDFLAGS@
#advance_p_LDADD = lib-advance_p_pipeline_spu.o ../lib/libvpic.la @EXTRA_LIBS@

simple_SOURCES = @top_srcdir@/src/simple/simple.c
simple_LDFLAGS = @EXTRA_LDFLAGS@
simple_LDADD = lib-simple.o @EXTRA_LIBS@

include @top_builddir@/config/rules.mk
