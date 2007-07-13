################################################################################
# SPE bin directory Makefile.am
################################################################################

bin_PROGRAMS = advance_p_pipeline_spu

INCLUDES = -I@top_srcdir@/src/include -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

advance_p_pipeline_spu_SOURCES = \
	@top_srcdir@/src/particle/advance_p_pipeline_spu.c
advance_p_pipeline_spu_LDFLAGS = @EXTRA_LDFLAGS@
advance_p_pipeline_spu_LDADD = @EXTRA_LIBS@
