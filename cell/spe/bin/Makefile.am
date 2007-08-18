################################################################################
# SPE bin directory Makefile.am
################################################################################

clean-local:
	@(rm -rf *.timing *.asm; rm -rf `find @HOST_ACDIR@/src -name "*.timing"`)

bin_PROGRAMS = advance_p_pipeline_spu

INCLUDES = -I@top_srcdir@/src/include -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

advance_p_pipeline_spu_SOURCES = \
	@HOST_ACDIR@/src/particle/pipelines_cell/advance_p_pipeline_spu.c
advance_p_pipeline_spu_LDFLAGS = @EXTRA_LDFLAGS@
advance_p_pipeline_spu_LDADD = @EXTRA_LIBS@

if HAS_SPU_TIMING
advance_p_pipeline_spu_DEPENDENCIES = \
	@HOST_ACDIR@/src/particle/pipelines_cell/advance_p_pipeline_spu.timing

include @top_builddir@/config/rules.mk
endif
