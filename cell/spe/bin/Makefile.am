################################################################################
# SPE bin directory Makefile.am
################################################################################

clean-local:
	@(rm -rf `find @HOST_ACDIR@/src -name "*.timing"`)

bin_PROGRAMS = root_segment

INCLUDES = -I@top_srcdir@/src/include -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

root_segment_SOURCES = \
	@HOST_ACDIR@/src/util/pipelines/root_segment.c \
	@HOST_ACDIR@/src/species_advance/standard/pipelines_cell/advance_p_pipeline_spu.c \
	@HOST_ACDIR@/src/species_advance/standard/pipelines_cell/sort_p_pipeline_spu.c \
	@HOST_ACDIR@/src/sf_interface/pipelines_cell/clear_accumulators_pipeline_spu.c \
	@HOST_ACDIR@/src/sf_interface/pipelines_cell/reduce_accumulators_pipeline_spu.c \
	@HOST_ACDIR@/src/sf_interface/pipelines_cell/load_interpolator_pipeline_spu.c \
	@HOST_ACDIR@/src/sf_interface/pipelines_cell/unload_accumulator_pipeline_spu.c
# Use -Wl,--auto-overlay=spu_overlay.ld,--auto-relink to generate a skeleton
# overlay file.  It won't be good but it will work.
root_segment_LDFLAGS = @EXTRA_LDFLAGS@ -Wl,-T,spu_overlay.ld
#root_segment_LDFLAGS = @EXTRA_LDFLAGS@ -Wl,--auto-overlay=auto_overlay.ld,--auto-relink
root_segment_LDADD = @EXTRA_LIBS@
root_segment_DEPENDENCIES = @HOST_ACDIR@/cell/spe/bin/spu_overlay.ld

# FIXME: MAKE DISTCLEAN DOESN'T SEEM TO CLEAN THIS
if HAS_SPU_TIMING
root_segment_DEPENDENCIES += \
	@HOST_ACDIR@/src/species_advance/standard/pipelines_cell/advance_p_pipeline_spu.timing \
	@HOST_ACDIR@/src/species_advance/standard/pipelines_cell/sort_p_pipeline_spu.timing \
	@HOST_ACDIR@/src/sf_interface/pipelines_cell/clear_accumulators_pipeline_spu.timing \
	@HOST_ACDIR@/src/sf_interface/pipelines_cell/reduce_accumulators_pipeline_spu.timing \
	@HOST_ACDIR@/src/sf_interface/pipelines_cell/load_interpolator_pipeline_spu.timing \
	@HOST_ACDIR@/src/sf_interface/pipelines_cell/unload_accumulator_pipeline_spu.timing
include @top_builddir@/config/rules.mk
endif
