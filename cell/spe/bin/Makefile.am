################################################################################
# SPE bin directory Makefile.am
################################################################################

clean-local:
	@(rm -rf `find @HOST_ACDIR@/src -name "*.timing"`)

bin_PROGRAMS = root_segment

INCLUDES = -I@top_srcdir@/src/include -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

SPECIES_SPE_DIR = species_advance/standard/pipelines_cell

root_segment_SOURCES = \
	@HOST_ACDIR@/src/util/overlays/root_segment.c \
	@HOST_ACDIR@/src/${SPECIES_SPE_DIR}/advance_p_overlay.c \
	@HOST_ACDIR@/src/${SPECIES_SPE_DIR}/sort_p_overlay.c
root_segment_LDFLAGS = @EXTRA_LDFLAGS@ -Wl,-T,@HOST_ACDIR@/src/util/overlays/root_segment.ld
root_segment_LDADD = @EXTRA_LIBS@

if HAS_SPU_TIMING
root_segment_DEPENDENCIES = \
	@HOST_ACDIR@/src/species_advance/standard/pipelines_cell/advance_p_overlay.timing
include @top_builddir@/config/rules.mk
endif
