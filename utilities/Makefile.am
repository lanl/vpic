################################################################################
# Server Makefile.am
################################################################################

bin_PROGRAMS = restart_remap

INCLUDES = -I${top_srcdir}/src/include -I${top_srcdir}/utilities \
	@EXTRA_CPPFLAGS@

restart_remap_SOURCES = @abs_top_srcdir@/utilities/restart_remap.cxx
restart_remap_LDFLAGS = @EXTRA_LDFLAGS@
restart_remap_LDADD = @EXTRA_LIBS@
