################################################################################
# Server Makefile.am
################################################################################

bin_PROGRAMS = accelerator

INCLUDES = -I@HOST_ACDIR@/src/include \
	@EXTRA_CPPFLAGS@

accelerator_SOURCES = @HOST_ACDIR@/accelerator/accelerator.cxx
accelerator_LDFLAGS = @EXTRA_LDFLAGS@
accelerator_LDADD = @EXTRA_LIBS@
