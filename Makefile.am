################################################################################
# top-level Makefile.am
################################################################################

SUBDIRS = ${CELLDIR} ${HOSTDIRS} test ${DOC}

#setup: setup-host copy-samples
setup: copy-samples

all-local: copy-build-script copy-extension-lib @COPY_LOCAL@

clean-local: clean-build-script clean-extension-lib clean-test @CLEAN_LOCAL@

#setup-host:
#	$(call symlinks, "@abs_top_srcdir@/src", '.*\.h+?x+?')

copy-build-script:
	@(cp buildscript build.@LIB_EXTENSION@ ; chmod 755 build.@LIB_EXTENSION@)

copy-extension-lib:
	@(cp @ARCH_LIB_BUILD_PATH@/lib/.libs/libvpic.a \
	@abs_top_builddir@/libvpic.@LIB_EXTENSION@.a)

copy-accelerator-exe:
	@(cp @ARCH_ACCEL_BUILD_PATH@/accelerator/accelerator \
	@abs_top_builddir@/accelerator.@ACCEL_EXTENSION@)

copy-mprelay-exe:
	@(cp @ARCH_MPRELAY_BUILD_PATH@/mprelay/mprelay \
	@abs_top_builddir@/mprelay.@MPRELAY_EXTENSION@)

copy-samples:
	@(if [ ! -d @abs_top_builddir@/sample ] ; then \
		cp -r @abs_top_srcdir@/sample @abs_top_builddir@ ; \
		rm -rf @abs_top_builddir@/sample/.svn ; \
	fi)

clean-build-script:
	@(rm -rf @abs_top_builddir@/build.@LIB_EXTENSION@)

clean-extension-lib:
	@(rm -rf @abs_top_builddir@/libvpic.@LIB_EXTENSION@.a)

clean-accelerator-exe:
	@(rm -rf @abs_top_builddir@/accelerator.@ACCEL_EXTENSION@)

clean-mprelay-exe:
	@(rm -rf @abs_top_builddir@/mprelay.@MPRELAY_EXTENSION@)

clean-samples:
	@(rm -rf @abs_top_builddir@/sample)

clean-test:
	@(cd test; ${MAKE} clean)

doc: FORCE
	@([ -d @DOC@ ] && cd @DOC@ && make doc > doc.log 2>&1)

#include config/utils.mk

FORCE:
