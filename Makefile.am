################################################################################
# top-level Makefile.am
################################################################################

SUBDIRS = ${CELLDIR} ${HOSTDIRS}

setup: setup-host copy-samples

all-local: copy-build-script copy-extension-lib

clean-local: clean-build-script clean-extension-lib

setup-host:
	$(call symlinks, "@abs_top_srcdir@/src", '.*\.h+?x+?')

copy-build-script:
	@(cp buildscript build.@EXTENSION@ ; chmod 755 build.@EXTENSION@)

copy-extension-lib:
	@(cp @ARCH_BUILD_PATH@/lib/.libs/libvpic.a @abs_top_builddir@/libvpic.@EXTENSION@.a)

copy-samples:
	@(if [ ! -d @abs_top_builddir@/sample ] ; then \
		cp -r @abs_top_srcdir@/sample @abs_top_builddir@ ; \
		rm -rf @abs_top_builddir@/sample/.svn ; \
	fi)

clean-build-script:
	@(rm -rf @abs_top_builddir@/build.@EXTENSION@)

clean-extension-lib:
	@(rm -rf @abs_top_builddir@/libvpic.@EXTENSION@.a)

clean-samples:
	@(rm -rf @abs_top_builddir@/sample)

################################################################################
# Function to make symbolic links
# arg1 - search directory
# arg2 - search wildcard
################################################################################
define symlinks
	@(cd $1; [ -d include ] && rm -rf include; \
		found=`find . -type f -regex $2 -print`; \
		[ -z "$$found" ] || (mkdir include && cd include && \
		ln -s `echo $$found | sed 's/\.\//\.\.\//g'` .))
endef

FORCE:
