################################################################################
# top-level Makefile.am
################################################################################

SUBDIRS = ${CELLDIR} ${HOSTDIRS}

setup: setup-host

all-local: copy-build-script copy-extension-lib copy-samples

clean-local: clean-extension-lib clean-samples clean-build-script

setup-host:
	$(call symlinks, "@top_srcdir@/src", '.*\.h+?x+?')

copy-build-script:
	@(cp build.script build.@EXTENSION@ ; chmod 755 build.@EXTENSION@)

copy-extension-lib:
	@(cp lib/.libs/libvpic.a @top_builddir@/libvpic.@EXTENSION@.a)

copy-samples:
	@([ ! -d @top_builddir@/sample ] && \
		cp -r @top_srcdir@/sample @top_builddir@)

clean-build-script:
	rm -rf @top_builddir@/build.@EXTENSION@
#	@(rm -rf @top_builddir@/build.@EXTENSION@)

clean-extension-lib:
	rm -rf @top_builddir@/libvpic.@EXTENSION@.a
#	@(rm -rf @top_builddir@/libvpic.@EXTENSION@.a)

clean-samples:
	rm -rf @top_builddir@/sample
#	@(rm -rf @top_builddir@/sample)

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
