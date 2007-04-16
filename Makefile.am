################################################################################
# top-level Makefile.am
################################################################################

SUBDIRS = ${CELLDIR} ${HOSTDIRS}

all-local: copy-extension-lib

setup: setup-host setup-build-script

setup-host:
	$(call symlinks, "@top_srcdir@/src", '.*\.h+?x+?')

setup-build-script:
	@(mv build.script build.@EXTENSION@ ; chmod 755 build.@EXTENSION@)

copy-extension-lib:
	@(cp lib/.libs/libvpic.a @top_srcdir@/libvpic.@EXTENSION@.a)

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
