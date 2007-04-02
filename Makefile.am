################################################################################
# top-level Makefile.am
################################################################################

SUBDIRS = ${CELLDIR} ${HOSTDIRS}

setup: setup-host

setup-host:
	$(call symlinks, "@top_srcdir@/src", '.*\.h+?x+?')

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
