################################################################################
# top-level Makefile.am
################################################################################

SUBDIRS = ${CELLDIR} ${HOSTDIR}

setup:
	$(call symlinks, '.*\.h+?p+?')

################################################################################
# Function to make symbolic links
# arg1 - search wildcard
################################################################################
define symlinks
	@(cd src; [ -d include ] && rm -rf include; \
		found=`find . -type f -regex $1 -print`; \
		[ -z "$$found" ] || (mkdir include && cd include && \
		ln -s `echo $$found | sed 's/\.\//\.\.\//g'` .))
endef

FORCE:
