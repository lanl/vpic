#------------------------------------------------------------------------------#
# Top-level Makefile.am
#------------------------------------------------------------------------------#

SUBDIRS = lib bin test

#------------------------------------------------------------------------------#
# Create links
#------------------------------------------------------------------------------#

links:
	@(cd ${top_srcdir}/src; [ -d include ] && rm -rf include; \
		mkdir include && cd include; \
		for file in `find .. -type f -regex '.*\.h*x*'` ; do \
			link=`echo $$file`; \
			ln -s $$link .; \
		done \
	)

#------------------------------------------------------------------------------#
# Extra install data
#------------------------------------------------------------------------------#

install-data-hook:
	@(([ -d ${prefix}/src ] || mkdir -p ${prefix}/src) && \
		([ -d ${prefix}/sample ] || mkdir -p ${prefix}/sample) && \
		cp ${top_srcdir}/deck/main.cxx ${prefix}/src && \
		cp ${top_srcdir}/deck/wrapper.cxx ${prefix}/src && \
		cp ${top_srcdir}/sample/harris ${prefix}/sample)

#------------------------------------------------------------------------------#
# Extra dist
#------------------------------------------------------------------------------#

dist-hook:
	@(mkdir -p ${distdir}/machine && \
		cp -r ${top_srcdir}/machine/* ${distdir}/machine)
	@(mkdir -p ${distdir}/test/integrated && \
		cp ${top_srcdir}/test/integrated/*.deck ${distdir}/test/integrated)
	@(mkdir -p ${distdir}/sample && \
		cp -r ${top_srcdir}/sample/* ${distdir}/sample)			
	@(mkdir -p ${distdir}/config && \
		cp ${top_srcdir}/config/local.at ${distdir}/config)
