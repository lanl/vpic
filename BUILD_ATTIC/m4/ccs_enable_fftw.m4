AC_DEFUN([CCS_ENABLE_FFTW], [
    #
    # User hints...
    #
    AC_ARG_VAR([ENABLE_FFTW], [Enable fftw linkage])
    AC_ARG_ENABLE([fftw],
        [AC_HELP_STRING([--enable-fftw],
        [enable OpenSSL linkage])],
        [
            if test -n "$ENABLE_FFTW" ; then
                enable_fftw=1
            elif test "$enableval" = "yes" ; then
                enable_fftw=1
            else
                enable_fftw=0
            fi
        ],
        [
            if test -n "$ENABLE_FFTW" ; then
                enable_fftw=1
            else
                enable_fftw=0
            fi
        ])

        AC_SUBST(ENABLE_FFTW, [$enable_fftw])

		if test "$enable_fftw" = "1" ; then
			AC_DEFINE(ENABLE_FFTW, [1],
				[define if you want to use fftw])
		fi
])
