AC_DEFUN([CCS_WITH_FFTW], [

    AC_ARG_VAR([FFTW],
        [User specified path to FFTW library top-level directory])
    AC_ARG_VAR([FFTW_CPPFLAGS], [User specified path to FFTW library includes])
    AC_ARG_VAR([FFTW_LDFLAGS], [User specified path to FFTW libraries])
    AC_ARG_VAR([FFTW_LIBS],
        [User specified link library names, e.g., -lfftw3])

    AC_ARG_WITH([fftw],
        [AC_HELP_STRING([--with-fftw],
        [User specified path to SPE libraries])],
        [
            if test -n "$FFTW" ; then
                if test -z "$FFTW_CPPFLAGS" ; then
                    fftw_CPPFLAGS="-I$FFTW/include"
                else
                    fftw_CPPFLAGS="$FFTW_CPPFLAGS"
                fi

                if test -z "$FFTW_LDFLAGS" ; then
                    fftw_LDFLAGS="-L$FFTW/lib"
                else
                    fftw_LDFLAGS="$FFTW_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$FFTW_CPPFLAGS" ; then
                    fftw_CPPFLAGS="-I$withval/include"
                else
                    fftw_CPPFLAGS="$FFTW_CPPFLAGS"
                fi

                if test -z "$FFTW_LDFLAGS" ; then
                    fftw_LDFLAGS="-L$withval/lib"
                else
                    fftw_LDFLAGS="$FFTW_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$FFTW" ; then
                if test -z "$FFTW_CPPFLAGS" ; then
                    fftw_CPPFLAGS="-I$FFTW/include"
                else
                    fftw_CPPFLAGS="$FFTW_CPPFLAGS"
                fi

                if test -z "$FFTW_LDFLAGS" ; then
                    fftw_LDFLAGS="-L$FFTW/lib"
                else
                    fftw_LDFLAGS="$FFTW_LDFLAGS"
                fi
            else
                if test -z "$FFTW_CPPFLAGS" ; then
                    fftw_CPPFLAGS="-I/usr/include"
                else
                    fftw_CPPFLAGS="$FFTW_CPPFLAGS"
                fi

                if test -z "$FFTW_LDFLAGS" ; then
                    fftw_LDFLAGS="-L/usr/lib"
                else
                    fftw_LDFLAGS="$FFTW_LDFLAGS"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $fftw_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $fftw_LDFLAGS"

        AC_CHECK_HEADER(fftw3.h, [fftw_h=yes], [fftw_h=no], /* check */)
        AC_CHECK_LIB(fftw3, fftw_transpose, [fftw_a=yes], [fftw_a=no], [])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$fftw_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid fftw3.h)
        else
            AC_SUBST(FFTW_CPPFLAGS, [$fftw_CPPFLAGS])
        fi

        if test "$fftw_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libfftw3.a)
        else
            AC_SUBST(FFTW_LDFLAGS, [$fftw_LDFLAGS])
        fi

        if test -z "$FFTW_LIBS" ; then
            AC_SUBST(FFTW_LIBS, ["-lfftw3"])
        fi
    ])
