AC_DEFUN([CCS_WITH_IBM_SIMDMATH], [
    AC_ARG_VAR([IBM_SIMDMATH],
        [User specified path to SPE SIMD math library top-level directory])
    AC_ARG_VAR([IBM_SIMDMATH_CPPFLAGS],
        [User specified path to SPE SIMD math library includes])
    AC_ARG_VAR([IBM_SIMDMATH_LDFLAGS],
        [User specified path to SPE SIMD math libraries])
    AC_ARG_VAR([IBM_SIMDMATH_LIBS],
        [User specified link library names, e.g., -lsimdmath])

    AC_ARG_WITH([ibm-simdmath],
        [AC_HELP_STRING([--with-ibm-simdmath],
        [User specified path to SPE SIMD math libraries])],
        [
            if test -n "$IBM_SIMDMATH" ; then
                if test -z "$IBM_SIMDMATH_CPPFLAGS" ; then
                    ibm_simdmath_CPPFLAGS="-I$IBM_SIMDMATH/include"
                else
                    ibm_simdmath_CPPFLAGS="$IBM_SIMDMATH_CPPFLAGS"
                fi

                if test -z "$IBM_SIMDMATH_LDFLAGS" ; then
                    ibm_simdmath_LDFLAGS="-L$IBM_SIMDMATH/lib"
                else
                    ibm_simdmath_LDFLAGS="$IBM_SIMDMATH_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$IBM_SIMDMATH_CPPFLAGS" ; then
                    ibm_simdmath_CPPFLAGS="-I$withval/include"
                else
                    ibm_simdmath_CPPFLAGS="$IBM_SIMDMATH_CPPFLAGS"
                fi

                if test -z "$IBM_SIMDMATH_LDFLAGS" ; then
                    ibm_simdmath_LDFLAGS="-L$withval/lib"
                else
                    ibm_simdmath_LDFLAGS="$IBM_SIMDMATH_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$IBM_SIMDMATH" ; then
                if test -z "$IBM_SIMDMATH_CPPFLAGS" ; then
                    ibm_simdmath_CPPFLAGS="-I$IBM_SIMDMATH/include"
                else
                    ibm_simdmath_CPPFLAGS="$IBM_SIMDMATH_CPPFLAGS"
                fi

                if test -z "$IBM_SIMDMATH_LDFLAGS" ; then
                    ibm_simdmath_LDFLAGS="-L$IBM_SIMDMATH/lib"
                else
                    ibm_simdmath_LDFLAGS="$IBM_SIMDMATH_LDFLAGS"
                fi
            else
                DEFAULT_PATH="/opt/cell/sysroot/usr/spu"
                ibm_simdmath_CPPFLAGS="-I$DEFAULT_PATH/include"
                ibm_simdmath_LDFLAGS="-L$DEFAULT_PATH/lib"
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $ibm_simdmath_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $ibm_simdmath_LDFLAGS"

        AC_CHECK_HEADER(simdmath.h,
            [ibm_simdmath_h=yes], [ibm_simdmath_h=no], /* check */)
        AC_CHECK_LIB(simdmath, divf4, [ibm_simdmath_a=yes], [ibm_simdmath_a=no])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$ibm_simdmath_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid simdmath.h)
        else
            AC_SUBST(IBM_SIMDMATH_CPPFLAGS, [$ibm_simdmath_CPPFLAGS])
        fi

        if test "$ibm_simdmath_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libsimdmath.a)
        else
            AC_SUBST(IBM_SIMDMATH_LDFLAGS, [$ibm_simdmath_LDFLAGS])
        fi

        if test -z "$IBM_SIMDMATH_LIBS" ; then
            AC_SUBST(IBM_SIMDMATH_LIBS, [-lsimdmath])
        fi
    ])
