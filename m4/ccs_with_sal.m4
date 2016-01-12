AC_DEFUN([CCS_WITH_SAL], [
    AC_ARG_VAR([SAL],
        [User specified path to SAL installation])
    AC_ARG_VAR([SAL_CPPFLAGS], [User specified path to SAL headers])

    AC_ARG_WITH([sal],
        [AC_HELP_STRING([--with-sal],
        [User specified path to SPE libraries])],
        [
            if test -n "$SAL" ; then
                if test -z "$SAL_CPPFLAGS" ; then
                    sal_CPPFLAGS="-I$SAL/include"
                else
                    sal_CPPFLAGS="$SAL_CPPFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$SAL_CPPFLAGS" ; then
                    sal_CPPFLAGS="-I$withval/include"
                else
                    sal_CPPFLAGS="$SAL_CPPFLAGS"
                fi
            fi
        ],
        [
            if test -n "$SAL" ; then
                if test -z "$SAL_CPPFLAGS" ; then
                    sal_CPPFLAGS="-I$SAL/include"
                else
                    sal_CPPFLAGS="$SAL_CPPFLAGS"
                fi
            else
                sal_CPPFLAGS="-I/opt/sal/include"
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $sal_CPPFLAGS"

        AC_CHECK_HEADER(sal.h,
            [sal_h=yes], [sal_h=no], /* check */)

        CPPFLAGS="$old_CPPFLAGS"

        if test "$sal_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid sal.h)
        else
            AC_SUBST(SAL_CPPFLAGS, [$sal_CPPFLAGS])
        fi
    ])
