AC_DEFUN([CCS_WITH_TOOLCHAIN], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([TOOLCHAIN],
        [User specified path to Toolchain library top-level directory])
    AC_ARG_VAR([TOOLCHAIN_CPPFLAGS],
        [User specified path to Toolchain library includes])

    AC_ARG_WITH([toolchain],
        [AC_HELP_STRING([--with-toolchain],
        [User specified path to Toolchain libraries])],
        [
            if test -n "$TOOLCHAIN" ; then
                if test -z "$TOOLCHAIN_CPPFLAGS" ; then
                    toolchain_CPPFLAGS="-I$TOOLCHAIN/include"
                else
                    toolchain_CPPFLAGS="$TOOLCHAIN_CPPFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$TOOLCHAIN_CPPFLAGS" ; then
                    toolchain_CPPFLAGS="-I$withval/include"
                else
                    toolchain_CPPFLAGS="$TOOLCHAIN_CPPFLAGS"
                fi
            fi
        ],
        [
            if test -n "$TOOLCHAIN" ; then
                if test -z "$TOOLCHAIN_CPPFLAGS" ; then
                    toolchain_CPPFLAGS="-I$TOOLCHAIN/include"
                else
                    toolchain_CPPFLAGS="$TOOLCHAIN_CPPFLAGS"
                fi
            else
                toolchain_CPPFLAGS="-I/opt/cell/toolchain/include"
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $toolchain_CPPFLAGS"

        AC_CHECK_HEADER(pthread.h, [toolchain_h=yes],
            [toolchain_h=no], /* check */)

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$toolchain_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid toolchain.h)
        else
            AC_SUBST(TOOLCHAIN_CPPFLAGS, [$toolchain_CPPFLAGS])
        fi
    ])
