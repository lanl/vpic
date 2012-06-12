AC_DEFUN([CCS_ENABLE_HOST], [
    #
    # User hints...
    #
    AC_ARG_VAR([ENABLE_HOST], [Enable host build])
    AC_ARG_ENABLE([host],
        [AC_HELP_STRING([--enable-host],
        [enable IBM CBEA build])],
        [
            if test -n "$ENABLE_HOST" ; then
                enable_host=1
            elif test "$enableval" = "yes" ; then
                enable_host=1
            else
                enable_host=0
            fi
        ],
        [
            if test -n "$ENABLE_HOST" ; then
                enable_host=1
            else
                enable_host=1
            fi
        ])

        if test "$enable_host" = "1" ; then
            AC_DEFINE(ENABLE_HOST, [1],
                [define is you want IBM CBEA build enabled])
        fi

        AC_SUBST(ENABLE_HOST, [$enable_host])
])
