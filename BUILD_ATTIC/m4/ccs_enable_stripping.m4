AC_DEFUN([CCS_ENABLE_STRIPPING], [
    #
    # User hints...
    #
    AC_ARG_VAR([STRIPPING], [Enable IBM CBEA build])
    AC_ARG_ENABLE([stripping],
        [AC_HELP_STRING([--enable-stripping],
        [enable IBM CBEA build])],
        [
            if test -n "$STRIPPING" ; then
                enable_stripping=1
            elif test "$enableval" = "yes" ; then
                enable_stripping=1
            else
                enable_stripping=0
            fi

            if test "$enable_stripping" = "1" ; then
                AC_DEFINE(STRIPPING, [1],
                    [define is you want IBM CBEA build enabled])
            fi

            AC_SUBST(ENABLE_STRIPPING, [$enable_stripping])
        ],
        [
            if test -n "$STRIPPING" ; then
                enable_stripping=1
            else
                enable_stripping=0
            fi

            if test "$enable_stripping" = "1" ; then
                AC_DEFINE(STRIPPING, [1],
                    [define is you want IBM CBEA build enabled])
            fi

            AC_SUBST(ENABLE_STRIPPING, [$enable_stripping])
    ])
])
