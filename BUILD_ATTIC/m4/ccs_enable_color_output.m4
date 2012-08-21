AC_DEFUN([CCS_ENABLE_COLOR_OUTPUT], [
    #
    # User hints...
    #
    AC_ARG_VAR([HAVE_COLOR_OUTPUT], [Enable colorized output])
    AC_ARG_ENABLE([color-output],
        [AC_HELP_STRING([--enable-color-output],
        [enable colorized output, this is useful for debugging if your shell supports color output])],
        [
            if test -n "$HAVE_COLOR_OUTPUT" ; then
                enable_color_output=1
            elif test "$enableval" = "yes" ; then
                enable_color_output=1
            else
                enable_color_output=0
            fi

            if test "$enable_color_output" = "1" ; then
                AC_DEFINE(HAVE_COLOR_OUTPUT, [1],
                    [define if you want colorized output])
            fi

            AC_SUBST(ENABLE_COLOR_OUTPUT, [$enable_color_output])
        ],
        [
            if test -n "$HAVE_COLOR_OUTPUT" ; then
                enable_color_output=1
            else
                enable_color_output=0
            fi

            if test "$enable_color_output" = "1" ; then
                AC_DEFINE(HAVE_COLOR_OUTPUT, [1],
                    [define if you want colorized output])
            fi

            AC_SUBST(ENABLE_COLOR_OUTPUT, [$enable_color_output])
    ])
])
