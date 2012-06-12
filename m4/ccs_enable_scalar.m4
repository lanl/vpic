AC_DEFUN([CCS_ENABLE_SCALAR], [
    AC_ARG_VAR([SCALAR], [Enable scalar mode])

    AC_ARG_ENABLE([scalar],
        [AC_HELP_STRING([--enable-scalar],
        [enable scalar mode])],
        [
            if test -n "$SCALAR" ; then
                enable_scalar=1
            elif test "$enableval" = "yes" ; then
                enable_scalar=1
            else
                enable_scalar=0
            fi
        ],
        [
            if test -n "$SCALAR" ; then
                enable_scalar=1
            else
                enable_scalar=0
            fi
        ])

        if test "$enable_scalar" = "1" ; then
            AC_DEFINE(SCALAR, [1],
                [define if you want SCALAR mode])
            AC_SUBST(ENABLE_SCALAR, [1])
        else
            AC_SUBST(ENABLE_SIMD, [1])
        fi
])
