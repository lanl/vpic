AC_DEFUN([CCS_WITH_FPTYPE], [
    AC_ARG_VAR([FPTYPE], [Floating Point Type e.g. double, float, ...])

    AC_ARG_WITH([fptype],
        [AC_HELP_STRING([--with-fptype],
        [floating point type e.g. double, float, ...])],
        [
            if test -n "$FPTYPE" ; then
                AC_DEFINE(FPTYPE, $FPTYPE)
            elif test "$withval" != no ; then
                AC_DEFINE(FPTYPE, $withval)
            fi
        ],
        [
            if test -n "$FPTYPE" ; then
                AC_DEFINE(FPTYPE, $FPTYPE)
            else
                AC_DEFINE(FPTYPE, double)
            fi
    ])
])
