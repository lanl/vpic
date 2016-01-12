AC_DEFUN([CCS_ENABLE_INLINE], [
    #
    # User hints...
    #
    AC_ARG_VAR([ENABLE_INLINE], [Force or disable inlining.])
    AC_ARG_ENABLE([inline],
        [AC_HELP_STRING([--enable-inline], [Force or disable inlining.
            By default, forced inlining is enabled using
            __attribute__((always_inline))])],
        [
            if test -n "$ENABLE_INLINE" ; then
                AC_DEFINE(SAL_INLINE)
            elif test "$enableval" = "yes" ; then
                AC_DEFINE(SAL_INLINE)
            else
                AC_DEFINE(SAL_NOINLINE)
            fi
        ],
        [
            if test -n "$ENABLE_INLINE" ; then
                AC_DEFINE(SAL_INLINE)
            else
                AC_DEFINE(SAL_NOINLINE)
            fi
        ])
])
