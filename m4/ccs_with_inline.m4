dnl ----------------------------------------------------------------------------
dnl CCS_WITH_INLINE
dnl ----------------------------------------------------------------------------
AC_DEFUN([CCS_WITH_INLINE], [
    AC_ARG_VAR([CCS_INLINE], [Specify function inlinging strategy <always,never,compiler>])

    AC_ARG_WITH([inline],
        [AC_HELP_STRING([--with-inline],
        [Specify function inlinging strategy <always,never,compiler>])],
        [
            if test -n "$CCS_INLINE" ; then
                AC_SUBST(CCS_INLINE, $CCS_INLINE)
            elif test "$withval" != no ; then
                AC_SUBST(CCS_INLINE, $withval)
            fi
        ],
        [
            if test -n "$CCS_INLINE" ; then
                AC_SUBST(CCS_INLINE, $CCS_INLINE)
            else
                AC_SUBST(CCS_INLINE, "compiler")
            fi
    ])

    AC_MSG_CHECKING(inlining strategy)
    AC_MSG_RESULT($CCS_INLINE)
])
