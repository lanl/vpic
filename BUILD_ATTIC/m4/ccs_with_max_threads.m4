dnl ----------------------------------------------------------------------------
dnl CCS_WITH_MAX_THREADS
dnl ----------------------------------------------------------------------------
AC_DEFUN([CCS_WITH_MAX_THREADS], [
    AC_ARG_VAR([CCS_MAX_THREADS], [Specify maximum number of threads])

    AC_ARG_WITH([max-threads],
        [AC_HELP_STRING([--with-max-threads],
        [Specify maximum number of threads])],
        [
            if test -n "$CCS_MAX_THREADS" ; then
                AC_SUBST(CCS_MAX_THREADS, $CCS_MAX_THREADS)
            elif test "$withval" != no ; then
                AC_SUBST(CCS_MAX_THREADS, $withval)
            fi
        ],
        [
            if test -n "$CCS_MAX_THREADS" ; then
                AC_SUBST(CCS_MAX_THREADS, $CCS_MAX_THREADS)
            else
                AC_SUBST(CCS_MAX_THREADS, 1)
            fi
    ])

    AC_MSG_CHECKING(maximum number of threads)
    AC_MSG_RESULT($CCS_MAX_THREADS)
])
