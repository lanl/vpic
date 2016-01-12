dnl ----------------------------------------------------------------------------
dnl CCS_WITH_PTHREAD_AFFINITY
dnl ----------------------------------------------------------------------------
AC_DEFUN([CCS_WITH_PTHREAD_AFFINITY], [
    AC_ARG_VAR([CCS_PTHREAD_AFFINITY], [Specify CPU affinity policy for pthread creation <system, round_robin, cpu_zero>])

    AC_ARG_WITH([pthread-affinity],
        [AC_HELP_STRING([--with-pthread-affinity],
        [Specify CPU affinity policy for pthread creation])],
        [
            if test -n "$CCS_PTHREAD_AFFINITY" ; then
                AC_SUBST(CCS_PTHREAD_AFFINITY, $CCS_PTHREAD_AFFINITY)
            elif test "$withval" != no ; then
                AC_SUBST(CCS_PTHREAD_AFFINITY, $withval)
            fi
        ],
        [
            if test -n "$CCS_PTHREAD_AFFINITY" ; then
                AC_SUBST(CCS_PTHREAD_AFFINITY, $CCS_PTHREAD_AFFINITY)
            else
                AC_SUBST(CCS_PTHREAD_AFFINITY, "system")
            fi
    ])

    AC_MSG_CHECKING(pthread affinity policy)
    AC_MSG_RESULT($CCS_PTHREAD_AFFINITY)
])
