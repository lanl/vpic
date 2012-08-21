dnl ----------------------------------------------------------------------------
dnl CCS_WITH_MEMORY_MODEL
dnl ----------------------------------------------------------------------------
AC_DEFUN([CCS_WITH_MEMORY_MODEL], [
    AC_ARG_VAR([MEMORY_MODEL], [])

    AC_MSG_CHECKING(memory model)

    AC_ARG_WITH([memory-model],
        [AC_HELP_STRING([--with-memory-model],
        [Specify memory model <local,spu>.])],
        [
            if test -n "$MEMORY_MODEL" ; then
                AC_SUBST(MEMORY_MODEL, $MEMORY_MODEL)
            elif test "$withval" != no ; then
                AC_SUBST(MEMORY_MODEL, $withval)
            fi
        ],
        [
            if test -n "$MEMORY_MODEL" ; then
                AC_SUBST(MEMORY_MODEL, $MEMORY_MODEL)
            else
                AC_SUBST(MEMORY_MODEL, "local")
            fi
    ])


    AC_MSG_RESULT($MEMORY_MODEL)
])
