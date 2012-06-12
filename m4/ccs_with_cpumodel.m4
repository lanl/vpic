AC_DEFUN([CCS_WITH_CPUMODEL], [
    AC_ARG_VAR([CPUMODEL], [User specified cpu type])

    AC_ARG_WITH([cpumodel],
        [AC_HELP_STRING([--with-cpumodel],
        [User specified cpu type, this affects some compiler options])],
        [
            if test -n "$CPUMODEL" ; then
                AC_SUBST(CPUMODEL, [$CPUMODEL])
            elif test "$withval" != no ; then
                AC_SUBST(CPUMODEL, [$withval])
            else
                AC_SUBST(CPUMODEL, [generic])
            fi
        ],
        [
            if test -n "$CPUMODEL" ; then
                AC_SUBST(CPUMODEL, [$CPUMODEL])
            else
                AC_SUBST(CPUMODEL, [generic])
            fi
        ])
])
