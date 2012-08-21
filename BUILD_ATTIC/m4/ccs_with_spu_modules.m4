dnl ----------------------------------------------------------------------------
dnl CCS_WITH_SPU_MODULES
dnl ----------------------------------------------------------------------------
AC_DEFUN([CCS_WITH_SPU_MODULES], [
    AC_ARG_VAR([SPU_MODULES], [SPU modules to link])

    AC_ARG_WITH([spe-modules],
        [AC_HELP_STRING([--with-spe-modules],
        [Specify SPU modules to link])],
        [
            if test -n "$SPU_MODULES" ; then
                AC_SUBST(SPU_MODULES, $SPU_MODULES)
            elif test "$withval" != no ; then
                AC_SUBST(SPU_MODULES, $withval)
            fi
        ],
        [
            if test -n "$SPU_MODULES" ; then
                AC_SUBST(SPU_MODULES, $SPU_MODULES)
            fi
    ])
])
