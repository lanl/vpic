AC_DEFUN([CCS_WITH_OPT_LEVEL], [
    AC_ARG_VAR([OPT_LEVEL], [Specify optimization level <-2,-1,0,1,2,3>])

    AC_ARG_WITH([opt-level],
        [AC_HELP_STRING([--with-opt-level],
        [specify optimization level <-2,-1,0,1,2,3>])],
        [
            AC_MSG_CHECKING(optimization level)
            if test -n "$OPT_LEVEL" ; then
                AC_SUBST(OPT_LEVEL, [$OPT_LEVEL])
            elif test "$withval" != no ; then
                AC_SUBST(OPT_LEVEL, [$withval])
            else
                AC_SUBST(OPT_LEVEL, [0])
            fi

            AC_MSG_RESULT($OPT_LEVEL)
        ],
        [
            AC_MSG_CHECKING(optimization level)
            if test -n "$OPT_LEVEL" ; then
                AC_SUBST(OPT_LEVEL, [$OPT_LEVEL])
            else
                AC_SUBST(OPT_LEVEL, [0])
            fi

            AC_MSG_RESULT($OPT_LEVEL)
        ])
])
