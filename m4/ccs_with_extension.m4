AC_DEFUN([CCS_WITH_EXTENSION], [
    AC_ARG_VAR([EXTENSION], [User specified cpu type])

    AC_ARG_WITH([extension],
        [AC_HELP_STRING([--with-extension],
        [User specified cpu type, this affects some compiler options])],
        [
            if test -n "$EXTENSION" ; then
                AC_SUBST(EXTENSION, [$EXTENSION])
            elif test "$withval" != no ; then
                AC_SUBST(EXTENSION, [$withval])
            else
                AC_SUBST(EXTENSION, [generic])
            fi
        ],
        [
            if test -n "$EXTENSION" ; then
                AC_SUBST(EXTENSION, [$EXTENSION])
            else
                AC_SUBST(EXTENSION, [generic])
            fi
        ])
])
