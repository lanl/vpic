AC_DEFUN([CCS_WITH_CELL_INSTDIR], [
    AC_ARG_VAR([CELL_INSTDIR], [Cell install target])

    AC_ARG_WITH([cell-instdir],
        [AC_HELP_STRING([--with-cell-instdir],
        [Cell install target])],
        [
            AC_MSG_CHECKING(cell install target)
            if test -n "$CELL_INSTDIR" ; then
                AC_SUBST(CELL_INSTDIR, [$CELL_INSTDIR])
                AC_MSG_RESULT([$CELL_INSTDIR])
            elif test "$withval" != no ; then
                AC_SUBST(CELL_INSTDIR, [$withval])
                AC_MSG_RESULT([$withval])
            else
                AC_MSG_RESULT([no])
            fi
        ],
        [
            AC_MSG_CHECKING(cell install target)
            if test -n "$CELL_INSTDIR" ; then
                AC_SUBST(CELL_INSTDIR, [$CELL_INSTDIR])
                AC_MSG_RESULT([$CELL_INSTDIR])
            else
                AC_MSG_RESULT([no])
            fi
    ])
])
