AC_DEFUN([CCS_PROG_DOXYGEN], [
    AC_ARG_VAR([DOXYGEN],
        [User specified path to the doxygen executable])

    AC_PATH_PROG(DOXYGEN, doxygen)

    if test -n "$DOXYGEN" ; then
        AC_SUBST(HAS_DOXYGEN, [yes])
    else
        AC_SUBST(HAS_DOXYGEN, [no])
    fi
])
