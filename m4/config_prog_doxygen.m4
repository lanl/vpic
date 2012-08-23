dnl ------------------------------------------------------------------------ dnl
dnl Copyright (c) 2012 Los Alamos National Security, LLC
dnl All rights reserved.
dnl ------------------------------------------------------------------------ dnl

dnl ------------------------------------------------------------------------ dnl
dnl Detect location of the doxygen program
dnl ------------------------------------------------------------------------ dnl

AC_DEFUN([CONFIG_PROG_DOXYGEN], [
    AC_ARG_VAR([DOXYGEN],
        [User specified path to the doxygen executable])

    AC_PATH_PROG(DOXYGEN, doxygen)

    if test -n "$DOXYGEN" ; then
        AC_SUBST(HAS_DOXYGEN, [yes])
    else
        AC_SUBST(HAS_DOXYGEN, [no])
    fi
])
