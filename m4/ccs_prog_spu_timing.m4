AC_DEFUN([CCS_PROG_SPU_TIMING], [
    AC_ARG_VAR([SPU_TIMING],
        [User specified path to the doxygen executable])

    AC_PATH_PROG(SPU_TIMING, spu_timing)

    if test -n "$SPU_TIMING" ; then
        AC_SUBST(HAS_SPU_TIMING, [yes])
    else
        AC_SUBST(HAS_SPU_TIMING, [no])
    fi
])
