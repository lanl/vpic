AC_DEFUN([CCS_ENABLE_PROFILING], [
    #
    # User hints...
    #
    AC_ARG_VAR([PROFILING], [Enable Profiling])
    AC_ARG_ENABLE([profiling],
        [AC_HELP_STRING([--enable-profiling],
        [enable profiling [ARG=profiler]])],
        [
            if test -n "$PROFILING" ; then
                enable_profiling=$PROFILING
            elif test "$enableval" = "yes" ; then
                enable_profiling=native
            elif test -n "$enableval" -a "$enableval" != "no"; then
                enable_profiling=$enableval
            else
                AC_MSG_ERROR([You must specify which profiler to use])
            fi

            if test -n "$enable_profiling" ; then
                AC_DEFINE(PROFILING, [1],
                    [define if you want to enable profiling])
            fi

            AC_SUBST(ENABLE_PROFILING, [$enable_profiling])
        ],
        [
            if test -n "$PROFILING" ; then
                enable_profiling=$PROFILING
            else
                enable_profiling=0
            fi

            if test -n "$enable_profiling" ; then
                AC_DEFINE(PROFILING, [1],
                    [define if you want distributed memory parallelism])
            fi

            AC_SUBST(ENABLE_PROFILING, [$enable_profiling])
    ])
])
