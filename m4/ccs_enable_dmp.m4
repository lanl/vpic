AC_DEFUN([CCS_ENABLE_DMP], [
    #
    # User hints...
    #
    AC_ARG_VAR([ENABLE_DMP], [Enable distributed memory parallelism])
    AC_ARG_ENABLE([dmp],
        [AC_HELP_STRING([--enable-dmp],
        [enable distributed memory parallelism])],
        [
            if test -n "$ENABLE_DMP" ; then
                enable_dmp=1
            elif test "$enableval" = "yes" ; then
                enable_dmp=1
            else
                enable_dmp=0
            fi

            if test "$enable_dmp" = "1" ; then
                AC_DEFINE(ENABLE_DMP, [1],
                    [define is you want IBM CBEA build enabled])
            fi

            AC_SUBST(ENABLE_DMP, [$enable_dmp])
        ],
        [
            if test -n "$ENABLE_DMP" ; then
                enable_dmp=1
            else
                enable_dmp=0
            fi

            if test "$enable_dmp" = "1" ; then
                AC_DEFINE(ENABLE_DMP, [1],
                    [define is you want distributed memory parallelism enabled])
            fi

            AC_SUBST(ENABLE_DMP, [$enable_dmp])
    ])
])
