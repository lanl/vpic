AC_DEFUN([CCS_ENABLE_UNIT_TESTS], [
    #
    # User hints...
    #
    AC_ARG_VAR([UNIT_TESTS], [Enable distributed memory parallelism])
    AC_ARG_ENABLE([unit-tests],
        [AC_HELP_STRING([--enable-unit-tests],
        [enable distributed memory parallelism])],
        [
            if test -n "$UNIT_TESTS" ; then
                enable_unit_tests=1
            elif test "$enableval" = "yes" ; then
                enable_unit_tests=1
            else
                enable_unit_tests=0
            fi

            if test "$enable_unit_tests" = "1" ; then
                AC_DEFINE(UNIT_TESTS, [1],
                    [define is you want IBM CBEA build enabled])
            fi

            AC_SUBST(ENABLE_UNIT_TESTS, [$enable_unit_tests])
        ],
        [
            if test -n "$UNIT_TESTS" ; then
                enable_unit_tests=1
            else
                enable_unit_tests=0
            fi

            if test "$enable_unit_tests" = "1" ; then
                AC_DEFINE(UNIT_TESTS, [1],
                    [define is you want distributed memory parallelism enabled])
            fi

            AC_SUBST(ENABLE_UNIT_TESTS, [$enable_unit_tests])
    ])
])
