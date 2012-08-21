AC_DEFUN([CCS_ENABLE_CELL], [
    #
    # User hints...
    #
    AC_ARG_VAR([ENABLE_CELL], [Enable IBM CBEA build])
    AC_ARG_ENABLE([cell],
        [AC_HELP_STRING([--enable-cell],
        [enable IBM CBEA build])],
        [
            if test -n "$ENABLE_CELL" ; then
                enable_cell=1
            elif test "$enableval" = "yes" ; then
                enable_cell=1
            else
                enable_cell=0
            fi
        ],
        [
            if test -n "$ENABLE_CELL" ; then
                enable_cell=1
            else
                enable_cell=0
            fi
        ])

        AC_SUBST(ENABLE_CELL, [$enable_cell])
])
