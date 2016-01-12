AC_DEFUN([CCS_ENABLE_DEBUG_SYMBOLS], [
    #
    # User hints...
    #
    AC_ARG_VAR([DEBUG_SYMBOLS], [Force generation of debug symbols])
    AC_ARG_ENABLE([debug-symbols],
        [AC_HELP_STRING([--enable-debug-symbols],
        [Force generation of debug symbols])],
        [
            if test -n "$DEBUG_SYMBOLS" ; then
                enable_debug_symbols=1
            elif test "$enableval" = "yes" ; then
                enable_debug_symbols=1
            else
                enable_debug_symbols=0
            fi

            if test "$enable_debug_symbols" = "1" ; then
                AC_DEFINE(DEBUG_SYMBOLS, [1],
                    [define is you want debug symbols enabled])
            fi

            AC_SUBST(ENABLE_DEBUG_SYMBOLS, [$enable_debug_symbols])
        ],
        [
            if test -n "$DEBUG_SYMBOLS" ; then
                enable_debug_symbols=1
            else
                enable_debug_symbols=0
            fi

            if test "$enable_debug_symbols" = "1" ; then
                AC_DEFINE(DEBUG_SYMBOLS, [1],
                    [define is you want debug symbols enabled])
            fi

            AC_SUBST(ENABLE_DEBUG_SYMBOLS, [$enable_debug_symbols])
    ])
])
