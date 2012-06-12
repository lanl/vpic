AC_DEFUN([CCS_FLAGS_IBM_SCC], [
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([CCS_ENABLE_DEBUG_SYMBOLS])
    AC_REQUIRE([AC_PROG_CC])

    if test "$OPT_LEVEL" != "0" ; then
        AC_MSG_CHECKING(SPU xlc tuning)

        if test "$ENABLE_DEBUG_SYMBOLS" = "1" ; then
            CFLAGS="-g $CFLAGS"
        fi

        case "$OPT_LEVEL" in
            -2 | -1)
                CFLAGS="$CFLAGS -g -qnoopt -O0 -qnoinline"
                AC_MSG_RESULT($CFLAGS)
            ;;
            1)
                CFLAGS="$CFLAGS -O1"
                AC_MSG_RESULT($CFLAGS)
            ;;
            2)
                CFLAGS="$CFLAGS -O2"
                AC_MSG_RESULT($CFLAGS)
            ;;
            3)
                CFLAGS="$CFLAGS -O3"
                AC_MSG_RESULT($CFLAGS)
            ;;
            *)
                CFLAGS="$CFLAGS -O3 -qunroll"
                AC_MSG_RESULT($CFLAGS)
            ;;
        esac
    fi
])
