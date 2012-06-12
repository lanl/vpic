AC_DEFUN([CCS_FLAGS_GNU_SCXX], [
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([CCS_ENABLE_DEBUG_SYMBOLS])
    AC_REQUIRE([AC_PROG_CXX])

    if test "$OPT_LEVEL" != "0" ; then
        AC_MSG_CHECKING(SPU g++ tuning)

        if test "$ENABLE_DEBUG_SYMBOLS" = "1" ; then
            CXXFLAGS="-g $CXXFLAGS"
        fi

        case "$OPT_LEVEL" in
            -2 | -1)
                CXXFLAGS="$CXXFLAGS -g -O0 -fno-inline"
                AC_MSG_RESULT($CXXFLAGS)
            ;;
            1)
                CXXFLAGS="$CXXFLAGS -O1"
                AC_MSG_RESULT($CXXFLAGS)
            ;;
            2)
                CXXFLAGS="$CXXFLAGS -O2"
                AC_MSG_RESULT($CXXFLAGS)
            ;;
            3)
                CXXFLAGS="$CXXFLAGS -O3"
                AC_MSG_RESULT($CXXFLAGS)
            ;;
            *)
                CXXFLAGS="$CXXFLAGS -O3 -funroll-loops"
                AC_MSG_RESULT($CXXFLAGS)
            ;;
        esac
    fi
])
