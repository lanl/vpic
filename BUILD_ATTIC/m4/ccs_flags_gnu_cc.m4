AC_DEFUN([CCS_FLAGS_GNU_CC], [
    AC_REQUIRE([CCS_WITH_CPUMODEL])
    AC_REQUIRE([CCS_ENABLE_PROFILING])
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([CCS_ENABLE_DEBUG_SYMBOLS])
    AC_REQUIRE([CCS_WITH_ADDRESSING])
    AC_REQUIRE([AC_PROG_CC])

    AC_MSG_CHECKING(gcc tuning)

    if test "$BITS" = "32" ; then
        CFLAGS="-m32 $CFLAGS"
    else
        CFLAGS="-m64 $CFLAGS"
    fi

    CFLAGS="$CFLAGS"

    if test "$OPT_LEVEL" != "0" ; then

        proc_flag=""
        if test "$CPUMODEL" != "generic" ; then
            case "$CPUMODEL" in
                opteron | k8 | athlon64 | athlon-fx)
                    proc_flag="-march=k8"
                ;;
                pentium4)
                    proc_flag="-march=pentium4"
                ;;
                prescott)
                    proc_flag="-march=prescott"
                ;;
                nocona)
                    proc_flag="-march=nocona"
                ;;
                pentium-m)
                    proc_flag="-march=pentium-m"
                ;;
            esac
        fi

        if test "$ENABLE_PROFILING" = "native" ; then
            CFLAGS="$CFLAGS -pg"
        fi

        if test "$ENABLE_DEBUG_SYMBOLS" = "1" ; then
            CFLAGS="-g $CFLAGS"
        fi

        case "$OPT_LEVEL" in
            -2 | -1)
                CFLAGS="$CFLAGS -g -O0 -fno-inline"
            ;;
            1)
                CFLAGS="$CFLAGS -O1"
            ;;
            2)
                CFLAGS="$CFLAGS -O2"
            ;;
            *)
                CFLAGS="$CFLAGS -O3 -funroll-loops $proc_flag"
            ;;
        esac
    fi

    AC_MSG_RESULT($CFLAGS)
])
