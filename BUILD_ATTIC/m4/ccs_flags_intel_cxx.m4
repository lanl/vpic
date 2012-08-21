AC_DEFUN([CCS_FLAGS_INTEL_CXX], [
    AC_REQUIRE([CCS_WITH_CPUMODEL])
    AC_REQUIRE([CCS_ENABLE_PROFILING])
    AC_REQUIRE([CCS_ENABLE_DEBUG_SYMBOLS])
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([AC_PROG_CXX])

    AC_MSG_CHECKING(Intel icpc tuning)

    if test "$OPT_LEVEL" != "0" ; then

        proc_flag=""
        if test "$CPUMODEL" != "generic" ; then
            case "$CPUMODEL" in
                pentium)
                    proc_flag="-mcpu=pentium"
                ;;
                pentium2)
                    proc_flag="-xi"
                ;;
                pentium3)
                    proc_flag="-xK"
                ;;
                pentium4)
                    proc_flag="-xW"
                ;;
                pentium-m)
                    proc_flag="-xB"
                ;;
                itanium)
                    proc_flag="-mcpu=itanium"
                ;;
                itanium2)
                    proc_flag="-mcpu=itanium2"
                ;;
            esac
        fi

        if test "$ENABLE_PROFILING" = "native" ; then
            CXXFLAGS="$CXXFLAGS -pg"
        fi

        if test "$ENABLE_DEBUG_SYMBOLS" = "1" ; then
            CXXFLAGS="-g $CXXFLAGS"
        fi

        case "$OPT_LEVEL" in
            -2 | -1)
                CXXFLAGS="$CXXFLAGS -g -O0"
            ;;
            1)
                CXXFLAGS="$CXXFLAGS -O1"
            ;;
            2)
                CXXFLAGS="$CXXFLAGS -O2"
            ;;
            *)
    dnl         CXXFLAGS="$CXXFLAGS -O3 -fno-alias $proc_flag -opt_report -opt_report_file $"*".log -opt_report_level max"
                CXXFLAGS="$CXXFLAGS -O3 -fno-alias $proc_flag"
            ;;
        esac
    fi

    AC_MSG_RESULT($CXXFLAGS)
])
