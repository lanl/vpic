AC_DEFUN([CCS_TUNE_CXX], [
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([AC_PROG_CXX])

    cxx=`echo $CXX | sed 's%^.*[[/$]]%%g'`

    case "$cxx" in
        g++*)
            CCS_FLAGS_GNU_CXX
        ;;
        icpc)
            CCS_FLAGS_INTEL_CXX
        ;;
        pathCC)
            CCS_FLAGS_PATHSCALE_CXX
        ;;
    esac
])
