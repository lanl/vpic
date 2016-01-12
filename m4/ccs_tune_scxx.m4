AC_DEFUN([CCS_TUNE_SCXX], [
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([AC_PROG_CXX])

    cxx=`echo $CXX | sed 's%^.*[[/$]]%%g'`

    case "$cxx" in
        spu-g++)
            CCS_FLAGS_GNU_SCXX
        ;;
    esac
])
