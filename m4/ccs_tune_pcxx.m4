AC_DEFUN([CCS_TUNE_PCXX], [
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([AC_PROG_CXX])

    cxx=`echo $CXX | sed 's%^.*[[/$]]%%g'`

    case "$cxx" in
        ppu-g++ | ppu32-g++)
            CCS_FLAGS_GNU_PCXX
        ;;
        ppuxlc++)
            CCS_FLAGS_IBM_PCXX
        ;;
    esac
])
