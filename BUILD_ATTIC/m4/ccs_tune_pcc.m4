AC_DEFUN([CCS_TUNE_PCC], [
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([AC_PROG_CC])

    cc=`echo $CC | sed 's%^.*[[/$]]%%g'`

    case "$cc" in
        ppu-gcc | ppu32-gcc)
            CCS_FLAGS_GNU_PCC
        ;;
        ppuxlc)
            CCS_FLAGS_IBM_PCC
        ;;
    esac
])
