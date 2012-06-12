AC_DEFUN([CCS_TUNE_SCC], [
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([AC_PROG_CC])

    cc=`echo $CC | sed 's%^.*[[/$]]%%g'`

    case "$cc" in
        spu-gcc)
            CCS_FLAGS_GNU_SCC
        ;;
        spuxlc)
            CCS_FLAGS_IBM_SCC
        ;;
    esac
])
