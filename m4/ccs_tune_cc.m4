AC_DEFUN([CCS_TUNE_CC], [
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])
    AC_REQUIRE([AC_PROG_CC])

    cc=`echo $CC | sed 's%^.*[[/$]]%%g'`

    case "$cc" in
        *-gcc)
            CCS_FLAGS_GNU_PCC
        ;;
        gcc*)
            CCS_FLAGS_GNU_CC
        ;;
        icc)
            CCS_FLAGS_INTEL_CC
        ;;
    esac
])
