AC_DEFUN([CCS_WITH_V4], [
    AC_ARG_VAR([V4], [Short vector ISA selection])

    AC_ARG_WITH([v4],
        [AC_HELP_STRING([--with-v4],
        [Short vector ISA selection])],
        [
            if test -n "$V4" ; then
                v4 = $V4
            elif test "$withval" != no ; then
                v4 = $withval
            fi
        ],
        [
            if test -n "$V4" ; then
                v4 = $V4
            else
                v4 = portable
            fi
    ])

    AC_MSG_CHECKING(v4 selection)
    case $v4 in
        portable)
            AC_DEFINE(USE_V4_PORTABLE)
            AC_MSG_RESULT(portable)
        ;;
        sse)
            AC_DEFINE(USE_V4_SSE)
            AC_MSG_RESULT(SSE)
        ;;
        altivec)
            AC_DEFINE(USE_V4_ALTIVEC)
            AC_MSG_RESULT(Altivec)
        ;;
        spu)
            AC_DEFINE(USE_V4_SPU)
            AC_MSG_RESULT(SPU)
        ;;
    esac
])
