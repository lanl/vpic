dnl ----------------------------------------------------------------------------
dnl CCS_WITH_ISA
dnl ----------------------------------------------------------------------------
AC_DEFUN([CCS_WITH_ISA], [
    AC_ARG_VAR([SAL_ISA], [])

    AC_ARG_WITH([isa],
        [AC_HELP_STRING([--with-isa],
        [Specify instruction set architecture
            <scalar,sse,sse2,sse3,sse4a,sse5,altivec,spu>.
            Specifying 'sse' will configure for the highest level
            of SSE supported by the GNU preprocessor.])],
        [
            if test -n "$SAL_ISA" ; then
                AC_SUBST(SAL_ISA, $SAL_ISA)
            elif test "$withval" != no ; then
                AC_SUBST(SAL_ISA, $withval)
            fi
        ],
        [
            if test -n "$SAL_ISA" ; then
                AC_SUBST(SAL_ISA, $SAL_ISA)
            else
                AC_SUBST(SAL_ISA, "scalar")
            fi
    ])

    AC_MSG_CHECKING(isa for vectorization)

    case "$SAL_ISA" in

        scalar)
            AC_DEFINE(SAL_SCALAR)
        ;;

        sse)
            dnl FIXME: add support for automatic detection of SSE level
            sselist=`(touch test.h; cpp -dM test.h; rm -f test.h)`

            dnl query preprocessor for SSE support
            sse5=`echo $sselist | grep "__SSE5__ 1"`
            sse4=`echo $sselist | grep "__SSE4__ 1"`
            sse3=`echo $sselist | grep "__SSE3__ 1"`
            sse2=`echo $sselist | grep "__SSE2__ 1"`

            if test -n "$sse5" ; then
                AC_DEFINE(SAL_SSE)
                AC_DEFINE(SAL_SSE2)
                AC_DEFINE(SAL_SSE3)
                AC_DEFINE(SAL_SSE5)
            elif test -n "$sse4" ; then
                AC_DEFINE(SAL_SSE)
                AC_DEFINE(SAL_SSE2)
                AC_DEFINE(SAL_SSE3)
                AC_DEFINE(SAL_SSE4)
            elif test -n "$sse3" ; then
                AC_DEFINE(SAL_SSE)
                AC_DEFINE(SAL_SSE2)
                AC_DEFINE(SAL_SSE3)
            elif test -n "$sse2" ; then
                AC_DEFINE(SAL_SSE)
                AC_DEFINE(SAL_SSE2)
            else
                AC_DEFINE(SAL_SSE)
                AC_DEFINE(SAL_SSE2)
            fi
        ;;

        sse2)
            AC_DEFINE(SAL_SSE)
            AC_DEFINE(SAL_SSE2)
        ;;

        sse3)
            AC_DEFINE(SAL_SSE)
            AC_DEFINE(SAL_SSE2)
            AC_DEFINE(SAL_SSE3)
        ;;

        sse4a) dnl Intel only !!!
            AC_DEFINE(SAL_SSE)
            AC_DEFINE(SAL_SSE2)
            AC_DEFINE(SAL_SSE3)
            AC_DEFINE(SAL_SSE4)
        ;;

        sse5) dnl AMD only !!! FIXME: does AMD support SSE4??? 
            AC_DEFINE(SAL_SSE)
            AC_DEFINE(SAL_SSE2)
            AC_DEFINE(SAL_SSE3)
            AC_DEFINE(SAL_SSE5)
        ;;

        spu)
            AC_DEFINE(SAL_SPU)
        ;;

        altivec)
            AC_DEFINE(SAL_ALTIVEC)
        ;;

    esac

    AC_MSG_RESULT($SAL_ISA)
])
