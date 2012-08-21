AC_DEFUN([OCL_WITH_BUILDSTYLE], [
    AC_ARG_VAR([BUILDSTYLE], [Specify build style])
    AC_ARG_WITH([buildstyle],
        [AC_HELP_STRING([--with-buildstyle],
        [Specify build style])],
        [
            AC_MSG_CHECKING(Build Style)
            if test -n "$BUILDSTYLE" ; then
                OCL_SET_BUILDSTYLE($BUILDSTYLE)
            elif test "$withval" != no ; then
                OCL_SET_BUILDSTYLE($withval)
            else
                OCL_SET_BUILDSTYLE(standard)
            fi
        ],
        [
            if test -n "$BUILDSTYLE" ; then
                OCL_SET_BUILDSTYLE($BUILDSTYLE)
            else
                OCL_SET_BUILDSTYLE(standard)
            fi
    ])
])

AC_DEFUN([OCL_SET_BUILDSTYLE], [
    AC_MSG_CHECKING(Build Style)
    case "$1" in

        standard)
            AC_DEFINE([BUILDSTYLE], [standard],
                [Standard build])
            AC_SUBST(buildstyle, [standard])
            AC_MSG_RESULT([standard])
        ;;

        ocl_apple)
            AC_DEFINE([BUILDSTYLE], [ocl_apple], [Apple OpenCL])
            AC_SUBST(buildstyle, [ocl_apple])
            AC_MSG_RESULT([ocl_apple])
        ;;

        *)
            AC_MSG_ERROR(Invalid Build Style Option)
        ;;

    esac
])
