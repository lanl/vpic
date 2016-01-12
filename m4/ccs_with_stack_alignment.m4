dnl ----------------------------------------------------------------------------
dnl CCS_WITH_STACK_ALIGNMENT
dnl
dnl ARG1: stack-alignment file directory ($1)
dnl ----------------------------------------------------------------------------
AC_DEFUN([CCS_WITH_STACK_ALIGNMENT], [
    AC_ARG_VAR([STACK_ALIGNMENT], [Stack alignment in bytes])

    AC_ARG_WITH([stack-alignment],
        [AC_HELP_STRING([--with-stack-alignment],
        [Specify stack-alignment of target system])],
        [
            if test -n "$STACK_ALIGNMENT" ; then
                AC_SUBST(stack_alignment, [$STACK_ALIGNMENT])
            elif test "$withval" != no ; then
                AC_SUBST(stack_alignment, [$withval])
            fi
        ],
        [
            if test -n "$STACK_ALIGNMENT" ; then
                AC_SUBST(stack_alignment, [$STACK_ALIGNMENT])
            else
                AC_SUBST(stack_alignment, [4])
            fi
    ])

	case "$stack_alignment" in
		4)
			AC_DEFINE(STACK_ALIGNMENT, 4)
		;;
		8)
			AC_DEFINE(STACK_ALIGNMENT, 8)
		;;
		16)
			AC_DEFINE(STACK_ALIGNMENT, 16)
		;;
		32)
			AC_DEFINE(STACK_ALIGNMENT, 32)
		;;
		64)
			AC_DEFINE(STACK_ALIGNMENT, 64)
		;;
		128)
			AC_DEFINE(STACK_ALIGNMENT, 128)
		;;
		256)
			AC_DEFINE(STACK_ALIGNMENT, 256)
		;;
		512)
			AC_DEFINE(STACK_ALIGNMENT, 512)
		;;
	esac
])
