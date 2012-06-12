# should be extended beyond single and double precision
AC_DEFUN([CCS_WITH_PRECISION], [
    AC_ARG_VAR([PRECISION], [Specify precision <single, double>])

    AC_ARG_WITH([precision],
        [AC_HELP_STRING([--with-precision],
        [Specifying addressing <single, double>])],
        [
            AC_MSG_CHECKING(precision)
            if test -n "$PRECISION" ; then
                case "$PRECISION" in
                    single)
			AC_DEFINE(__SINGLEPREC__)
                        AC_MSG_RESULT(single precision)
                    ;;
                    double)
			AC_DEFINE(__DOUBLEPREC__)
                        AC_MSG_RESULT(double precision)
                    ;;
                    *)
                        AC_MSG_ERROR($PRECISION is an invalid precision option)
                    ;;
                esac
            elif test "$withval" != no ; then
                case "$withval" in
                    single)
			AC_DEFINE(__SINGLEPREC__)
                        AC_MSG_RESULT(single precision)
                    ;;
                    double)
			AC_DEFINE(__DOUBLEPREC__)
                        AC_MSG_RESULT(double precision)
                    ;;
                    *)
                        AC_MSG_ERROR($withval is an invalid precision option)
                    ;;
                esac
            else
		AC_DEFINE(__DOUBLEPREC__)
                AC_MSG_RESULT(double precision)
            fi
        ],
        [
            AC_MSG_CHECKING(precision)

    	    # allow usage CCS_WITH_PRECISION(value) from within configure.ac
	    if test "$1" != "" ; then
               PRECISION="$1"
            fi

            if test -n "$PRECISION" ; then
                case "$PRECISION" in
                    single)
			AC_DEFINE(__SINGLEPREC__)
                        AC_MSG_RESULT(single precision)
                    ;;
                    double)
			AC_DEFINE(__DOUBLEPREC__)
                        AC_MSG_RESULT(double precision)
                    ;;
                    *)
                        AC_MSG_ERROR($PRECISION is an invalid precision option)
                    ;;
                esac
            else
		AC_DEFINE(__DOUBLEPREC__)
                AC_MSG_RESULT(double precision)
            fi
    ])
])
