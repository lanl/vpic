AC_DEFUN([CCS_WITH_ADDRESSING], [
    AC_ARG_VAR([ADDRESSING], [Specify addressing <32, 64>])
    AC_ARG_WITH([addressing],
        [AC_HELP_STRING([--with-addressing],
        [Specifying addressing <32, 64>])],
        [
            AC_MSG_CHECKING(Addressing)
            if test -n "$ADDRESSING" ; then
                case "$ADDRESSING" in
                    32)
                        AC_DEFINE([ADDRESSING_32], [1],
                            [Using 32-bit addressing])
                        AC_SUBST(BITS, [32])
                        AC_MSG_RESULT(32 bit)
                    ;;
                    64)
                        AC_DEFINE([ADDRESSING_64], [1],
                            [Using 64-bit addressing])
                        AC_SUBST(BITS, [64])
                        AC_MSG_RESULT(64 bit)
                    ;;
                    *)
                        AC_MSG_ERROR(Invalid Addressing Option)
                    ;;
                esac
            elif test "$withval" != no ; then
                case "$withval" in
                    32)
                        AC_DEFINE([ADDRESSING_32], [1],
                            [Using 32-bit addressing])
                        AC_SUBST(BITS, [32])
                        AC_MSG_RESULT(32 bit)
                    ;;
                    64)
                        AC_DEFINE([ADDRESSING_64], [1],
                            [Using 64-bit addressing])
                        AC_SUBST(BITS, [64])
                        AC_MSG_RESULT(64 bit)
                    ;;
                    *)
                        AC_MSG_ERROR(Invalid Addressing Option)
                    ;;
                esac
            else
                AC_DEFINE([ADDRESSING_64], [1], [Using 64-bit addressing])
                AC_SUBST(BITS, [64])
                AC_MSG_RESULT(64 bit)
            fi
        ],
        [
            AC_MSG_CHECKING(Addressing)
            if test -n "$ADDRESSING" ; then
                case "$ADDRESSING" in
                    32)
                        AC_DEFINE([ADDRESSING_32], [1],
                            [Using 32-bit addressing])
                        AC_SUBST(BITS, [32])
                        AC_MSG_RESULT(32 bit)
                    ;;
                    64)
                        AC_DEFINE([ADDRESSING_64], [1],
                            [Using 64-bit addressing])
                        AC_SUBST(BITS, [64])
                        AC_MSG_RESULT(64 bit)
                    ;;
                    *)
                        AC_MSG_ERROR(Invalid Addressing Option)
                    ;;
                esac
            else
                AC_DEFINE([ADDRESSING_64], [1], [Using 64-bit addressing])
                AC_SUBST(BITS, [64])
                AC_MSG_RESULT(64 bit)
            fi
    ])
])
