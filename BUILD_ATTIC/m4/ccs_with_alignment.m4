AC_DEFUN([CCS_WITH_ALIGNMENT], [
    AC_ARG_VAR([ALIGNMENT], [Memory byte alignment <Natural, 16, 128>])

    AC_ARG_WITH([alignment],
        [AC_HELP_STRING([--with-alignment],
        [Memory byte alignment <Natural, 16, 128>])],
        [
            AC_MSG_CHECKING(Memory alignment)
            if test -n "$ALIGNMENT" ; then
                case "$ALIGNMENT" in
                    Natural)
                        AC_DEFINE([ALIGNMENT_NATURAL],[1],[Use native alignment])
                        AC_MSG_RESULT(Natural)
                    ;;
                    128)
                        AC_DEFINE([ALIGNMENT_128],[1],[Use 128-byte alignment])
                        AC_MSG_RESULT(128 Byte)
                    ;;
                    16)
                        AC_DEFINE([ALIGNMENT_16],[1],[Use 16-byte alignment])
                        AC_MSG_RESULT(16 Byte)
                    ;;
                    *)
                        AC_MSG_ERROR(Invalid Alignment Option)
                    ;;
                esac
            elif test "$withval" != no ; then
                case "$withval" in
                    Natural)
                        AC_DEFINE([ALIGNMENT_NATURAL],[1],[Use native alignment])
                        AC_MSG_RESULT(Natural)
                    ;;
                    128)
                        AC_DEFINE([ALIGNMENT_128],[1],[Use 128-byte alignment])
                        AC_MSG_RESULT(128 Byte)
                    ;;
                    16)
                        AC_DEFINE([ALIGNMENT_16],[1],[Use 16-byte alignment])
                        AC_MSG_RESULT(16 Byte)
                    ;;
                    *)
                        AC_MSG_ERROR(Invalid Alignment Option)
                    ;;
                esac
            fi
        ],
        [
            AC_MSG_CHECKING(Memory alignment)
            if test -n "$ALIGNMENT" ; then
                case "$ALIGNMENT" in
                    Natural)
                        AC_DEFINE([ALIGNMENT_NATURAL],[1],[Use native alignment])
                        AC_MSG_RESULT(Natural)
                    ;;
                    128)
                        AC_DEFINE([ALIGNMENT_128],[1],[Use 128-byte alignment])
                        AC_MSG_RESULT(128 Byte)
                    ;;
                    16)
                        AC_DEFINE([ALIGNMENT_16],[1],[Use 16-byte alignment])
                        AC_MSG_RESULT(16 Byte)
                    ;;
                    *)
                        AC_MSG_ERROR(Invalid Alignment Option)
                    ;;
                esac
            else
                AC_DEFINE([ALIGNMENT_NATURAL],[1],[Use native alignment])
                AC_MSG_RESULT(Natural)
            fi
    ])
])
